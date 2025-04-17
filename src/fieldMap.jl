"""
    EMField

Represents the electromagnetic field with electric (E) and magnetic (B) components.
"""
struct EMField{T}
    E::SVector{3, T}
    B::SVector{3, T}
end

"""
    FieldMap{T}

Holds the raw data and coordinate grids used for building interpolation objects
that return E, B. You can store these as arrays or any structure needed
by your interpolation approach.
"""
struct FieldMap{T}
    Edata::Union{Nothing, Array{T,4}}   # e.g. E[i,j,k,component]
    Bdata::Union{Nothing, Array{T,4}}   # e.g. B[i,j,k,component]
    metadata::Dict{String, Any}         # Metadata dictionary
    # Field interpolators
    Eitp::Union{Nothing,NTuple{3,Interpolations.GriddedInterpolation}}
    Bitp::Union{Nothing,NTuple{3,Interpolations.GriddedInterpolation}}
end

"""
    load_fieldmap(path::String)::FieldMap

Load a FieldMap from an HDF5 file following the format used in Bmad.
The function reads electric and magnetic field data along with metadata.
"""
function load_fieldmap(path::String)
    h5open(path, "r") do file
        # Read metadata
        basePath = ""
        if haskey(attributes(file), "externalFieldPath")
            basePath = read_attribute(file, "externalFieldPath")
        end
        
        # Handle the case with or without %T in the path
        root_group = file
        if occursin("/%T/", basePath)
            # For multiple grids, we just take the first one for now
            base_dir = basePath[1:findfirst("/%T/", basePath)[1]-1]
            grid_groups = keys(file[base_dir])
            # Find the first valid grid group (should be an integer)
            grid_id = ""
            for name in grid_groups
                if tryparse(Int, name) !== nothing
                    grid_id = name
                    break
                end
            end
            if grid_id == ""
                error("No valid grid found in HDF5 file")
            end
            root_group = file[join([base_dir, grid_id],"/")]
        elseif basePath != ""
            root_group = file[basePath]
        end
        
        # Extract metadata attributes into a dictionary
        metadata = Dict(
            "lower_bound" => haskey(attributes(root_group), "gridLowerBound") ? 
                             read_attribute(root_group, "gridLowerBound") : [1, 1, 1],
            "origin_offset" => haskey(attributes(root_group), "gridOriginOffset") ? 
                               read_attribute(root_group, "gridOriginOffset") : [0.0, 0.0, 0.0],
            "spacing" => haskey(attributes(root_group), "gridSpacing") ? 
                         read_attribute(root_group, "gridSpacing") : [1.0, 1.0, 1.0],
            "grid_size" => haskey(attributes(root_group), "gridSize") ? 
                           read_attribute(root_group, "gridSize") : [1, 1, 1],
            "geometry" => haskey(attributes(root_group), "gridGeometry") ? 
                          read_attribute(root_group, "gridGeometry") : "rectangular",
            "field_scale" => haskey(attributes(root_group), "fieldScale") ? 
                             only(read_attribute(root_group, "fieldScale")) : 1.0,
            "harmonic" => haskey(attributes(root_group), "harmonic") ? 
                          only(read_attribute(root_group, "harmonic")) : 0,
            "phi0_fieldmap" => haskey(attributes(root_group), "RFphase") ? 
                               only(read_attribute(root_group, "RFphase")) : 0.0,
            "ele_anchor_pt" => haskey(attributes(root_group), "eleAnchorPt") ? 
                               read_attribute(root_group, "eleAnchorPt") : "beginning",
            "interpolation_order" => haskey(attributes(root_group), "interpolationOrder") ? 
                                     only(read_attribute(root_group, "interpolationOrder")) : 1,
            "fundamental_frequency" => haskey(attributes(root_group), "fundamentalFrequency") ? 
                                       only(read_attribute(root_group, "fundamentalFrequency")) : 0.0,
            "axisLabels" => haskey(attributes(root_group), "axisLabels") ? 
                            read_attribute(root_group, "axisLabels") : ["x", "y", "z"]
        )

        # Use metadata dictionary
        origin_offset = metadata["origin_offset"]
        spacing = metadata["spacing"]
        grid_size = metadata["grid_size"]
        field_scale = metadata["field_scale"]
        geometry = metadata["geometry"]
        axis_labels = metadata["axisLabels"]
        
        # Create coordinate grids
        xgrid = [origin_offset[1] + i * spacing[1] for i in 0:(grid_size[1] - 1)]
        ygrid = [origin_offset[2] + i * spacing[2] for i in 0:(grid_size[2] - 1)]
        zgrid = [origin_offset[3] + i * spacing[3] for i in 0:(grid_size[3] - 1)]
        
        # Initialize field arrays
        Edata = zeros(ComplexF64, grid_size[1], grid_size[2], grid_size[3], 3)
        Bdata = zeros(ComplexF64, grid_size[1], grid_size[2], grid_size[3], 3)
        
        # Define component names based on geometry
        logical_labels = geometry == "cylindrical" ? ["r", "theta", "z"] : ["x", "y", "z"]
        
        # Determine data ordering (C or Fortran)
        data_order = ""
        if haskey(attributes(root_group), "axisLabels")
            if all(logical_labels .== axis_labels)
                data_order = "C"
            elseif all(logical_labels .== reverse(axis_labels))
                data_order = "F"
            end
        end
        
        # Function to read field component and handle pseudo-datasets
        function read_field_component(group, component_name)
            if !haskey(group, component_name)
                return zeros(ComplexF64, grid_size...)
            end
            
            component = group[component_name]
            
            # Check if this is a pseudo-dataset (constant value)
            if typeof(component) <: HDF5.Group && haskey(attributes(component), "value")
                value = read_attribute(component, "value")
                if !isa(value, Complex)
                    value = complex(value)
                end
                return fill(value, grid_size...)
            end
            
            # Regular dataset
            return read(component)
        end
        
        # Read electric field data
        if haskey(root_group, "electricField")
            e_group = root_group["electricField"]
            for (i, label) in enumerate(logical_labels)
                if haskey(e_group, label)
                    data = read_field_component(e_group, label)
                    
                    # Apply data ordering
                    if data_order == "C"
                        for ix in 1:grid_size[1], iy in 1:grid_size[2], iz in 1:grid_size[3]
                            Edata[ix, iy, iz, i] = data[ix, iy, iz] * field_scale
                        end
                    elseif data_order == "F"
                        for ix in 1:grid_size[1], iy in 1:grid_size[2], iz in 1:grid_size[3]
                            Edata[ix, iy, iz, i] = data[iz, iy, ix] * field_scale
                        end
                    else
                        # Default to C ordering if not specified
                        for ix in 1:grid_size[1], iy in 1:grid_size[2], iz in 1:grid_size[3]
                            Edata[ix, iy, iz, i] = data[ix, iy, iz] * field_scale
                        end
                    end
                end
            end
        end
        
        # Read magnetic field data
        if haskey(root_group, "magneticField")
            b_group = root_group["magneticField"]
            for (i, label) in enumerate(logical_labels)
                if haskey(b_group, label)
                    data = read_field_component(b_group, label)
                    
                    # Apply data ordering
                    if data_order == "C"
                        for ix in 1:grid_size[1], iy in 1:grid_size[2], iz in 1:grid_size[3]
                            Bdata[ix, iy, iz, i] = data[ix, iy, iz] * field_scale
                        end
                    elseif data_order == "F"
                        for ix in 1:grid_size[1], iy in 1:grid_size[2], iz in 1:grid_size[3]
                            Bdata[ix, iy, iz, i] = data[iz, iy, ix] * field_scale
                        end
                    else
                        # Default to C ordering if not specified
                        for ix in 1:grid_size[1], iy in 1:grid_size[2], iz in 1:grid_size[3]
                            Bdata[ix, iy, iz, i] = data[ix, iy, iz] * field_scale
                        end
                    end
                end
            end
        end
        
        return Edata, Bdata, metadata
    end
end

"""
    _make_itps(knots, A)

Return a tuple `(itp₁,itp₂,itp₃)` of gridded interpolators
for the three vector components stored in the 4‑D field array `A`.
Singleton axes are treated as constants (NoInterp), other axes are linear.
"""
function _make_itps(knots::NTuple{3,Vector},
                    A::Array{<:Number,4})

    # choose a scheme for every spatial axis
    schemes = ntuple(d -> size(A,d)==1 ? NoInterp() : Gridded(Linear()), 3)

    # build one interpolator per vector component
    ntuple(c -> begin
        base_itp = interpolate(knots,
                               view(A,:,:,:,c),
                               Gridded(schemes))
    end, 3)
end

"""
    _interpolants(knots, E, B)

Generates interpolators for the electric (E) and magnetic (B) field components
using the provided coordinate grids (`knots`). Returns a tuple of interpolators
for E and B, or `nothing` if the respective field data is not provided.
"""
function _interpolants(knots, E, B)
    Eitp = E === nothing ? nothing : _make_itps(knots, E)
    Bitp = B === nothing ? nothing : _make_itps(knots, B)
    return Eitp, Bitp
end

"""
    FieldMap(path::String)::FieldMap

Constructs a `FieldMap` object by loading field data and metadata from the
specified HDF5 file. The function reads electric and magnetic field data,
coordinate grids, and other metadata required for interpolation.
"""
function FieldMap(path::String)
    Edata, Bdata, metadata = load_fieldmap(path)
    
    knots = (collect(metadata["lower_bound"][1]:metadata["lower_bound"][1]+metadata["grid_size"][1]-1), 
             collect(metadata["lower_bound"][2]:metadata["lower_bound"][2]+metadata["grid_size"][2]-1), 
             collect(metadata["lower_bound"][3]:metadata["lower_bound"][3]+metadata["grid_size"][3]-1))
    
    Eitp, Bitp = _interpolants(knots, Edata, Bdata)
    
    return FieldMap(Edata, Bdata, metadata, Eitp, Bitp)
end

"""
    FieldMap(x::Vector{T}, y::Vector{T}, z::Vector{T};
             E::Union{Nothing,Array{T,4}} = nothing,
             B::Union{Nothing,Array{T,4}} = nothing,
             geometry::String = "rectangular") where T

Constructs a `FieldMap` object using the provided coordinate grids (`x`, `y`, `z`)
and optional electric (E) and magnetic (B) field data. The geometry of the grid
can be specified as "rectangular" or "cylindrical".
"""
function FieldMap(x::Vector{T}, y::Vector{T}, z::Vector{T};
                  E::Union{Nothing,Array{T,4}} = nothing,
                  B::Union{Nothing,Array{T,4}} = nothing,
                  geometry::String     = "rectangular") where T
    knots = (x, y, z)
    Eitp, Bitp = _interpolants(knots, E, B)
    metadata = Dict{String,Any}(
        "geometry" => geometry
    )
    return FieldMap(E, B, metadata, Eitp, Bitp)
end


"""
    get_fields_index(fm::FieldMap, i::Float64, j::Float64, k::Float64)::EMField

Returns an `EMField` object representing the electric and magnetic field
components at the specified grid indices `(i, j, k)`. The indices correspond
to Cartesian or cylindrical coordinates based on the grid geometry.
"""
function get_fields_index(fm::FieldMap, i::Float64, j::Float64, k::Float64)
    # Small helper closes over the point and handles the “maybe‑tuple”.
    field_at(itp::NTuple{3}, i, j, k) = SVector{3}(map(t -> t(i, j, k), itp)...)
    field_at(::Nothing,      i, j, k) = SVector{3, Float64}(0.0, 0.0, 0.0)    # fallback

    E = field_at(fm.Eitp, i, j, k)
    B = field_at(fm.Bitp, i, j, k)

    return EMField(E, B)
end

"""
    get_fields(fm::FieldMap, x::Float64, y::Float64, z::Float64)::EMField

Returns an `EMField` object representing the electric and magnetic field
components at the specified Cartesian coordinates `(x, y, z)`. If the grid
geometry is cylindrical, the coordinates are converted to cylindrical before
interpolation.
"""
function get_fields(fm::FieldMap, x::Float64, y::Float64, z::Float64)
    if lowercase(fm.metadata["geometry"]) == "rectangular"
        # Already in Cartesian coordinates
        return get_fields_index(fm, x, y, z)
    elseif lowercase(fm.metadata["geometry"]) == "cylindrical"
        # Convert Cartesian (x,y,z) to cylindrical (r,θ,z)
        r = sqrt(x^2 + y^2)
        θ = atan(y, x)  # atan2 equivalent
        
        EM = get_fields_index(fm, r, θ, z)
        
        # Convert field components from cylindrical to Cartesian
        cos_θ, sin_θ = cos(θ), sin(θ)
        
        Ex = EM.E[1] * cos_θ - EM.E[2] * sin_θ
        Ey = EM.E[1] * sin_θ + EM.E[2] * cos_θ
        Bx = EM.B[1] * cos_θ - EM.B[2] * sin_θ
        By = EM.B[1] * sin_θ + EM.B[2] * cos_θ
        
        return EMField(SVector{3}(Ex, Ey, EM.E[3]), SVector{3}(Bx, By, EM.B[3]))
    else
        error("Unsupported geometry: $(fm.metadata["geometry"])")
    end
end
