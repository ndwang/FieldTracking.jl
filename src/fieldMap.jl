"""
    EMField{T}

Represents the electromagnetic field with electric (E) and magnetic (B) components.
"""
struct EMField{T}
    E::SVector{3, T}
    B::SVector{3, T}
end

"""
    FieldMap{E,B}

Holds metadata and interpolation objects for electric (`E`)
and magnetic (`B`) fields. The raw data is stored in interpolation objects.
"""
struct FieldMap{E,B}
    # Required Attributes from openPMD-beamphysics standard
    eleAnchorPt::String
    gridGeometry::String
    gridSpacing::Vector{Float64} # Grid spacing along each axis
    gridLowerBound::Vector{Int}  # Lower bound index for each grid dimension
    gridSize::Vector{Int}        # Number of grid points along each axis
    gridOriginOffset::Vector{Float64} # Physical coordinates of the grid origin (index [1,1,1] or equivalent)
    harmonic::Int                # Harmonic number (0 for static fields)

    # Optional Metadata
    metadata::Dict{String, Any}         # Dictionary for optional metadata

    # Field interpolators
    Eitp::E     # Interpolator for Electric field (e.g., NTuple{3, Interpolations.GriddedInterpolation} or Nothing)
    Bitp::B     # Interpolator for Magnetic field (e.g., NTuple{3, Interpolations.GriddedInterpolation} or Nothing)
end

"""
    load_fieldmap(path::String)

Load field data and attributes from an HDF5 file specified by `path`.
Returns Edata, Bdata, required attributes, and optional metadata.
"""
function load_fieldmap(path::String)
    h5open(path, "r") do file
        # Determine the base path within the HDF5 file for field data
        basePath = ""
        if haskey(attributes(file), "externalFieldPath")
            basePath = read_attribute(file, "externalFieldPath")
        end

        # Handle paths potentially containing '%T'
        root_group = file
        if occursin("/%T/", basePath)
            # For multiple grids, select the first one found
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
            root_group = file[joinpath(base_dir, grid_id)] # Use joinpath for robustness
        elseif basePath != ""
            root_group = file[basePath]
        end

        # Read Required Attributes (provide defaults for robustness, though the spec requires them)
        eleAnchorPt = haskey(attributes(root_group), "eleAnchorPt") ?
                      read_attribute(root_group, "eleAnchorPt") : "beginning"
        gridGeometry = haskey(attributes(root_group), "gridGeometry") ?
                       read_attribute(root_group, "gridGeometry") : "rectangular"
        gridSpacing = haskey(attributes(root_group), "gridSpacing") ?
                      convert(Vector{Float64}, read_attribute(root_group, "gridSpacing")) : [1.0, 1.0, 1.0]
        gridLowerBound = haskey(attributes(root_group), "gridLowerBound") ?
                         convert(Vector{Int}, read_attribute(root_group, "gridLowerBound")) : [1, 1, 1]
        gridSize = haskey(attributes(root_group), "gridSize") ?
                   convert(Vector{Int}, read_attribute(root_group, "gridSize")) : [1, 1, 1]
        gridOriginOffset = haskey(attributes(root_group), "gridOriginOffset") ?
                           convert(Vector{Float64}, read_attribute(root_group, "gridOriginOffset")) : [0.0, 0.0, 0.0]
        harmonic = haskey(attributes(root_group), "harmonic") ?
                   convert(Int, only(read_attribute(root_group, "harmonic"))) : 0

        # Extract optional metadata attributes into a dictionary
        metadata = Dict{String, Any}()
        optional_attrs = [
            "fieldScale", "fundamentalFrequency", "gridCurvatureRadius",
            "name", "RFphase", "interpolationOrder", "axisLabels"
        ]
        # Default values for optional attributes if not present in the file
        defaults = Dict(
            "fieldScale" => 1.0, "fundamentalFrequency" => 0.0, "gridCurvatureRadius" => 0.0,
            "name" => "", "RFphase" => 0.0, "interpolationOrder" => 1, # Default to linear interpolation
            "axisLabels" => ["x", "y", "z"] # Default axis labels
        )

        for attr in optional_attrs
            if haskey(attributes(root_group), attr)
                val = read_attribute(root_group, attr)
                # Ensure scalar attributes read from HDF5 (which might be 1-element arrays) are stored as scalars
                if attr in ["fieldScale", "fundamentalFrequency", "gridCurvatureRadius", "RFphase", "interpolationOrder", "harmonic"]
                    metadata[attr] = only(val)
                else
                    metadata[attr] = val
                end
            end
        end

        # Use required attributes and optional metadata where needed
        field_scale = metadata["fieldScale"] # Apply scaling factor during data reading
        geometry = gridGeometry # Use the required attribute directly

        # Initialize field data arrays (using ComplexF64 as per standard)
        Edata = zeros(ComplexF64, gridSize[1], gridSize[2], gridSize[3], 3)
        Bdata = zeros(ComplexF64, gridSize[1], gridSize[2], gridSize[3], 3)

        # Define logical component names based on grid geometry
        logical_labels = geometry == "cylindrical" ? ["r", "theta", "z"] : ["x", "y", "z"]

        # Determine data ordering (C or Fortran)
        data_order = ""
        if haskey(metadata, "axisLabels") # Check metadata now
            labels_in_meta = metadata["axisLabels"]
            if all(logical_labels .== labels_in_meta)
                data_order = "C"
            elseif all(logical_labels .== reverse(labels_in_meta))
                data_order = "F"
            end
        end

        # Helper function to read a field component dataset or handle pseudo-datasets (constant value)
        function read_field_component(group, component_name)
            if !haskey(group, component_name)
                # Return zeros if the component dataset doesn't exist
                return zeros(ComplexF64, gridSize...)
            end

            component = group[component_name]

            # Check for pseudo-dataset (represented by an HDF5 group with a 'value' attribute)
            if typeof(component) <: HDF5.Group && haskey(attributes(component), "value")
                value = read_attribute(component, "value")
                # Ensure the value is complex
                if !isa(value, Complex)
                    value = complex(Float64(value))
                end
                # Return an array filled with the constant value. 
                # TODO: store a constant value instead of a full array
                return fill(convert(ComplexF64, value), gridSize...)
            end

            return read(component)
        end

        # Read electric field data
        if haskey(root_group, "electricField")
            e_group = root_group["electricField"]
            for (i, label) in enumerate(logical_labels)
                if haskey(e_group, label)
                    data = read_field_component(e_group, label)

                    # Apply data ordering
                    if data_order == "F"
                        for ix in 1:gridSize[1], iy in 1:gridSize[2], iz in 1:gridSize[3]
                            Edata[ix, iy, iz, i] = data[iz, iy, ix] * field_scale
                        end
                    else
                        for ix in 1:gridSize[1], iy in 1:gridSize[2], iz in 1:gridSize[3]
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
                    if data_order == "F"
                        for ix in 1:gridSize[1], iy in 1:gridSize[2], iz in 1:gridSize[3]
                            Bdata[ix, iy, iz, i] = data[iz, iy, ix] * field_scale
                        end
                    else
                        for ix in 1:gridSize[1], iy in 1:gridSize[2], iz in 1:gridSize[3]
                            Bdata[ix, iy, iz, i] = data[ix, iy, iz] * field_scale
                        end
                    end
                                end
            end
                end

        # Return loaded data, required attributes, and the metadata dictionary
        return Edata, Bdata, eleAnchorPt, gridGeometry, gridSpacing, gridLowerBound, gridSize, gridOriginOffset, harmonic, metadata
    end
end

"""
    _make_itps(knots, A)

Return a tuple `(itp₁, itp₂, itp₃)` of gridded interpolators for the three vector
components stored in the N+1 dimensional field array `A`, using the N-dimensional
`knots` (tuple of coordinate vectors for each dimension).
Uses linear interpolation by default.
TODO: Add support for other interpolation methods.
"""
function _make_itps(knots::NTuple{N, AbstractVector},
                    A::AbstractArray{<:Number}) where N

    ndims(A) == N + 1 || error("Dimension mismatch: knots imply $(N)D grid, but data has $(ndims(A)-1) spatial dimensions.")
    size(A)[1:N] == length.(knots) || error("Size mismatch between knots and data array dimensions.")
    size(A, N + 1) == 3 || error("Data array must have size 3 in the last dimension (components).")

    # build one interpolator per vector component
    ntuple(c ->
            interpolate(knots, view(A, ntuple(_ -> Colon(), N)..., c), Gridded(Linear()))
           , 3)
end

"""
    _interpolants(knots, E, B)

Generate interpolators for the electric (`E`) and magnetic (`B`) field component arrays
using the provided coordinate grid vectors (`knots`).
Returns a tuple `(Eitp, Bitp)`, where each element is either a tuple of
interpolators (one per vector component) or `nothing` if the corresponding
field data array (`E` or `B`) is `nothing`.
"""
function _interpolants(knots, E, B)
    # Create interpolators only if the data array is not nothing
    Eitp = E === nothing ? nothing : _make_itps(knots, E)
    Bitp = B === nothing ? nothing : _make_itps(knots, B)
    return Eitp, Bitp
end

"""
    FieldMap(path::String)::FieldMap

Construct a `FieldMap` object by loading field data and metadata from the
HDF5 file specified by `path`. It reads the data, determines the grid structure,
and creates the necessary interpolators.
"""
function FieldMap(path::String)
    # Load raw data and attributes from the HDF5 file
    Edata, Bdata, eleAnchorPt, gridGeometry, gridSpacing, gridLowerBound, gridSize, gridOriginOffset, harmonic, metadata = load_fieldmap(path)

    # Create coordinate vectors (knots) and interpolators based on grid geometry
    if lowercase(gridGeometry) == "cylindrical"
        # Cylindrical geometry: Interpolate on the r-z plane (assuming azimuthal symmetry)
        knots = (gridOriginOffset[1] .+ (gridLowerBound[1]:gridLowerBound[1]+gridSize[1]-1) .* gridSpacing[1],
                gridOriginOffset[3] .+ (gridLowerBound[3]:gridLowerBound[3]+gridSize[3]-1) .* gridSpacing[3])

        # Extract the r-z plane data
        E_slice = Edata === nothing ? nothing : view(Edata, :, 1, :, :)
        B_slice = Bdata === nothing ? nothing : view(Bdata, :, 1, :, :)

        # Create 2D interpolators for the r-z plane
        Eitp, Bitp = _interpolants(knots, E_slice, B_slice)

    elseif lowercase(gridGeometry) == "rectangular"
        # Rectangular geometry: Interpolate in 3D Cartesian space
        knots = (gridOriginOffset[1] .+ (gridLowerBound[1]:gridLowerBound[1]+gridSize[1]-1) .* gridSpacing[1],
                gridOriginOffset[2] .+ (gridLowerBound[2]:gridLowerBound[2]+gridSize[2]-1) .* gridSpacing[2],
                gridOriginOffset[3] .+ (gridLowerBound[3]:gridLowerBound[3]+gridSize[3]-1) .* gridSpacing[3])
        # Create 3D interpolators for the full grid
        Eitp, Bitp = _interpolants(knots, Edata, Bdata)
    else
         error("Unsupported grid geometry: $gridGeometry")
    end

    # Construct and return the FieldMap object
    return FieldMap(eleAnchorPt, gridGeometry, gridSpacing, gridLowerBound, gridSize, gridOriginOffset, harmonic,
                    metadata, Eitp, Bitp)
end

"""
    FieldMap(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64};
             eleAnchorPt::String = "beginning", gridGeometry::String = "rectangular",
             gridLowerBound::Vector{Int} = [1, 1, 1], gridOriginOffset::Vector{Float64} = [0.0, 0.0, 0.0],
             harmonic::Int = 0,
             E::Union{Nothing,Array{T,4}} = nothing,
             B::Union{Nothing,Array{T,4}} = nothing,
             metadata::Dict{String, Any} = Dict{String,Any}()) where T

Construct a `FieldMap` object directly from coordinate vectors and field data arrays.
Assumes a rectangular grid geometry by default. Calculates `gridSize` and `gridSpacing`
from the provided coordinate vectors.
"""
function FieldMap(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64};
                  eleAnchorPt::String = "beginning",
                  gridGeometry::String = "rectangular",
                  gridLowerBound::Vector{Int} = [1, 1, 1],
                  gridOriginOffset::Vector{Float64} = [0.0, 0.0, 0.0],
                  harmonic::Int = 0,
                  E::Union{Nothing,Array{T,4}} = nothing,
                  B::Union{Nothing,Array{T,4}} = nothing,
                  metadata::Dict{String, Any} = Dict{String,Any}()) where T

    # Calculate default gridSize and gridSpacing
    gridSize = [length(x), length(y), length(z)]
    gridSpacing = [x[2] - x[1], y[2] - y[1], z[2] - z[1]]

    # Create interpolators
    knots = (x, y, z)
    Eitp, Bitp = _interpolants(knots, E, B)

    # Construct and return the FieldMap object
    return FieldMap(eleAnchorPt, gridGeometry, gridSpacing, gridLowerBound, gridSize, gridOriginOffset, harmonic,
                    metadata, Eitp, Bitp)
end

"""
    get_fields(fm::FieldMap, x::Float64, y::Float64, z::Float64)::EMField

Return an `EMField` object containing the interpolated electric and magnetic field
vectors (`E`, `B`) at the specified Cartesian coordinates `(x, y, z)`.
Handles coordinate transformations for cylindrical grids.
"""
function get_fields(fm::FieldMap, x::Float64, y::Float64, z::Float64)
    # Helper function to evaluate field components using the interpolator tuple
     function field_at(itp_tuple::Union{Nothing, NTuple{N, Interpolations.AbstractInterpolation}}, coords...) where N
        # Return zero vector if no interpolator exists for this field (E or B)
        itp_tuple === nothing && return SVector{3, ComplexF64}(0.0, 0.0, 0.0)
        # TODO: Ensure coordinates are within bounds or handle extrapolation
        # Rely on Interpolations.jl default behavior for now
        return SVector{3}(map(itp -> itp(coords...), itp_tuple)...)
    end

    if lowercase(fm.gridGeometry) == "rectangular"
        # Rectangular grid: Interpolate directly using Cartesian coordinates (x, y, z)
        E = field_at(fm.Eitp, x, y, z)
        B = field_at(fm.Bitp, x, y, z)
        return EMField(E, B)

    elseif lowercase(fm.gridGeometry) == "cylindrical"
        # Cylindrical grid:
        # 1. Convert input Cartesian coordinates (x, y, z) to cylindrical coordinates (r, θ, z).
        r = sqrt(x^2 + y^2)
        θ = atan(y, x)

        # 2. Interpolate using the physical (r, z) coordinates on the 2D grid.
        # The result gives field components in the *cylindrical basis* (Er, Eθ, Ez) / (Br, Bθ, Bz) at the requested (r, z).
        EM_cyl = EMField(
            field_at(fm.Eitp, r, z), # Interpolate using r and z
            field_at(fm.Bitp, r, z)
        )

        # 3. Convert the interpolated field vectors from the cylindrical basis back to the Cartesian basis.
        cos_θ, sin_θ = cos(θ), sin(θ)
        
        Ex = EM_cyl.E[1] * cos_θ - EM_cyl.E[2] * sin_θ
        Ey = EM_cyl.E[1] * sin_θ + EM_cyl.E[2] * cos_θ
        Ez = EM_cyl.E[3]

        Bx = EM_cyl.B[1] * cos_θ - EM_cyl.B[2] * sin_θ
        By = EM_cyl.B[1] * sin_θ + EM_cyl.B[2] * cos_θ
        Bz = EM_cyl.B[3] # Bz component remains the same

        # Construct the final EMField object with Cartesian components.
        E = SVector(Ex, Ey, Ez)
        B = SVector(Bx, By, Bz)

        return EMField(E, B)
    else
        error("Unsupported grid geometry for field evaluation: $(fm.gridGeometry)")
    end
end
