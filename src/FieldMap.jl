"""
    FieldMap

Holds the raw data and coordinate grids used for building interpolation objects
that return E, B. You can store these as arrays or any structure needed
by your interpolation approach.
"""
struct FieldMap
    xgrid::Vector{Float64}
    ygrid::Vector{Float64}
    zgrid::Vector{Float64}
    Edata::Array{ComplexF64,4}   # e.g. E[i,j,k,component]
    Bdata::Array{ComplexF64,4}   # e.g. B[i,j,k,component]
    # Metadata
    geometry::String          # "rectangular" or "cylindrical"
    field_scale::Float64      # Scale factor for field values
    harmonic::Int             # RF harmonic number
    phi0_fieldmap::Float64    # RF phase
    master_parameter::String  # Master parameter name
    ele_anchor_pt::String     # Element anchor point
    curved_ref_frame::Bool    # Whether using curved reference frame
    interpolation_order::Int  # Interpolation order
end

"""
    load_fieldmap(path::String)::FieldMap

Load a FieldMap from an HDF5 file following the format used in Bmad.
The function reads electric and magnetic field data along with metadata.
"""
function load_fieldmap(path::String)::FieldMap
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
        
        # Read grid geometry and metadata
        geometry = haskey(attributes(root_group), "gridGeometry") ? 
                   read_attribute(root_group, "gridGeometry") : "rectangular"
        
        field_scale = 1.0
        if haskey(attributes(root_group), "fieldScale")
            field_scale = only(read_attribute(root_group, "fieldScale"))
        elseif haskey(attributes(root_group), "componentFieldScale")
            field_scale = only(read_attribute(root_group, "componentFieldScale"))
        end
        
        harmonic = haskey(attributes(root_group), "harmonic") ? 
                   only(read_attribute(root_group, "harmonic")) : 0
        
        phi0 = haskey(attributes(root_group), "RFphase") ? 
               only(read_attribute(root_group, "RFphase")) : 0.0
        
        # Handle RF phase conversion as in Fortran code
        if harmonic != 0
            if haskey(attributes(root_group), "masterParameter") && 
               lowercase(read_attribute(root_group, "masterParameter")) == "lcavity"
                phi0 = phi0 / harmonic
            else
                phi0 = 0.25 - phi0 / harmonic
            end
        end
        
        master_param = ""
        if haskey(attributes(root_group), "masterParameter")
            master_param = read_attribute(root_group, "masterParameter")
        end
        
        # Read additional metadata
        ele_anchor_pt = haskey(attributes(root_group), "eleAnchorPt") ? 
                        read_attribute(root_group, "eleAnchorPt") : "beginning"
        
        interpolation_order = haskey(attributes(root_group), "interpolationOrder") ? 
                             only(read_attribute(root_group, "interpolationOrder")) : 1
        
        # Check for curved reference frame
        curved_ref_frame = false
        if haskey(attributes(root_group), "gridCurvatureRadius")
            rho = read_attribute(root_group, "gridCurvatureRadius")
            curved_ref_frame = rho != 0
        elseif haskey(attributes(root_group), "curvedRefFrame")
            rho = read_attribute(root_group, "curvedRefFrame")
            curved_ref_frame = rho != 0
        end
        
        # Read grid information
        lb = read_attribute(root_group, "gridLowerBound")
        grid_size = read_attribute(root_group, "gridSize")
        origin = read_attribute(root_group, "gridOriginOffset")
        spacing = read_attribute(root_group, "gridSpacing")
        
        # Create coordinate grids
        xgrid = [origin[1] + (i-lb[1])*spacing[1] for i in lb[1]:(lb[1]+grid_size[1]-1)]
        ygrid = [origin[2] + (i-lb[2])*spacing[2] for i in lb[2]:(lb[2]+grid_size[2]-1)]
        zgrid = [origin[3] + (i-lb[3])*spacing[3] for i in lb[3]:(lb[3]+grid_size[3]-1)]
        
        # Initialize field arrays
        Edata = zeros(ComplexF64, grid_size[1], grid_size[2], grid_size[3], 3)
        Bdata = zeros(ComplexF64, grid_size[1], grid_size[2], grid_size[3], 3)
        
        # Determine axis labels and field component names
        axis_labels = ["x", "y", "z"]
        if haskey(attributes(root_group), "axisLabels")
            axis_labels = read_attribute(root_group, "axisLabels")
        end
        
        # Define component names based on geometry
        logical_labels = geometry == "cylindrical" ? ["r", "theta", "z"] : ["x", "y", "z"]
        E_names = geometry == "cylindrical" ? ["Er", "Etheta", "Ez"] : ["Ex", "Ey", "Ez"]
        B_names = geometry == "cylindrical" ? ["Br", "Btheta", "Bz"] : ["Bx", "By", "Bz"]
        
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
        
        # Check for old-style electric field components
        for (i, name) in enumerate(E_names)
            if haskey(root_group, name)
                if typeof(root_group[name]) <: HDF5.Group && haskey(attributes(root_group[name]), "value")
                    # Pseudo-dataset
                    value = read_attribute(root_group[name], "value")
                    if !isa(value, Complex)
                        value = complex(value)
                    end
                    for ix in 1:grid_size[1], iy in 1:grid_size[2], iz in 1:grid_size[3]
                        Edata[ix, iy, iz, i] = value * field_scale
                    end
                else
                    # Regular dataset
                    data = read(root_group[name])
                    
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
        
        # Check for old-style magnetic field components
        for (i, name) in enumerate(B_names)
            if haskey(root_group, name)
                if typeof(root_group[name]) <: HDF5.Group && haskey(attributes(root_group[name]), "value")
                    # Pseudo-dataset
                    value = read_attribute(root_group[name], "value")
                    if !isa(value, Complex)
                        value = complex(value)
                    end
                    for ix in 1:grid_size[1], iy in 1:grid_size[2], iz in 1:grid_size[3]
                        Bdata[ix, iy, iz, i] = value * field_scale
                    end
                else
                    # Regular dataset
                    data = read(root_group[name])
                    
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
        
        return FieldMap(xgrid, ygrid, zgrid, Edata, Bdata, 
                        geometry, field_scale, harmonic, phi0, master_param,
                        ele_anchor_pt, curved_ref_frame, interpolation_order)
    end
end
