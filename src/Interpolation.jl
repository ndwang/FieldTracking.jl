using Interpolations
using StaticArrays

"""
    InterpolatedFieldMap

Wraps your `FieldMap` plus the actual interpolation objects (for E and B).
You may create separate interpolants for Ex, Ey, Ez, or store them in one
4D interpolation object, etc. For demonstration, let's store them as 
six separate interpolants (one for each field component).
"""
struct InterpolatedFieldMap
    Eitp::NTuple{3, Interpolations.GriddedInterpolation}
    Bitp::NTuple{3, Interpolations.GriddedInterpolation}
end

"""
    make_fieldmap_interpolant(fm::FieldMap)

Builds a set of 3D interpolations for E and B from the raw grid data in `fm`.
Returns an `InterpolatedFieldMap`. Supports both rectangular (cartesian) and 
cylindrical coordinate systems.

For axes with dimension equal to 1, the interpolation will drop that axis.
"""
function make_fieldmap_interpolant(fm::FieldMap)::InterpolatedFieldMap
    # For 3D grids, Interpolations.jl wants something like:
    # itp = interpolate((fm.xgrid, fm.ygrid, fm.zgrid), fm.Edata[...], Gridded(Linear()))

    # Determine which axes to include based on their dimensions
    axes = []
    if length(fm.xgrid) > 1
        push!(axes, fm.xgrid)
    end
    if length(fm.ygrid) > 1
        push!(axes, fm.ygrid)
    end
    if length(fm.zgrid) > 1
        push!(axes, fm.zgrid)
    end
    
    # For the 0D case, we need to wrap our functions in a callable object
    # that matches the interface expected by get_fields
    if length(axes) == 0
        Ex_func = (i,j,k) -> fm.Edata[1,1,1,1]
        Ey_func = (i,j,k) -> fm.Edata[1,1,1,2]
        Ez_func = (i,j,k) -> fm.Edata[1,1,1,3]
        
        Bx_func = (i,j,k) -> fm.Bdata[1,1,1,1]
        By_func = (i,j,k) -> fm.Bdata[1,1,1,2]
        Bz_func = (i,j,k) -> fm.Bdata[1,1,1,3]
        
        return InterpolatedFieldMap((Ex_func, Ey_func, Ez_func), (Bx_func, By_func, Bz_func))
    end

    # Create interpolation objects for each field component
    function create_interpolant(data)
        if length(axes) == 3
            # Full 3D interpolation
            return interpolate((fm.xgrid, fm.ygrid, fm.zgrid), data, Gridded(Linear()))
        elseif length(axes) == 2
            # 2D interpolation - determine which axes to use
            if length(fm.xgrid) == 1
                return interpolate((fm.ygrid, fm.zgrid), data[1,:,:], Gridded(Linear()))
            elseif length(fm.ygrid) == 1
                return interpolate((fm.xgrid, fm.zgrid), data[:,1,:], Gridded(Linear()))
            else # length(fm.zgrid) == 1
                return interpolate((fm.xgrid, fm.ygrid), data[:,:,1], Gridded(Linear()))
            end
        elseif length(axes) == 1
            # 1D interpolation - determine which axis to use
            if length(fm.xgrid) > 1
                return interpolate((fm.xgrid,), data[:,1,1], Gridded(Linear()))
            elseif length(fm.ygrid) > 1
                return interpolate((fm.ygrid,), data[1,:,1], Gridded(Linear()))
            else # length(fm.zgrid) > 1
                return interpolate((fm.zgrid,), data[1,1,:], Gridded(Linear()))
            end
        else
            # 0D case (all axes have dimension 1) - just return constant values
            return (i,j,k) -> data[1,1,1]
        end
    end

    # Create interpolants for each field component
    Ex = create_interpolant(fm.Edata[:,:,:,1])
    Ey = create_interpolant(fm.Edata[:,:,:,2])
    Ez = create_interpolant(fm.Edata[:,:,:,3])
    
    Bx = create_interpolant(fm.Bdata[:,:,:,1])
    By = create_interpolant(fm.Bdata[:,:,:,2])
    Bz = create_interpolant(fm.Bdata[:,:,:,3])


    
    return InterpolatedFieldMap((Ex, Ey, Ez), (Bx, By, Bz))
end

"""
    get_fields(ifm::InterpolatedFieldMap, i::Float64, j::Float64, k::Float64)

Returns (E, B) as SVector{3, Float64} at position (i, j, k).
For cartesian coordinates, i is x, j is y, k is z.
For cylindrical coordinates, x is r, y is θ, and z is z.

Handles interpolation objects with reduced dimensionality when some axes have dimension 1.
"""
function get_fields(ifm::InterpolatedFieldMap, i::Float64, j::Float64, k::Float64)
    Ex, Ey, Ez = ifm.Eitp
    Bx, By, Bz = ifm.Bitp
    
    # The interpolants are already set up to handle the correct dimensionality
    return SVector{3}(Ex(i,j,k), Ey(i,j,k), Ez(i,j,k)),
           SVector{3}(Bx(i,j,k), By(i,j,k), Bz(i,j,k))
end

"""
    get_fields_cartesian(ifm::InterpolatedFieldMap, fm::FieldMap, x::Float64, y::Float64, z::Float64)

Returns (E, B) as SVector{3, Float64} at position (x, y, z) in Cartesian coordinates,
automatically converting from cylindrical if needed.
"""
function get_fields_cartesian(ifm::InterpolatedFieldMap, fm::FieldMap, x::Float64, y::Float64, z::Float64)
    if lowercase(fm.geometry) == "rectangular"
        # Already in Cartesian coordinates
        return get_fields(ifm, x, y, z)
    elseif lowercase(fm.geometry) == "cylindrical"
        # Convert Cartesian (x,y,z) to cylindrical (r,θ,z)
        r = sqrt(x^2 + y^2)
        θ = atan(y, x)  # atan2 equivalent
        
        E, B = get_fields(ifm, r, θ, z)
        
        # Convert field components from cylindrical to Cartesian
        cos_θ, sin_θ = cos(θ), sin(θ)
        
        Ex = E[1] * cos_θ - E[2] * sin_θ
        Ey = E[1] * sin_θ + E[2] * cos_θ
        
        Bx = B[1] * cos_θ - B[2] * sin_θ
        By = B[1] * sin_θ + B[2] * cos_θ
        
        return SVector{3}(Ex, Ey, E[3]),
               SVector{3}(Bx, By, B[3])
    else
        error("Unsupported geometry: $(fm.geometry)")
    end
end
