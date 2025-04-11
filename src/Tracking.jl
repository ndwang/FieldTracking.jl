module Tracking

using StaticArrays
using DifferentialEquations
using ..Interpolation: InterpolatedFieldMap, get_fields

# For a single charged particle, the ODE is:
#   dpos/dt = p/m
#   dp/dt = q * (E + (p/m) × B)
# We'll assemble these into a system of 6 variables: [x, px, y, py, z, pz].

"""
    dynamics!(du, u, p, t)

Computes the time derivative of (pos, momentum) at state `u = [x, px, y, py, z, pz]`. 
Parameters `p` is a NamedTuple or struct containing:
- `q::Float64` (charge)
- `m::Float64` (mass)
- `fieldmap::InterpolatedFieldMap`
"""
function dynamics!(du, u, p, t)
    # Unpack state
    x, px, y, py, z, pz = u[1], u[2], u[3], u[4], u[5], u[6]

    # Unpack parameters
    q = p.q
    m = p.m
    fieldmap = p.fieldmap

    # Get fields at (x,y,z)
    E, B = get_fields(fieldmap, x, y, z)

    # Calculate velocities from momenta
    vx, vy, vz = px/m, py/m, pz/m
    
    # dp/dt = q*(E + v×B)
    momentum = SVector{3}(px, py, pz)
    vel = momentum / m
    force = q * (E + cross(vel, B))

    # du/dt
    du[1] = vx
    du[2] = force[1]
    du[3] = vy
    du[4] = force[2]
    du[5] = vz
    du[6] = force[3]
end

"""
    track_particle!(pos, momentum, q, m, ifm::InterpolatedFieldMap; tspan=(0.0,1e-6), solver=Tsit5())

Track a single particle with initial position `pos` and momentum `momentum`,
given charge `q`, mass `m`, and an `InterpolatedFieldMap`. 
We create an ODEProblem in DifferentialEquations.jl, solve it, and return
the solution. The user can specify solver, timesteps, etc.
"""
function track_particle!(pos::SVector{3,Float64},
                         momentum::SVector{3,Float64},
                         q::Float64,
                         m::Float64,
                         ifm::InterpolatedFieldMap;
                         tspan=(0.0, 1e-6),
                         solver=Tsit5(),
                         save_everystep=true,
                         kwargs...)
    # Build initial condition
    u0 = [pos[1], momentum[1], pos[2], momentum[2], pos[3], momentum[3]]

    # Parameters
    p = (q=q, m=m, fieldmap=ifm)

    # Define ODE problem
    prob = ODEProblem(dynamics!, u0, tspan, p)

    # Solve with chosen method (default: Tsit5())
    sol = solve(prob, solver; save_everystep=save_everystep, kwargs...)

    return sol
end

end # module
