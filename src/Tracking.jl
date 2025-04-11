module Tracking

using StaticArrays
using DifferentialEquations
using ..Interpolation: InterpolatedFieldMap, get_fields

# For a single charged particle, the ODE is:
#   dpos/dt = vel
#   dvel/dt = (q/m) * (E + v × B)
# We'll assemble these into a system of 6 variables: [x, y, z, vx, vy, vz].

"""
    dynamics!(du, u, p, t)

Computes the time derivative of (pos, vel) at state `u = [x, y, z, vx, vy, vz]`. 
Parameters `p` is a NamedTuple or struct containing:
- `q::Float64` (charge)
- `m::Float64` (mass)
- `fieldmap::InterpolatedFieldMap`
"""
function dynamics!(du, u, p, t)
    # Unpack state
    x, y, z = u[1], u[2], u[3]
    vx, vy, vz = u[4], u[5], u[6]

    # Unpack parameters
    q = p.q
    m = p.m
    fieldmap = p.fieldmap

    # Get fields at (x,y,z)
    E, B = get_fields(fieldmap, x, y, z)

    # dv/dt = (q/m)*(E + v×B)
    vel = SVector{3}(vx, vy, vz)
    accel = (q/m) * (E + cross(vel, B))

    # du/dt
    du[1] = vx
    du[2] = vy
    du[3] = vz
    du[4] = accel[1]
    du[5] = accel[2]
    du[6] = accel[3]
end

"""
    track_particle!(pos, vel, q, m, ifm::InterpolatedFieldMap; tspan=(0.0,1e-6), solver=Tsit5())

Track a single particle with initial position `pos` and velocity `vel`,
given charge `q`, mass `m`, and an `InterpolatedFieldMap`. 
We create an ODEProblem in DifferentialEquations.jl, solve it, and return
the solution. The user can specify solver, timesteps, etc.
"""
function track_particle!(pos::SVector{3,Float64},
                         vel::SVector{3,Float64},
                         q::Float64,
                         m::Float64,
                         ifm::InterpolatedFieldMap;
                         tspan=(0.0, 1e-6),
                         solver=Tsit5(),
                         save_everystep=true,
                         kwargs...)
    # Build initial condition
    u0 = [pos[1], pos[2], pos[3], vel[1], vel[2], vel[3]]

    # Parameters
    p = (q=q, m=m, fieldmap=ifm)

    # Define ODE problem
    prob = ODEProblem(dynamics!, u0, tspan, p)

    # Solve with chosen method (default: Tsit5())
    sol = solve(prob, solver; save_everystep=save_everystep, kwargs...)

    return sol
end

end # module
