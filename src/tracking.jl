# For a single charged particle, the ODE is:
#   dq/dt = p/m
#   dp/dt = charge * (E + (p/m) × B)
# We'll assemble these into a system of 6 variables: [x, px, y, py, z, pz].

"""
    dynamics!(du, u, p, t)

Computes the time derivative of (pos, momentum) at state `u = [x, px, y, py, z, pz]`. 
Parameters `p` is a NamedTuple or struct containing:
- `charge::Float64` (charge)
- `gamma_m::Float64` (mass times relativistic gamma)
- `fieldmap::FieldMap`
"""
function dynamics!(du, u, p, t)
    # Unpack state
    x, px, y, py, z, pz = u[1], u[2], u[3], u[4], u[5], u[6]

    # Unpack parameters
    charge = p.charge
    gamma_m = p.gamma_m
    fieldmap = p.fieldmap

    # Get fields at (x,y,z)
    EM = get_fields(fieldmap, x, y, z)

    # Calculate velocities from momenta
    vx, vy, vz = px/gamma_m, py/gamma_m, pz/gamma_m
    
    # dp/dt = charge*(E + v×B)
    momentum = SVector{3}(px, py, pz)
    vel = momentum / gamma_m
    force = charge * (EM.E + cross(vel, EM.B))

    # du/dt
    du[1] = vx / vz
    du[2] = force[1] / vz
    du[3] = vy / vz
    du[4] = force[2] / vz
    du[5] = 1.0
    du[6] = force[3] / vz
end

"""
    track_particle!(pos, momentum, q, m, fm::FieldMap; tspan=(0.0,1e-6), solver=Tsit5())

Track a single particle with initial position `pos` and momentum `momentum`,
given charge , mass `m`, and an `InterpolatedFieldMap`. 
We create an ODEProblem in DifferentialEquations.jl, solve it, and return
the solution. The user can specify solver, timesteps, etc.
"""
function track_particle!(particle::AbstractVector,
                         charge::Float64,
                         gamma_m::Float64,
                         fm::FieldMap;
                         tspan=(0.0, 1e-6),
                         solver=Tsit5(),
                         save_everystep=true,
                         kwargs...)
    # Parameters
    p = (charge=charge, gamma_m=gamma_m, fieldmap=fm)

    # Define ODE problem
    prob = ODEProblem(dynamics!, particle, tspan, p)

    # Solve with chosen method (default: DPRKN6())
    sol = solve(prob, solver; save_everystep=save_everystep, kwargs...)

    return sol
end
