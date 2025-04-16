using Interpolations, FieldTracking, Plots

# Define a simple grid
x = y = z = collect(-10.0:1.0:10.0)  # Grid from 0 to 1 in each direction

# Create field values: zero everywhere for E and for Bx, By
B0 = fill(0.0, length(x), length(y), length(z))
Bz = fill(1.0, length(x), length(y), length(z))  # Uniform Bz = 1.0
B = cat(B0, B0, Bz, dims=4)

# Create the field map
fieldmap = FieldTracking.FieldMap(x,y,z;B=B)

# Create a particle with initial position and momentum
particle = [0, 1.0, 0, 0, 0, 1.0]

# Track the particle
sol = FieldTracking.track_particle!(particle, 1.0, 0.1, fieldmap; tspan=(0.0, 1))

# Plot the results
p = plot(sol, idxs=[(1,3,5)])
display(p)
println("Press Enter to exit...")
readline()