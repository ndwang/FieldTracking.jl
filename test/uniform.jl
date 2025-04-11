using Interpolations, FieldTracking, Plots

# Define a simple grid
x = y = z = -10.0:1.0:10.0  # Grid from 0 to 1 in each direction

# Create field values: zero everywhere for E and for Bx, By
E0 = fill(0.0, length(x), length(y), length(z))
B0 = fill(0.0, length(x), length(y), length(z))
Bz = fill(1.0, length(x), length(y), length(z))  # Uniform Bz = 1.0

# Create interpolations
eitp = (
    interpolate((x, y, z), E0, Gridded(Linear())),
    interpolate((x, y, z), E0, Gridded(Linear())),
    interpolate((x, y, z), E0, Gridded(Linear()))
)
bitp = (
    interpolate((x, y, z), B0, Gridded(Linear())),
    interpolate((x, y, z), B0, Gridded(Linear())),
    interpolate((x, y, z), Bz, Gridded(Linear()))
)

# Create the field map
fieldmap = FieldTracking.InterpolatedFieldMap(eitp, bitp)

# Create a particle with initial position and momentum
particle = [0, 1.0, 0, 0, 0, 1.0]

# Track the particle
sol = FieldTracking.track_particle!(particle, 1.0, 1.0, fieldmap; tspan=(0.0, 1))

# Plot the results
plotly()
p = plot(sol, idxs=[(1,3,5)])
display(p)
println("Press Enter to exit...")
readline()