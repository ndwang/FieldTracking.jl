using FieldTracking

test_fieldmap_path = raw"C:\Users\wdong\Downloads\lcls2_solenoid_fieldmesh.h5"
fm = FieldTracking.FieldMap(test_fieldmap_path)

# Create a particle with initial position and momentum
particle = [0., 0., 0., 0., 0., 1.0]

# Track the particle
sol = FieldTracking.track_particle!(particle, 1.0, 0.1, fm; tspan=(0.0, 0.48))
print(sol[end])