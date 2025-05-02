using FieldTracking # Assuming your module is FieldTracking
using Plots         # Add Plots to your test environment if needed: using Pkg; Pkg.add("Plots")
# using GR          # Example backend, add if needed: using Pkg; Pkg.add("GR")

# --- Configuration ---
# !!! IMPORTANT: Replace this with the actual path to your HDF5 fieldmap file !!!
test_fieldmap_path = raw"C:\Users\wdong\Downloads\lcls2_solenoid_fieldmesh.h5"
# Set to true to generate and display the plot
enable_plotting = true # Set to true to see the plot
# --- End Configuration ---

if !isfile(test_fieldmap_path)
    error("Fieldmap file not found at '$test_fieldmap_path'. Please provide a valid path.")
else
    # 1. Load the FieldMap
    println("Loading FieldMap from: ", test_fieldmap_path)
    fm = FieldTracking.FieldMap(test_fieldmap_path)
    println("FieldMap loaded successfully.")
    println("  Grid Geometry: ", fm.gridGeometry)
    println("  Grid Origin Offset: ", fm.gridOriginOffset)
    println("  Grid Size: ", fm.gridSize)
    println("  Grid Spacing: ", fm.gridSpacing)
    println("  Grid Lower Bound: ", fm.gridLowerBound)


    # 2. Determine the physical bounds of the z-axis from metadata
    # Grid indices start from gridLowerBound (often [1,1,1])
    # Physical coordinate = origin + (index - 1) * spacing
    z_origin = fm.gridOriginOffset[3]
    z_spacing = fm.gridSpacing[3]
    z_points = fm.gridSize[3]
    z_lower_bound_index = fm.gridLowerBound[3] # Starting index for z

    # Calculate the physical z-coordinates corresponding to the grid points
    z_grid_coords = ( z_lower_bound_index .+ (0:z_points-1) ) .* z_spacing
    z_min = minimum(z_grid_coords)
    z_max = maximum(z_grid_coords)

    println("Calculated Z-axis bounds: [", z_min, ", ", z_max, "]")

    # 3. Generate points along the z-axis for evaluation (e.g., 100 points)
    num_plot_points = 100
    z_eval_coords = range(z_min, stop=z_max, length=num_plot_points)
    println("Generating ", num_plot_points, " evaluation points along Z-axis.")

    # 4. Calculate the Bz component at each point (x=0, y=0)
    bz_values = zeros(Float64, num_plot_points) # Store real part
    println("Calculating Bz values...")

    if fm.Bitp === nothing
         @warn "Magnetic field interpolator (Bitp) is nothing. Cannot calculate Bz."
         # Fill with NaN or skip plotting/testing Bz
         fill!(bz_values, NaN)
    else
        for (i, z) in enumerate(z_eval_coords)
            # get_fields expects absolute coordinates
            em_field = FieldTracking.get_fields(fm, 0.0, 0.0, z)
            # Extract the real part of the Bz component
            bz_values[i] = real(em_field.B[3])
        end
        println("Bz calculation complete.")
    end


    # 5. Plot Bz vs z (optional)
    if enable_plotting && fm.Bitp !== nothing
        println("Generating plot...")
        # gr() # Select backend (optional, depends on setup)
        p = plot(z_eval_coords, bz_values,
                 xlabel="Z Coordinate (m)",
                 ylabel="Bz (Real Part) (T)", # Assuming Tesla, adjust if needed based on fieldScale
                 title="Bz along Z-axis (x=0, y=0)",
                 label="Bz",
                 legend=:outertopright,
                 linewidth=2)
        display(p) # Show the plot in the default plot pane
        # Or save it: savefig(p, "bz_along_z.png")
        println("Plot displayed.")
        readline()
    elseif enable_plotting
         @warn "Plotting enabled but Bz could not be calculated (Bitp was nothing)."
    else
        println("Plotting is disabled.")
    end

end

println("Script finished.")

# Reminder: To run this script, you might need to:
# 1. Ensure FieldTracking module is accessible.
# 2. Add Plots (and potentially a backend like GR) to your environment.
# 3. Update `test_fieldmap_path` to point to a valid HDF5 file.
# 4. Run this file directly using Julia: julia path/to/plot_fieldmap.jl