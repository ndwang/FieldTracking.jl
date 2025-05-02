module FieldTracking

using StaticArrays, HDF5, LinearAlgebra
using DifferentialEquations: ODEProblem, solve, Tsit5 # Import specific components
using Interpolations: interpolate, Gridded, Linear, AbstractInterpolation # Import specific components

# Load submodules
include("fieldMap.jl")
include("tracking.jl")

# Optionally export what you want the user to access
# export FieldMap, make_fieldmap, track_particle!, ...

"""
FieldTracking.jl

A small package demonstrating how to store 3D field data, interpolate it
using Interpolations.jl, and track a charged particle using ODE solvers
from DifferentialEquations.jl.
"""
end # module