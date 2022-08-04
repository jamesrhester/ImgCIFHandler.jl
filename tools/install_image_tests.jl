# Run this script whenever installing/updating image_tests.jl
# From a prompt, type:

# julia install_image_tests.jl
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
