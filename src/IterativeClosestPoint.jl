module IterativeClosestPoint

using NearestNeighbors
using Distances
using CoordinateTransformations
using Rotations
using StaticArrays
using Statistics
using LinearAlgebra

export kabsch
export iterate_kabsch

include("kabsch.jl")
include("correspondence.jl")
include("iterate.jl")
include("goicp/uncertaintyregion.jl")

include("tmalign/align.jl")
include("tmalign/tmscore.jl")
include("tmalign/needlemanwunsch.jl")


end
