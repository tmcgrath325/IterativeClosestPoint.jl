struct AbstractPointSet{N,T} <: AbstractVector{T}

struct PointSet{N,T} <: AbstractPointSet{N,T}
    coords::Vector{SVector{N,T}}
end

getindex(ps::PointSet) = getindex(ps.coords)
length(ps::PointSet) = length(ps.coords)