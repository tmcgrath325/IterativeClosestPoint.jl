d0(Lmin) = 1.24 * (Lmin - 15)^(1/3) - 1.8

function similarityscore(x::AbstractMatrix, y::AbstractMatrix; metric=SqEuclidean())
    @assert size(x,1) == size(y,1)
    Lmin = min(size(x,2), size(y,2))
    d₀ = d0(Lmin) 
    sqdist = pairwise(metric, x, y)
    return 1 ./ (1 .+ sqdist ./ (d₀^2))
end

function needlemanwunsch_align(x,y; gapopen=-0.6, gapextend=0.0)
    similarity = similarityscore(x,y)
    nw = DistanceNeedlemanWunsch(similarity, gapopen,gapextend)
    fill!(nw)
    return traceback(nw)[2]
end

function tm_align(x,y; kwargs...)
    return iterate_kabsch(x, y; correspondence=needlemanwunsch_align, kwargs...)
end