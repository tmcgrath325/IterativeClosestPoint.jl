## Adapted from Protein3DAlign.jl (Tim Holy)
@enum NWTrace::UInt8 UNDET LEFT TOP DIAG

Base.show(io::IO, t::NWTrace) = show(io, t == UNDET ? 'U' :
                                         t == LEFT ? 'L' :
                                         t == TOP ? 'T' : 'D')
Base.show(io::IO, ::typeof(MIME("text/plain")), t::NWTrace) = show(io, t)

"""
    PairIterator(P, start)

Construct an iterator over the aligned pairs, visiting only `DIAG` elements of `P`.
`start` is always visited regardless of the value of `P`, and it is implicitly treated as
`DIAG`.

# Examples

```jldoctest
julia> using Protein3DAlign: UNDET, LEFT, TOP, DIAG

julia> P = [UNDET LEFT LEFT LEFT;
            TOP   LEFT LEFT LEFT;
            TOP   DIAG LEFT LEFT;
            TOP   TOP  LEFT UNDET];

julia> for pair in Protein3DAlign.PairIterator(P, (4, 4))
           @show pair
       end
pair = (4, 4)
pair = (3, 2)
```

This returned sequence arises from the following:

- the iterator starts at (4, 4) and returns it regardless of the value of `P`
- the iterator moves to (3, 3) (treating the starting point as if it were `DIAG`).
  Since the value of `P` here is `LEFT`, it moves on to (3, 2).
  This is `DIAG`, so it returns this pair.
- The iterator then progresses to (2, 1) `TOP` followed by (1, 1) `UNDET`, which causes termination.
  None of these items are returned.
"""
struct PairIterator{M<:AbstractMatrix{NWTrace}}
    P::M
    start::Tuple{Int,Int}
end

Base.IteratorSize(::Type{<:PairIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:PairIterator}) = Tuple{Int,Int}

Base.iterate(iter::PairIterator) = iterate(iter, iter.start)

function Base.iterate(iter::PairIterator, (i,j)::Tuple{Int,Int})
    item = (i, j)
    P = iter.P
    ax1f, ax2f = map(first, axes(P))
    i -= 1; j -= 1    # treat the current item as DIAG (on entry, it might not have been set yet)
    while i >= ax1f && j >= ax2f
        p = P[i, j]
        p == UNDET && break
        p == DIAG && return item, (i, j)
        if p == TOP
            i -= 1
        else # p == LEFT
            j -= 1
        end
    end
    return item, nothing
end
Base.iterate(iter::PairIterator, ::Nothing) = nothing


mutable struct DistanceNeedlemanWunsch{T}
    similarity::Matrix{T}
    score::Matrix{T}
    trace::Matrix{NWTrace}
    gapopen::T
    gapextend::T
end

DistanceNeedlemanWunsch(similarity, gapopen=-0.6, gapextend=0.0
    ) = DistanceNeedlemanWunsch(similarity, 
                                 zeros(size(similarity).+1...), 
                                 similar(zeros(size(similarity).+1...), NWTrace), 
                                 gapopen, 
                                 gapextend)

function fill!(nw::DistanceNeedlemanWunsch)
    similarity = nw.similarity
    score = nw.score
    trace = nw.trace
    gapopen, gapextend = nw.gapopen, nw.gapextend

    m,n = size(similarity)

    # First column & row
    for i in 1:m
        score[i+1,1] = score[i,1] + (i == 1 ? gapopen : gapextend)
        trace[i+1,1] = TOP
    end

    for j in 1:n
        score[1,j+1] = score[1,j] + (j == 1 ? gapopen : gapextend)
        trace[1,j+1] = LEFT
    end

    biasleft = m > n
    for j in 1:n, i in 1:m
        scorediag = score[i,j] + similarity[i,j]     # only distance between points is considered (penalized)
        scoreleft = score[i+1,j] + (i == m ? gapextend : gappenalty(gapopen, gapextend, trace[i+1,j], LEFT))
        scoretop  = score[i,j+1] + (j == n ? gapextend : gappenalty(gapopen, gapextend, trace[i,j+1], TOP))

        if scorediag > max(scoreleft, scoretop)
            score[i+1,j+1] = scorediag
            trace[i+1,j+1] = DIAG
        elseif scoreleft > scoretop || (scoreleft == scoretop && biasleft)
            score[i+1,j+1] = scoreleft
            trace[i+1,j+1] = LEFT
        else
            score[i+1,j+1] = scoretop
            trace[i+1,j+1] = TOP
        end
    end
end

gappenalty(gap_open::Real, gap_extend::Real, p, pnew) = p == pnew ? gap_extend : gap_open

function traceback(F::AbstractMatrix, P::AbstractMatrix, i=nothing, j=nothing)
    ax1, ax2 = axes(P)
    i === nothing && (i = last(ax1))
    j === nothing && (j = last(ax2))
    score = F[i,j]
    # Trace back
    matchedpairs = Tuple{Int,Int}[]
    while i > first(ax1) || j > first(ax2)
        p = P[i,j]
        if p == DIAG
            i -= 1
            j -= 1
            push!(matchedpairs, (i,j))
        elseif p == LEFT
            j -= 1
        else
            i -= 1
        end
    end
    reverse!(matchedpairs)
    return score, matchedpairs
end

traceback(nw::DistanceNeedlemanWunsch, i=nothing, j=nothing) = traceback(nw.score, nw.trace, i, j)