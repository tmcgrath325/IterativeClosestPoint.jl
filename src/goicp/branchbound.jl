function global_align(x, y;
                      nsplits=2, initblock=nothing,
                      spacefun=UncertaintyRegion, objfun=alignment_objective, tformfun=AffineMap,
                      atol=0.1, rtol=0, maxblocks=5e8, maxsplits=Inf, maxevals=Inf, maxstagnant=Inf, threads=false)
    if isodd(nsplits)
        throw(ArgumentError("`nsplits` must be even"))
    end
    if dims(gmmx) != dims(gmmy)
        throw(ArgumentError("Dimensionality of the point sets must be equal"))
    end

    # initialization
    if isnothing(initblock)
        
        initblock = blockfun()
    end
    ndims = dims(initblock)
    ub, bestloc = initblock.upperbound, initblock.center    # best-so-far objective value and transformation
    t = promote_type(numbertype(gmmx), numbertype(gmmy))
    lb = typemin(t)
    pq = PriorityQueue{Block{ndims, t}, Tuple{t,t}}()
    enqueue!(pq, initblock, (initblock.lowerbound, initblock.upperbound))
    
    # split cubes until convergence
    ndivisions = 0
    sinceimprove = 0
    evalsperdiv = length(gmmx)*length(gmmy)*nsplits^ndims
    while !isempty(pq)
        if (length(pq) > maxblocks) || (ndivisions*evalsperdiv > maxevals) || (sinceimprove > maxstagnant) || (ndivisions > maxsplits)
            break
        end
        ndivisions += 1
        sinceimprove += 1

        # take the block with the lowest lower bound
        bl, (lb, blub) = dequeue_pair!(pq)

        # if the best solution so far is close enough to the best possible solution, end
        if abs((ub - lb)/lb) < rtol || abs(ub-lb) < atol
            return GMAlignmentResult(gmmx, gmmy, ub, lb, tformfun(bestloc...), bestloc, ndivisions*evalsperdiv, ndivisions, length(pq), sinceimprove)
        end

        # split up the block into `nsplits` smaller blocks across each dimension
        subrngs = subranges(bl.ranges, nsplits)
        sblks = fill(Block{ndims,t}(), nsplits^ndims)
        for i=1:length(subrngs)
            sblks[i] = blockfun(gmmx, gmmy, subrngs[i], pσ, pϕ, rot, trl)
        end

        # reset the upper bound if appropriate
        minub, ubidx = findmin([sblk.upperbound for sblk in sblks])
        if minub < ub
            ub, bestloc = local_align(gmmx, gmmy, sblks[ubidx], pσ, pϕ; objfun=objfun, rot=rot, trl=trl)
            sinceimprove = 0
        end

        # only add sub-blocks to the queue if they present possibility for improvement
        for sblk in sblks
            if sblk.lowerbound < ub
                enqueue!(pq, sblk, (sblk.lowerbound, sblk.upperbound))
            end
        end
    end
    if isempty(pq)
        return GMAlignmentResult(gmmx, gmmy, ub, lb, tformfun(bestloc...), bestloc, ndivisions*evalsperdiv, ndivisions, length(pq), sinceimprove)
    else
        return GMAlignmentResult(gmmx, gmmy, ub, dequeue_pair!(pq)[2][1], tformfun(bestloc...), bestloc, ndivisions*evalsperdiv, ndivisions, length(pq), sinceimprove)
    end
end