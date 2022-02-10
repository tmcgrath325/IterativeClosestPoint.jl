function tm_score(P, Q, matches; metric=SqEuclidean())
    matchedP, matchedQ = matched_points(P, Q, matches)
    LN = size(Q,2)
    tform = kabsch(matchedP, matchedQ)
    return 1/LN * sum(1 ./ (1 .+ colwise(metric, tform(matchedP), matchedQ) ./ d0(LN)^2))
end

