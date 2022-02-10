
# All matrices are DxN, where N is the number of positions and D is the dimensionality

# Here, P is the probe (to be rotated) and Q is the refereence
# https://en.wikipedia.org/wiki/Kabsch_algorithm

# assuming P and Q are already centered at the origin
# returns the rotation for alignment
function kabsch_centered(P,Q)
    @assert size(P) == size(Q)
    H = P*Q'
    D = Matrix{Float64}(I,size(H,1), size(H,2))
    U,Σ,V = svd(H)
    D[end] = sign(det(V*U'))
    return LinearMap(V * D * U')
end

# translation moving centroid to origin
center_translation(A) = Translation(-mean(A, dims=2))

# transform DxN matrices
function (tform::Translation)(A::AbstractMatrix)
    return hcat([tform(A[:,i]) for i=1:size(A,2)]...)
end

# P and Q are not necessarily centered
# returns the transformation for alignment
function kabsch(P,Q)
    centerP, centerQ = center_translation(P), center_translation(Q)
    R = kabsch_centered(centerP(P), centerQ(Q))
    return inv(centerQ) ∘ R ∘ centerP
end
