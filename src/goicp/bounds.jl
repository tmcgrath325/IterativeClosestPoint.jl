const sqrt3 = sqrt(3.)

# upper bound on distance between x and y in the uncertainty region centered at R, T
function upperdist(x, y, R, T)
    return norm(R*x .+ T .- y)
end

# greatest distance that x can be moved within the uncertainty region
# this can be used as a less expensive, but looser, calculation of the lowerbound; ldist = uppderdist - uncertaintydist
function uncertaintydist(x, σᵣ, σₜ)
    γᵣ = 2*sin(min(sqrt3*σᵣ, π/2))*norm(x)
    γₜ = sqrt3 * σₜ
    return γᵣ + γₜ
end

function lowerdist(x, y, R::AngleAxis, T, σᵣ, σₜ)
    # return Inf for bounds if the rotation lies outside the π-sphere
    if R.theta > π
        inf = typemax(promote_type(eltype(x), eltype(y)))
        return inf, inf
    end

    # magnitudes of x and (y-T)
    xnorm, ynorm = norm(x), norm(y-T)

    # α is the angle between R*x (the initial rotation), the origin, and y-T
    if xnorm*ynorm == 0
        cosα = one(promote_type(eltype(x), eltype(y)))
    else
        cosα = dot(R*x, y-T)/(xnorm*ynorm)
    end

    # β is the largest allowed angle of rotation
    cosβ = cos(min(sqrt3*σᵣ, π))

    # if the α > β, then the two points can be rotated within the uncertainty region to be colinear with the origin
    if cosα >= cosβ
        lbdist = max(abs(xnorm-ynorm) - sqrt3*σₜ, 0)
    # otherwise, use the law of cosines to calculate the minimum distance 
    else
        lbdist = try max(√(xnorm^2 + ynorm^2 - 2*xnorm*ynorm*(cosα*cosβ+√((1-cosα^2)*(1-cosβ^2)))) - sqrt3*σₜ, 0)  # law of cosines
        catch e     # when the argument for the square root is negative (within machine precision of 0, usually)
            zero(promote_type(eltype(x), eltype(y)))
        end
    end

    # return the lowerbound on distance between x and y within the uncertainty region
    return lbdist
end

# these let you modify the objective function, which takes x, y, and dist as arguments
upperbound(x, y, R, T, objective=p->p[3]) = objective(x,y,upperdist(x,y,R,T))
lowerbound(x, y, R, T, σᵣ, σₜ, objective=p->p[3]) = objective(x,y,lowerdist(x,y,R,T,σᵣ,σₜ))

