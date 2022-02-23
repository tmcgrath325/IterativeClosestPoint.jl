using IterativeClosestPoint
using CoordinateTransformations
using Rotations
using Distances
using StaticArrays
using LinearAlgebra
using Random
using Test

function randomtform(; s=10, maxangle=π)
    randomaxis = normalize(rand(3) .- 0.5)
    randomangle = (rand() .- 0.5) * 2 * maxangle
    randomrotation = AngleAxis(randomangle, randomaxis...)
    randomtranslation = s * (rand(3) .- 0.5)
    return AffineMap(randomrotation, randomtranslation)
end

@testset "Kabsch" begin
    for i=1:10
        x = SMatrix{3,20}(20 * (rand(3, 20) .- 0.5))
        randtform = randomtform()
        y = randtform(x)

        tform = kabsch(x, y)
        @test tform.linear ≈ randtform.linear
        @test tform.translation ≈ randtform.translation
    end
end

@testset "point-to-point ICP" begin
     # pure shuffling of points (no transformation)
    for i=1:10
        x = SMatrix{3,20}(20 * (rand(3, 20) .- 0.5))
        perm = randperm(20)
        y = (hcat([x[:,i] for i in perm]...))

        matches = iterate_kabsch(y,x)
        @test perm == [m[2] for m in matches]
    end

    # small-rotation with translation
    for i=1:10
        x = SMatrix{3,20}(20 * (rand(3, 20) .- 0.5))
        randtform = randomtform(maxangle = π/4)
        y = randtform(x)
        perm = randperm(20)
        y = (hcat([y[:,i] for i in perm]...))

        matches = iterate_kabsch(y,x)
        @test perm == [m[2] for m in matches]
    end
end
