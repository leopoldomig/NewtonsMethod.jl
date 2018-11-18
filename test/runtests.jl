using NewtonsMethod
using Test
using LinearAlgebra, Statistics, Compat, ForwardDiff

@testset "NewtonsMethod.jl" begin
    f1(x) = (x-1.0)^3
    f1′(x) = 3.0*(x-1.0)^2
    x1₀ = 0.1
    @test newtonroot(f1, f1′, x₀ = x1₀) ≈ 1.0
    @test newtonroot(f1, x₀ = x1₀) ≈ 1.0;
end
