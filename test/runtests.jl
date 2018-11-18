using NewtonsMethod
using Test
using LinearAlgebra, Statistics, Compat, ForwardDiff

@testset "NewtonsMethod.jl" begin
    #Several @test for the root of a known functions
    #1
    f1(x) = (x-1.0)^3
    f1′(x) = 3.0*(x-1.0)^2
    x1₀ = 0.1
    @test newtonroot(f1, f1′, x₀ = x1₀).root ≈ 0.9999998168636032
    @test newtonroot(f1,      x₀ = x1₀).root ≈ 0.9999998168636032
    @test newtonroot(f1, f1′, x₀ = BigFloat(x1₀), tolerance = 1E-80).root ≈ BigFloat(1.0)
    @test newtonroot(f1,      x₀ = BigFloat(x1₀), tolerance = 1E-80).root ≈ BigFloat(1.0)
    #2
    f2(x) = log(x)
    f2′(x) = 1.0/x
    x2₀ = 2.0
    @test newtonroot(f2, f2′, x₀ = x2₀).root ≈ 1.0
    @test newtonroot(f2,      x₀ = x2₀).root ≈ 1.0
    @test newtonroot(f2, f2′, x₀ = BigFloat(x2₀), tolerance = 1E-50).root ≈ BigFloat(1.0)
    @test newtonroot(f2,      x₀ = BigFloat(x2₀), tolerance = 1E-50).root ≈ BigFloat(1.0)
    #3
    f3(x) = (sin(x)/x)
    f3′(x) = (x*cos(x) - sin(x))/(x^2)
    x3₀ = 1.0
    @test newtonroot(f3, f3′, x₀ = x3₀).root ≈ 3.141592653589793
    @test newtonroot(f3,      x₀ = x3₀).root ≈ 3.141592653589793
    @test newtonroot(f3, f3′, x₀ = BigFloat(x3₀), tolerance = 1E-50).root ≈ pi
    @test newtonroot(f3,      x₀ = BigFloat(x3₀), tolerance = 1E-50).root ≈ pi
    #4
    f4(x) = (x+10.0)*(x^2+2.0)
    f4′(x) = 3x^2 + 20x + 2.0
    x4₀ = 10.0
    @test newtonroot(f4, f4′, x₀ = x4₀).root ≈ -10.0
    @test newtonroot(f4,      x₀ = x4₀).root ≈ -10.0
    @test newtonroot(f4, f4′, x₀ = BigFloat(x4₀), tolerance = 1E-100).root ≈ BigFloat(-10.0)
    @test newtonroot(f4,      x₀ = BigFloat(x4₀), tolerance = 1E-100).root ≈ BigFloat(-10.0)

    #Some @test of non-convergence
    #5
    f5(x) = x^2 + 2.0
    f5′(x) = 2x
    x5₀ = 10.0
    @test newtonroot(f5, f5′, x₀ = x5₀).root == nothing
    @test newtonroot(f5,      x₀ = x5₀).root == nothing
    @test newtonroot(f5, f5′, x₀ = BigFloat(x5₀), tolerance = 1E-100).root == nothing
    @test newtonroot(f5,      x₀ = BigFloat(x5₀), tolerance = 1E-100).root == nothing
    #6
    f6(x) = 2.0
    f6′(x) = 0.0
    x6₀ = 10.0
    @test newtonroot(f6, f6′, x₀ = x6₀).root == nothing
    @test newtonroot(f6,      x₀ = x6₀).root == nothing
    @test newtonroot(f6, f6′, x₀ = BigFloat(x6₀), tolerance = 1E-100).root == nothing
    @test newtonroot(f6,      x₀ = BigFloat(x6₀), tolerance = 1E-100).root == nothing

    #Some @test for maxiter
    #1
    @test newtonroot(f1, f1′, x₀ = x1₀, maxiter = 10).root == nothing
    @test newtonroot(f1,      x₀ = x1₀, maxiter = 10).root == nothing
    @test newtonroot(f1, f1′, x₀ = BigFloat(x1₀), tolerance = 1E-80, maxiter = 100).root == nothing
    @test newtonroot(f1,      x₀ = BigFloat(x1₀), tolerance = 1E-80, maxiter = 100).root == nothing
    #4
    @test newtonroot(f4, f4′, x₀ = x4₀, maxiter = 10).root == nothing
    @test newtonroot(f4,      x₀ = x4₀, maxiter = 10).root == nothing
    @test newtonroot(f4, f4′, x₀ = BigFloat(x4₀), tolerance = 1E-100, maxiter = 100).root ≈ BigFloat(-10.0)
    @test newtonroot(f4,      x₀ = BigFloat(x4₀), tolerance = 1E-100, maxiter = 100).root ≈ BigFloat(-10.0)

    #Some @test for tolerance and BigFloat
    #2
    @test !(newtonroot(f2, f2′, x₀ = x2₀, tolerance = 1E-2).root ≈ 1.0)
    @test !(newtonroot(f2,      x₀ = x2₀, tolerance = 1E-2).root ≈ 1.0)
    #3
    @test !(BigFloat(0.0) + newtonroot(f3, f3′, x₀ = x3₀, tolerance = 1E-50).root ≈ pi)
    @test !(BigFloat(0.0) + newtonroot(f3,      x₀ = x3₀, tolerance = 1E-50).root ≈ pi)
    @test !(newtonroot(f3, f3′, x₀ = BigFloat(x3₀), tolerance = 1E-10).root ≈ pi)
    @test !(newtonroot(f3,      x₀ = BigFloat(x3₀), tolerance = 1E-10).root ≈ pi)
end
