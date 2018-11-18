module NewtonsMethod

using LinearAlgebra, Statistics, Compat, ForwardDiff

# Declaring the function that applies the newton method
function newtonroot(f, f′; x₀, tolerance = 1E-7, maxiter = 1000)
    x_old = x₀
    normdiff = Inf
    iter = 1
    while normdiff > tolerance && iter <= maxiter
        if f′(x_old) ≈ 0.0
            return (root = nothing, normdiff = nothing, iter = nothing)
        end
        x_new = x_old - f(x_old)/f′(x_old)
        normdiff = norm(x_new - x_old)
        x_old = x_new
        iter = iter + 1
    end
    if iter == maxiter+1
        return (root = nothing, normdiff = nothing, iter = nothing)
    else
        return (root = x_old, normdiff = normdiff, iter = iter)
    end
end

# Applying auto-differentiation
D(f) = x -> ForwardDiff.derivative(f, x)

# Using multiple dispatch
newtonroot(f; x₀, tolerance = 1E-7, maxiter = 1000) = newtonroot(f, D(f), x₀=x₀, tolerance = tolerance, maxiter = maxiter)

# Exporting
export newtonroot

end # module
