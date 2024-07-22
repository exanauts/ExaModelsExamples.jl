# Find the polygon of maximal area, among polygons with nv sides and    
#   diameter d <= 1

#   This is problem 1 in the COPS (Version 3) collection of 
#   E. Dolan and J. More'
#   see "Benchmarking Optimization Software with COPS"
#   Argonne National Labs Technical Report ANL/MCS-246 (2004)

function polygon_model(n; T = Float64, backend = nothing, kwargs...)
    N = div(n, 2)

    c = ExaModels.ExaCore(T; backend = backend)
    r = ExaModels.variable(c, N; lvar = 0, uvar = 1, start = fill(1,N))
    θ = ExaModels.variable(c, N; lvar = 0, uvar = π, start = [i * π / (N - 1) - π / (N - 1) for i in 1:N])
    tuples = Vector{Tuple{Int, Int}}()
    for i in 1:N-1
        for j in i+1:N
            push!(tuples, (i, j))
        end
    end    

    ExaModels.constraint(
        c,
        θ[N] - π
    )

    ExaModels.constraint(
        c,
        r[N] - 0
    )
    
    ExaModels.constraint(
        c,
        θ[i + 1] - θ[i] for i in 1:N-1;
        lcon = 0,
        ucon = Inf
    )

    ExaModels.constraint(
        c,
        r[i]^2 + r[j]^2 - 2 * r[i] * r[j] * cos(θ[i] - θ[j]) - 1 for (i,j) in tuples;
        lcon = -Inf,
        ucon = 0
    )

    ExaModels.objective(c, (r[i] * r[i + 1] * sin(θ[i + 1] - θ[i]) /(-2) for i = 1:Int(N - 1)))

    return ExaModels.ExaModel(c; kwargs...)
end

