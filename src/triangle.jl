# Minimize the time taken for a robot arm to travel between two points.

#  This is problem 18 in the COPS (Version 3) collection of 
#   E. Dolan and J. More
#   see "Benchmarking Optimization Software with COPS"
#   Argonne National Labs Technical Report ANL/MCS-246 (2004)

include(joinpath("..", "data", "triangle.jl"))

function triangle_model(x0 = xe, TRIS::Vector{Int64} = Tr, Const::Vector{Int64} = Constants; T = Float64, backend = nothing, kwargs...)
    τ = 0.0
    n = length(x0)
    N = Int(div(n, 2))
    E = Int(div(length(TRIS), 3))

    lvar = -Inf * ones(n)
    lvar[Const] = x0[Const]
    lvar[Const .+ N] = x0[Const .+ N]

    uvar = Inf * ones(n)
    uvar[Const] = x0[Const]
    uvar[Const .+ N] = x0[Const .+ N]

    c = ExaModels.ExaCore(T; backend = backend)

    x = ExaModels.variable(c, n; lvar = lvar, uvar = uvar, start = x0)

    ExaModels.objective(
        c,
        
        (
            (1 * x[TRIS_E] - x[TRIS])^2 +
            (2 * x[TRIS_2E] - x[TRIS_E] - x[TRIS])^2 / 3
        ) 
        / (
        2 * (
            2 * (
                (x[TRIS_E] - x[TRIS]) * (x[TRIS_2E_N] - x[TRIS_N]) -
                (x[TRIS_2E] - x[TRIS]) * (x[TRIS_E_N] - x[TRIS_N])
            ) / sqrt(3)
        )
    ) for (TRIS, TRIS_E, TRIS_2E, TRIS_N, TRIS_E_N, TRIS_2E_N) in [(TRIS[e], TRIS[e + E], TRIS[e + 2*E], TRIS[e]+N, TRIS[e + E]+N, TRIS[e + 2*E]+N) for e in 1:E]
    )

    ExaModels.objective(
        c,
        
        (
            (1 * x[TRIS_E_N] - x[TRIS_N])^2 +
            (2 * x[TRIS_2E_N] - x[TRIS_E_N]  - x[TRIS_N])^2 / 3
        ) 
        / (
        2 * (
            2 * (
                (x[TRIS_E] - x[TRIS]) * (x[TRIS_2E_N] - x[TRIS_N]) -
                (x[TRIS_2E] - x[TRIS]) * (x[TRIS_E_N] - x[TRIS_N])
            ) / sqrt(3)
        )
        ) for (TRIS, TRIS_E, TRIS_2E, TRIS_N, TRIS_E_N, TRIS_2E_N) in [(TRIS[e], TRIS[e + E], TRIS[e + 2*E], TRIS[e]+N, TRIS[e + E]+N, TRIS[e + 2*E]+N) for e in 1:E]
    )

    ExaModels.constraint(
        c,
        2 * (
                (x[TRIS_E] - x[TRIS]) * (x[TRIS_2E_N] - x[TRIS_N]) -
                (x[TRIS_2E] - x[TRIS]) * (x[TRIS_E_N] - x[TRIS_N])
            ) / sqrt(3) 
        for (TRIS, TRIS_E, TRIS_2E, TRIS_N, TRIS_E_N, TRIS_2E_N) in [(TRIS[e], TRIS[e + E], TRIS[e + 2*E], TRIS[e]+N, TRIS[e + E]+N, TRIS[e + 2*E]+N) for e in 1:E];
        lcon = τ,
        ucon = Inf
    )

    return ExaModels.ExaModel(c; kwargs...)
end

triangle_deer_model() = triangle_model(xe_deer, TRIS_deer, Const_deer)
triangle_pacman_model() = triangle_model(xe_pacman, TRIS_pacman, Const_pacman)
triangle_turtle_model() = triangle_model(xe_turtle, TRIS_turtle, Const_turtle)
