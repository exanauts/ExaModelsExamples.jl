function distillation_column_model(N = 3; T = Float64, backend = nothing, kwargs...)
    NT = 30
    FT = 17
    Ac = 0.5
    At = 0.25
    Ar = 1.0
    D = 0.2
    F = 0.4
    ybar = 0.8958
    ubar = 2.0
    alpha = 1.6
    dt = 10 / N
    xAf = 0.5

    c = ExaModels.ExaCore(T; backend = backend)

    xA = ExaModels.variable(c, 0:N, 0:NT+1; start = 0.5)
    yA = ExaModels.variable(c, 0:N, 0:NT+1; start = 0.5)
    u = ExaModels.variable(c, 0:N; start = 1.0)
    V = ExaModels.variable(c, 0:N; start = 1.0)
    L2 = ExaModels.variable(c, 0:N; start = 1.0)

    ExaModels.objective(c, (yA[t, 1] - ybar)^2 for t = 0:N)
    ExaModels.objective(c, (u[t] - ubar)^2 for t = 0:N)

    ExaModels.constraint(c, xA[0, i] - 0.5 for i in 0:NT+1)
    ExaModels.constraint(
        c,
        (xA[t, 0] - xA[t-1, 0]) / dt - (1 / Ac) * (yA[t, 1] - xA[t, 0]) for t = 1:N
    )
    ExaModels.constraint(
        c,
        (xA[t, i] - xA[t-1, i]) / dt -
        (1 / At) * (u[t] * D * (yA[t, i-1] - xA[t, i]) - V[t] * (yA[t, i] - yA[t, i+1])) for
        t in 1:N, i in 1:FT-1
    )
    ExaModels.constraint(
        c,
        (xA[t, FT] - xA[t-1, FT]) / dt -
        (1 / At) * (
            F * xAf + u[t] * D * xA[t, FT-1] - L2[t] * xA[t, FT] -
            V[t] * (yA[t, FT] - yA[t, FT+1])
        ) for t = 1:N
    )
    ExaModels.constraint(
        c,
        (xA[t, i] - xA[t-1, i]) / dt -
        (1 / At) * (L2[t] * (yA[t, i-1] - xA[t, i]) - V[t] * (yA[t, i] - yA[t, i+1])) for
        t in 1:N, i in FT+1:NT
    )
    ExaModels.constraint(
        c,
        (xA[t, NT+1] - xA[t-1, NT+1]) / dt -
        (1 / Ar) * (L2[t] * xA[t, NT] - (F - D) * xA[t, NT+1] - V[t] * yA[t, NT+1]) for
        t = 1:N
    )
    ExaModels.constraint(c, V[t] - u[t] * D - D for t = 0:N)
    ExaModels.constraint(c, L2[t] - u[t] * D - F for t = 0:N)
    ExaModels.constraint(
        c,
        yA[t, i] * (1 - xA[t, i]) - alpha * xA[t, i] * (1 - yA[t, i]) for t in 0:N, i in 0:NT+1
    )

    return ExaModels.ExaModel(c; kwargs...)
end
