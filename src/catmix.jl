# Catalyst Mixing Problem
# Collocation formulation
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function catmix_model(nh; T = Float64, backend = nothing, kwargs...)
    ne = 2
    nc = 3

    tf = 1
    h = tf / nh   # Final time

    rho = [
        0.11270166537926,
        0.50000000000000,
        0.88729833462074,
    ]
    bc = [1.0, 0.0]  # Boundary conditions for x
    alpha = 0.0      # Smoothing parameter
    rho_index = [(i, rho[i]) for i in 1:nc]
    j_fac = [(i, factorial(i)) for i in 1:nc]
    j_fac_1 = [(i, factorial(i-1)) for i in 1:nc]
    c = ExaModels.ExaCore(T; backend = backend)
    u = ExaModels.variable(c, nh, nc; lvar = 0, uvar = 1, start = zeros(nh, nc))
    v = ExaModels.variable(c, nh, ne; start = [mod(j, ne) for i in 1:nh, j in 1:ne])
    w = ExaModels.variable(c, nh, nc, ne; start = zeros(nh, nc, ne))
    pp = ExaModels.variable(c, nh, nc, ne; start = [mod(k, ne) for i in 1:nh, j in 1:nc, k in 1:ne])
    Dpp = ExaModels.variable(c, nh, nc, ne; start = zeros(nh, nc, ne))
    ppf = ExaModels.variable(c, ne; start = [mod(i,ne) for i in 1:ne])
    
    ExaModels.objective(c, -1.0 + ppf[1] + ppf[2])
    ExaModels.objective(c, alpha/h*(u[i+1, j] - u[i, j])^2 for i in 1:nh-1, j in 1:nc)
    
    c1 = ExaModels.constraint(
        c,
        pp[i, k, s] - v[i, s] for i=1:nh, k=1:nc, s=1:ne
    )

    ExaModels.constraint!(
        c,
        c1,
        (i-1)*nc*ne+(k-1)*ne+s => -h*(w[i, j, s]*(rho^j/j_fac)) for i in 1:nh, (k, rho) in rho_index, s=1:ne, (j, j_fac) in j_fac
    )

    c2 = ExaModels.constraint(
        c,
        Dpp[i, k, s] for i=1:nh, k=1:nc, s=1:ne
    )

    ExaModels.constraint!(
        c,
        c2,
        (i-1)*nc*ne+(k-1)*ne+s => -w[i, j, s]*(rho^(j-1)/j_fac) for i=1:nh, (k, rho) in rho_index, s=1:ne, (j, j_fac) in j_fac_1
    )

    c3 = ExaModels.constraint(
        c,
        ppf[s] - v[nh, s] for s in 1:ne
    )

    ExaModels.constraint!(
        c,
        c3,
        s => -h * w[nh, j, s] / j_fac for  s=1:ne, (j, j_fac) in j_fac
    )

    c4 = ExaModels.constraint(
        c,
        v[i, s] - v[i+1, s] for i in 1:nh-1, s in 1:ne
    )

    ExaModels.constraint!(
        c,
        c4,
        (i-1)*ne+s => w[i, j, s] * h / j_f for i in 1:nh-1, (j, j_f) in j_fac,s in 1:ne
    )

    ExaModels.constraint(
        c,
        Dpp[i,j,1] - u[i,j] * (10.0*pp[i,j,2] - pp[i,j,1]) for i=1:nh, j=1:nc
    )

    ExaModels.constraint(
        c, 
        Dpp[i,j,2] - u[i,j] * (pp[i,j,1] - 10.0*pp[i,j,2]) + (1 - u[i,j])*pp[i,j,2] for i=1:nh, j=1:nc
    )

    ExaModels.constraint(
        c,
        v[1, s] - bc for (s, bc) in [(i, bc[i]) for i in 1:ne]
    )
    
    ExaModels.ExaModel(c; kwargs...)
end


using NLPModelsIpopt, ExaModels

