# Flow in a Channel Problem
# Collocation formulation
# Alexander S. Bondarenko - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function channel_model(nh; T = Float64, backend = nothing, kwargs...)
    nc = 4
    nd = 4
    R = 10.0 # Reynolds number
    tf = 1.0
    h = tf / nh

    bc = [0.0 1.0; 0.0 0.0]
    rho = [0.06943184420297, 0.33000947820757, 0.66999052179243, 0.93056815579703]
    t = [(i-1)*h for i in 1:nh+1]
    tuples = Vector{Tuple{Int, Int, Int}}()
    for s in 1:nd
        for j in s:nd
            push!(tuples, (s, j, factorial(j-s)))
        end
    end    
    array = [(i,j, factorial(j+nd-i)) for i in 1:nd, j in 1:nc]
    rho_j = [(j, rho[j]) for j in 1:nc]
    # Initial value
    v0 = zeros(nh, nd)
    for i in 1:nh
        v0[i, 1] = t[i]^2*(3.0 - 2.0*t[i])
        v0[i, 2] = 6*t[i]*(1.0 - t[i])
        v0[i, 3] = 6*(1.0 - 2.0*t[i])
        v0[i, 4] = -12.0
    end

    core = ExaModels.ExaCore(T; backend= backend)

    v = ExaModels.variable(core, nh, nd;)
    w = ExaModels.variable(core, nh, nc; start=0.0)

    uc = ExaModels.variable(core, nh, nc, nd; start=[v0[i, s] for i=1:nh, j=1:nc, s=1:nd])
    Duc = ExaModels.variable(core, nh, nc, nd; start=0.0)
    y = ExaModels.variable(core, 1; lvar =1, uvar =1, start = 1)
    
    # Constant objective
    ExaModels.objective(core, y[1])

    # Collocation model
    ExaModels.constraint(
        core,
        uc[i, j, s] - v[i,s] - h*sum(w[i,k]*(rho^k/factorial(k)) for k in 1:nc) for i in 1:nh, (j,rho) in rho_j, s in 1:nd
    )
    c1 = ExaModels.constraint(
        core,
        Duc[i, j, s] for i in 1:nh, j in 1:nc, s in 1:nd
    )

    ExaModels.constraint!(
        core,
        c1,
        (s-1)*nc*nh + (j-1)*nh + i=> - v[i,k]*((rho*h)^(k-s))/fac_k for i in 1:nh, (j,rho) in rho_j, (s, k, fac_k) in tuples
    )

    ExaModels.constraint!(
        core,
        c1,
        (s-1)*nc*nh + (j-1)*nh + i => - h^(nd-s+1) * w[i, k]*(rho^(k+nd-s)/fac_nd) for i in 1:nh, (j,rho) in rho_j, (s, k, fac_nd) in array
    )

    # Boundary
    ExaModels.constraint(core, v[1, 1] - bc[1, 1])
    ExaModels.constraint(core, v[1, 2] -  bc[2, 1])
    ExaModels.constraint(
        core,
        sum(v[nh, k]*(h^(k-1)/factorial(k-1)) for k in 1:nd) +
        h^nd * sum(w[nh, k]/factorial(k+nd-1) for k in 1:nc) - bc[1, 2]
    )
    ExaModels.constraint(
        core,
        sum(v[nh, k]*(h^(k-2)/factorial(k-2)) for k in 2:nd) +
        h^(nd-1) * sum(w[nh, k]/factorial(k+nd-2) for k in 1:nc) - bc[2, 2]
    )
    c2 = ExaModels.constraint(
        core,
        - v[i+1, s] for i in 1:nh-1, s in 1:nd
    )
    ExaModels.constraint!(
        core,
        c2,
        (s-1)*(nh-1) + i => v[i, k]*(h^(k-s)/fac_k) for i in 1:nh-1, (s, k, fac_k) in tuples
    )
    ExaModels.constraint!(
        core,
        c2,
        (s-1)*(nh-1) + i =>  h^(nd-s+1)* w[i, k]/fac_nd for i in 1:nh-1, (s, k, fac_nd) in array

    )
    ExaModels.constraint(
        core,
        sum(w[i, k] * (rho^(k-1)/factorial(k-1)) for k in 1:nc) -
        R * (Duc[i, j, 2] * Duc[i, j, 3] - Duc[i, j, 1] * Duc[i, j, 4]) for i in 1:nh, (j,rho) in rho_j
    )

    return ExaModels.ExaModel(core; kwargs...)
end
