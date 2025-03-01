# Methanol-to-Hydrocarbons Problem
# Collocation formulation
# Michael Merritt - Summer 2000
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function methanol_model(nh; T = Float64, backend = nothing, kwargs...)
    ne = 3
    np = 5
    nc = 3
    nm = 17

    rho = [0.11270166537926, 0.5, 0.88729833462074]
    # times at which observations made
    tau = [
        0,
        0.050,
        0.065,
        0.080,
        0.123,
        0.233,
        0.273,
        0.354,
        0.397,
        0.418,
        0.502,
        0.553,
        0.681,
        0.750,
        0.916,
        0.937,
        1.122,
    ]
    tf = tau[nm]                   # ODEs defined in [0,tf]
    h = tf / nh                    # uniform interval length
    t = [(i-1)*h for i in 1:nh+1]  # partition

    # itau[i] is the largest integer k with t[k] <= tau[i]
    itau = Int[min(nh, floor(tau[i]/h)+1) for i in 1:nm]

    # Concentrations
    z = [
        1.0000         0         0;
        0.7085    0.1621    0.0811;
        0.5971    0.1855    0.0965;
        0.5537    0.1989    0.1198;
        0.3684    0.2845    0.1535;
        0.1712    0.3491    0.2097;
        0.1198    0.3098    0.2628;
        0.0747    0.3576    0.2467;
        0.0529    0.3347    0.2884;
        0.0415    0.3388    0.2757;
        0.0261    0.3557    0.3167;
        0.0208    0.3483    0.2954;
        0.0085    0.3836    0.2950;
        0.0053    0.3611    0.2937;
        0.0019    0.3609    0.2831;
        0.0018    0.3485    0.2846;
        0.0006    0.3698    0.2899;
    ]
    con1_matrix = [(j, s, itau[j], tau[j], z[j,s], t[itau[j]]) for j in 1:nm, s in 1:ne]
    bc = [1.0, 0.0, 0.0]

    # Starting-value
    v0 = zeros(nh, ne)
    for i in 1:itau[1], s in 1:ne
        v0[i, s] = bc[s]
    end
    for j in 2:nm, i in itau[j-1]+1:itau[j], s in 1:ne
        v0[i, s] = z[j, s]
    end
    for i in itau[nm]+1:nh, s in 1:ne
        v0[i, s] = z[nm, s]
    end
    v0 .= 0.001

    c = ExaModels.ExaCore(T; backend = backend)
    theta = ExaModels.variable(c, np; lvar = 0, start = fill(1, np))
    v = ExaModels.variable(c, nh, ne; start = v0)
    w = ExaModels.variable(c, nh, nc, ne; start = 0)
    uc = ExaModels.variable(c, nh, nc, ne; start = [v0[i,s] for i=1:nh, j=1:nc, s=1:ne])
    Duc = ExaModels.variable(c, nh, nc, ne; start = 0)

    ExaModels.objective(c, (v[itau,s] + sum(w[itau,k,s]*(tau-t)^k/(factorial(k)*h^(k-1)) for k in 1:nc) - z)^2 for (j,s,itau,tau,z,t) in con1_matrix)

    ExaModels.constraint(
        c,
        uc[i, j, s] - v[i,s] - h*sum(w[i,k,s]*(rho^k/factorial(k)) for k in 1:nc) for i=1:nh, (j,rho) in [(j, rho[j]) for j in 1:nc], s=1:ne
    )

    ExaModels.constraint(
        c,
        Duc[i, j, s] - sum(w[i,k,s]*(rho^(k-1)/factorial(k-1)) for k in 1:nc) for i=1:nh, (j,rho) in [(j, rho[j]) for j in 1:nc], s=1:ne
    )

    ExaModels.constraint(
        c,
        v[1, s] - bc for (s, bc) in [(s, bc[s]) for s in 1:ne]

    )

    ExaModels.constraint(
        c,
        v[i, s] + sum(w[i, j, s]*h/factorial(j) for j in 1:nc) - v[i+1, s] for i=1:nh-1, s=1:ne
    )

    ExaModels.constraint(
        c,
        Duc[i,j,1] + ((2*theta[2] - (theta[1]*uc[i,j,2])/((theta[2]+theta[5])*uc[i,j,1]+uc[i,j,2]) +
                         theta[3] + theta[4])*uc[i,j,1]) for i=1:nh, j=1:nc
    )    

    ExaModels.constraint(
        c,
        Duc[i,j,2] - ((theta[1]*uc[i,j,1]*(theta[2]*uc[i,j,1]-uc[i,j,2]))/ ((theta[2]+theta[5])*uc[i,j,1]+uc[i,j,2]) +
                     theta[3]*uc[i,j,1]) for i=1:nh, j=1:nc
    )   
    
    ExaModels.constraint(
        c,
        Duc[i,j,3] - ((theta[1]*uc[i,j,1]*(uc[i,j,2]+theta[5]*uc[i,j,1]))/ ((theta[2]+theta[5])*uc[i,j,1]+uc[i,j,2]) +
                        theta[4]*uc[i,j,1]) for i=1:nh, j=1:nc
    )   

    return ExaModels.ExaModel(c; kwargs...)
end
