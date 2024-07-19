
include("shapes/circle.jl")
include("shapes/circle_rec.jl")
include("shapes/rectangle.jl")

struct PDEProblem
    a::Float64
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
    p::Vector{Float64}
end

struct PDEDiscretizationDomain
    NODES::Int
    ELEM::Int
    DIMEN::Int
    BREAK::Int
    AREA::Vector{Float64}
    TRIANG::Matrix{Int}
    COORDS::Matrix{Float64}
    BNDRY::Vector{Int}
    EDGE::Array{Float64, 3}
    US::Vector{Float64}
    UE::Vector{Float64}
end

function PDEDiscretizationDomain(nh, domain::Dict)
    NODES = domain[:NODES]
    ELEM = domain[:ELEMS]
    COORDS = domain[:COORDS]
    BNDRY = domain[:BNDRY]
    DIMEN = 2
    BREAK = nh


    # Description of triangular elements
    TRIANG = domain[:TRIANG]
    # Edge lengths
    EDGE = [
        COORDS[TRIANG[e, mod(d1, DIMEN+1)+1], d2] - COORDS[TRIANG[e, d1], d2]
        for e in 1:ELEM, d1 in 1:DIMEN+1, d2 in 1:DIMEN
    ]
    # Area of element
    AREA = [(EDGE[e, 1, 1]*EDGE[e, 2, 2] - EDGE[e, 1, 2]*EDGE[e, 2, 1]) / 2.0 for e in 1:ELEM]
    US = domain[:US] # starting point
    UE = domain[:UE] # ending point

    return PDEDiscretizationDomain(
        NODES, ELEM, DIMEN, BREAK,
        AREA, TRIANG, COORDS, BNDRY, EDGE, US, UE,
    )
end

function _update_values!(integral, energy, u, problem::PDEProblem, dom::PDEDiscretizationDomain)
    # Unpack values
    a, b, c, d, p = problem.a, problem.b, problem.c, problem.d, problem.p
    TRIANG, EDGE = dom.TRIANG, dom.EDGE

    for b1 in 1:dom.BREAK+2, e1 in 1:dom.ELEM
        integral[b1, e1] =
            dom.AREA[e1]*(
            1 / (dom.DIMEN+1) *
                (sum((b[TRIANG[e1,c1]]*u[b1,TRIANG[e1,c1]]^2/2-
                            c[TRIANG[e1,c1]]*u[b1,TRIANG[e1,c1]]^(p[TRIANG[e1,c1]]+1)/(p[TRIANG[e1,c1]]+1)+
                            d[TRIANG[e1,c1]]*u[b1,TRIANG[e1,c1]]) for c1 in 1:dom.DIMEN+1)) +
            a / (8*dom.AREA[e1]^2)*(
                u[b1,TRIANG[e1,1]]^2*(EDGE[e1,2,1]^2 + EDGE[e1,2,2]^2) +
                u[b1,TRIANG[e1,2]]^2*(EDGE[e1,3,1]^2 + EDGE[e1,3,2]^2) +
                u[b1,TRIANG[e1,3]]^2*(EDGE[e1,1,1]^2 + EDGE[e1,1,2]^2) +
                2*u[b1,TRIANG[e1,1]]*u[b1,TRIANG[e1,2]]*(EDGE[e1,2,1]*EDGE[e1,3,1] + EDGE[e1,2,2]*EDGE[e1,3,2]) +
                2*u[b1,TRIANG[e1,1]]*u[b1,TRIANG[e1,3]]*(EDGE[e1,2,1]*EDGE[e1,1,1] + EDGE[e1,2,2]*EDGE[e1,1,2]) +
                2*u[b1,TRIANG[e1,2]]*u[b1,TRIANG[e1,3]]*(EDGE[e1,1,1]*EDGE[e1,3,1] + EDGE[e1,1,2]*EDGE[e1,3,2])
            )
        )
    end
    for b1 in 1:dom.BREAK+2
        energy[b1] = sum(integral[b1, e1] for e1 in 1:dom.ELEM)
    end
end

function _initial_position!(problem::PDEProblem, d::PDEDiscretizationDomain, niter)
    ubar = zeros(d.NODES)
    u0 = zeros(d.BREAK+2, d.NODES)
    energy0 = zeros(d.BREAK+2)
    integral0 = zeros(d.BREAK+2, d.ELEM)

    # Calculate an ending point on the other side of the barrier
    for i in 1:niter
        for b1 in 1:d.BREAK+2, n in 1:d.NODES
            u0[b1, n] = (1 - (b1-1)/(d.BREAK+1))*d.US[n] + ((b1-1)/(d.BREAK+1))*d.UE[n]
        end
        _update_values!(integral0, energy0, u0, problem, d)
        if energy0[end] < 0
            break
        end
        d.UE .*= 2.0
    end

    # Backtrack to the barrier to get a better representation
    ubar .= d.US
    for i in 1:niter
        for b1 in 1:d.BREAK+2, n in 1:d.NODES
            u0[b1, n] = (1 - (b1-1)/(d.BREAK+1))*ubar[n] + ((b1-1)/(d.BREAK+1))*d.UE[n]
        end
        _update_values!(integral0, energy0, u0, problem, d)
        ebar = maximum(energy0[b1] for b1 in 1:d.BREAK+2 if energy0[b1] < 0.0)
        bbar = minimum(b1 for b1 in 1:d.BREAK+2 if energy0[b1] == ebar)
        ubar .= u0[bbar-1, :]
        d.UE .= u0[bbar, :]
    end

    for b1 in 1:d.BREAK+2, n in 1:d.NODES
        u0[b1, n] = (1 - (b1-1)/(d.BREAK+1))*d.US[n] + ((b1-1)/(d.BREAK+1))*d.UE[n]
    end
    _update_values!(integral0, energy0, u0, problem, d)
    z0 = maximum(energy0[b1] for b1 in 1:d.BREAK+2)
    return (
        u=u0, energy=energy0, integral=integral0, z=z0,
    )
end

function _transition_state_model(problem, dom::PDEDiscretizationDomain; T = Float64, backend = nothing, kwargs...)
    a, b, c, d, p = problem.a, problem.b, problem.c, problem.d, problem.p
    x0 = _initial_position!(problem, dom, 10)
    array = [
        (
            e1, 
            dom.TRIANG[e1, 1], 
            dom.TRIANG[e1, 2], 
            dom.TRIANG[e1, 3], 
            b[dom.TRIANG[e1,1]], 
            c[dom.TRIANG[e1,1]], 
            d[dom.TRIANG[e1,1]], 
            p[dom.TRIANG[e1,1]],
            dom.AREA[e1],
            dom.EDGE[e1, 1, 1],
            dom.EDGE[e1, 1, 2],
            dom.EDGE[e1, 2, 1],
            dom.EDGE[e1, 2, 2],
            dom.EDGE[e1, 3, 1],
            dom.EDGE[e1, 3, 2]
        ) for e1 in 1:dom.ELEM
    ]
    #_proto_model()

    ALPHA = 2.0
    H = ALPHA / (dom.BREAK+1) * sqrt(sum((dom.US[n] - dom.UE[n])^2 for n in 1:dom.NODES))

    # Build optimization problem
    core = ExaModels.ExaCore(T; backend= backend)

    u = ExaModels.variable(core, 1:dom.BREAK+2, 1:dom.NODES; start=x0.u)
    integral = ExaModels.variable(core, 1:dom.BREAK+2, 1:dom.ELEM)
    z = ExaModels.variable(core, 1; start=x0.z)

    ExaModels.objective(core, z[1])

    c1 = ExaModels.constraint(
        core, 
        - z[1] for b1 in 2:dom.ELEM + 1;
        lcon = -Inf,
        ucon = 0.0,
    )

    ExaModels.constraint!(
    core,
    c1,
    b1 => integral[b1 + 1, e1] for b1 in 2:dom.ELEM + 1, e1 in 1:dom.NODES
    )

    c2 = ExaModels.constraint(
        core,
        dom.BREAK+1;
        lcon=-Inf,
        ucon=H^2,
    )

    ExaModels.constraint!(
        core,
        c2,
        b1 -1 => (u[b1+1, n] - u[b1, n])^2 for b1 in 1:dom.BREAK+1, n in 1:dom.NODES
    )

    ExaModels.constraint(
        core,
        AREA*(
            1 / (dom.DIMEN+1) *
                (((b*u[b1,TRIANG1]^2/2- c*u[b1,TRIANG1]^(p+1)/(p+1)+ d*u[b1, TRIANG1]))) +
            a / (8*AREA^2)*(
                u[b1,TRIANG1]^2*(EDGE_21^2 + EDGE_22^2) +
                u[b1,TRIANG2]^2*(EDGE_31^2 + EDGE_32^2) +
                u[b1,TRIANG3]^2*(EDGE_11^2 + EDGE_12^2) +
                2*u[b1,TRIANG1]*u[b1,TRIANG2]*(EDGE_21*EDGE_31 + EDGE_22*EDGE_32) +
                2*u[b1,TRIANG1]*u[b1,TRIANG3]*(EDGE_21*EDGE_11 + EDGE_22*EDGE_12) +
                2*u[b1,TRIANG2]*u[b1,TRIANG3]*(EDGE_11*EDGE_31 + EDGE_12*EDGE_32)
                )
            ) 
            - integral[b1, e1]
        for b1 in 1:dom.BREAK+2, (e1, TRIANG1, TRIANG2, TRIANG3, b, c, d, p, AREA, EDGE_11, EDGE_12, EDGE_21, EDGE_22, EDGE_31, EDGE_32) in array
    )

    # Boundary
    boundary_nodes = findall(isequal(1), dom.BNDRY)
    ExaModels.constraint(
        core,
        u[b1+1, n] for b1 in 1:dom.BREAK, n in boundary_nodes
    )
    ExaModels.constraint(
        core,
        u[1, n] for n in 1:dom.NODES;
        lcon=dom.US,
        ucon=dom.US,
    )
    ExaModels.constraint(
        core,
        u[dom.BREAK+2, n] for n in 1:dom.NODES;
        lcon=dom.UE,
        ucon=dom.UE,
    )
    return core
end

function dirichlet_model(nh)
    dom = PDEDiscretizationDomain(nh, CIRCLE_DOMAIN)
    println(dom.ELEM)
    pb = PDEProblem(
        0.01,
        fill(1.0, dom.NODES),
        fill(1.0, dom.NODES),
        fill(0.0, dom.NODES),
        fill(3.0, dom.NODES),
    )
    return _transition_state_model(pb, dom)
end

function henon_model(nh)
    dom = PDEDiscretizationDomain(nh, CIRCLE_REC_DOMAIN)
    pb = PDEProblem(
        1.0,
        fill(0.0, dom.NODES),
        sqrt.(dom.COORDS[:, 1].^2 .+ dom.COORDS[:, 2].^2),
        fill(0.0, dom.NODES),
        fill(3.0, dom.NODES),
    )
    return _transition_state_model(pb, dom)
end

function lane_emden_model(nh)
    dom = PDEDiscretizationDomain(nh, RECTANGLE_DOMAIN)
    pb = PDEProblem(
        1.0,
        fill(0.0, dom.NODES),
        fill(1.0, dom.NODES),
        fill(0.0, dom.NODES),
        fill(3.0, dom.NODES),
    )
    return _transition_state_model(pb, dom)
end


