convert_data(data::N, backend) where {names,N<:NamedTuple{names}} =
    NamedTuple{names}(ExaModels.convert_array(d, backend) for d in data)
parse_ac_power_data(filename, backend) =
    convert_data(parse_ac_power_data(filename), backend)


function parse_ac_power_data(filename)
    d, f = splitdir(filename)
    name, ext = splitext(f)

    if isfile(joinpath(TMPDIR, name) * ".jld2")
        @info "Loading cached JLD2 file"
        return JLD2.load(joinpath(TMPDIR, name) * ".jld2", "data")
    else
        ff = if isfile(filename)
            filename
        elseif isfile(joinpath(TMPDIR, name) * ".m")
            joinpath(TMPDIR, name) * ".m"
        else
            @info "Downloading $filename"
            Downloads.download(
                "https://raw.githubusercontent.com/power-grid-lib/pglib-opf/dc6be4b2f85ca0e776952ec22cbd4c22396ea5a3/$filename",
                joinpath(TMPDIR, name * ".m"),
            )
            joinpath(TMPDIR, name * ".m")
        end
        @info "Loading MATPOWER file"
        return process_ac_power_data(ff)
    end
end

function process_ac_power_data(filename)
    data = PowerModels.parse_file(filename)
    PowerModels.standardize_cost_terms!(data, order = 2)
    PowerModels.calc_thermal_limits!(data)

    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    arcdict = Dict(a => k for (k, a) in enumerate(ref[:arcs]))
    busdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:bus]))
    gendict = Dict(k => i for (i, (k, v)) in enumerate(ref[:gen]))
    branchdict = Dict(k => i for (i, (k, v)) in enumerate(ref[:branch]))

    data =  (
        bus = [
            begin
                bus_loads = [ref[:load][l] for l in ref[:bus_loads][k]]
                bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][k]]
                pd = sum(load["pd"] for load in bus_loads; init = 0.0)
                gs = sum(shunt["gs"] for shunt in bus_shunts; init = 0.0)
                qd = sum(load["qd"] for load in bus_loads; init = 0.0)
                bs = sum(shunt["bs"] for shunt in bus_shunts; init = 0.0)
                (i = busdict[k], pd = pd, gs = gs, qd = qd, bs = bs, bus_type = v["bus_type"])
            end for (k, v) in ref[:bus]
        ],
        gen = [
            (
                i = gendict[k],
                cost1 = v["cost"][1],
                cost2 = v["cost"][2],
                cost3 = v["cost"][3],
                bus = busdict[v["gen_bus"]],
            ) for (k, v) in ref[:gen]
        ],
        arc = [
            (i = k, rate_a = ref[:branch][l]["rate_a"], bus = busdict[i]) for
            (k, (l, i, j)) in enumerate(ref[:arcs])
                ],
        branch = [
            begin
                f_idx = arcdict[i, branch["f_bus"], branch["t_bus"]]
                t_idx = arcdict[i, branch["t_bus"], branch["f_bus"]]
                g, b = PowerModels.calc_branch_y(branch)
                tr, ti = PowerModels.calc_branch_t(branch)
                ttm = tr^2 + ti^2
                g_fr = branch["g_fr"]
                b_fr = branch["b_fr"]
                g_to = branch["g_to"]
                b_to = branch["b_to"]
                c1 = (-g * tr - b * ti) / ttm
                c2 = (-b * tr + g * ti) / ttm
                c3 = (-g * tr + b * ti) / ttm
                c4 = (-b * tr - g * ti) / ttm
                c5 = (g + g_fr) / ttm
                c6 = (b + b_fr) / ttm
                c7 = (g + g_to)
                c8 = (b + b_to)
                (
                    i = branchdict[i],
                    j = 1,
                    f_idx = f_idx,
                    t_idx = t_idx,
                    f_bus = busdict[branch["f_bus"]],
                    t_bus = busdict[branch["t_bus"]],
                    c1 = c1,
                    c2 = c2,
                    c3 = c3,
                    c4 = c4,
                    c5 = c5,
                    c6 = c6,
                    c7 = c7,
                    c8 = c8,
                    rate_a_sq = branch["rate_a"]^2,
                )
            end for (i, branch) in ref[:branch]
                ],
        ref_buses = [busdict[i] for (i, k) in ref[:ref_buses]],
        vmax = [v["vmax"] for (k, v) in ref[:bus]],
        vmin = [v["vmin"] for (k, v) in ref[:bus]],
        pmax = [v["pmax"] for (k, v) in ref[:gen]],
        pmin = [v["pmin"] for (k, v) in ref[:gen]],
        qmax = [v["qmax"] for (k, v) in ref[:gen]],
        qmin = [v["qmin"] for (k, v) in ref[:gen]],
        rate_a = [ref[:branch][l]["rate_a"] for (k, (l, i, j)) in enumerate(ref[:arcs])],
        angmax = [b["angmax"] for (i, b) in ref[:branch]],
        angmin = [b["angmin"] for (i, b) in ref[:branch]],
    )

    @info "Saving JLD2 cache file"
    d, f = splitdir(filename)
    name,ext = splitext(f)
    JLD2.save(joinpath(TMPDIR, name * ".jld2"), "data", data)
    return data
end

function ac_power_model(
    filename = "pglib_opf_case3_lmbd.m";
    backend = nothing,
    T = Float64,
    kwargs...,
)

    data = parse_ac_power_data(filename, backend)

    w = ExaModels.ExaCore(T; backend = backend)

    va = ExaModels.variable(w, length(data.bus);)

    vm = ExaModels.variable(
        w,
        length(data.bus);
        start = fill!(similar(data.bus, Float64), 1.0),
        lvar = data.vmin,
        uvar = data.vmax,
    )
    pg = ExaModels.variable(w, length(data.gen); lvar = data.pmin, uvar = data.pmax)

    qg = ExaModels.variable(w, length(data.gen); lvar = data.qmin, uvar = data.qmax)

    p = ExaModels.variable(w, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    q = ExaModels.variable(w, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    o = ExaModels.objective(
        w,
        g.cost1 * pg[g.i]^2 + g.cost2 * pg[g.i] + g.cost3 for g in data.gen
    )

    c1 = ExaModels.constraint(w, va[i] for i in data.ref_buses)

    c2 = ExaModels.constraint(
        w,
        p[b.f_idx] - b.c5 * vm[b.f_bus]^2 -
        b.c3 * (vm[b.f_bus] * vm[b.t_bus] * cos(va[b.f_bus] - va[b.t_bus])) -
        b.c4 * (vm[b.f_bus] * vm[b.t_bus] * sin(va[b.f_bus] - va[b.t_bus])) for
        b in data.branch
    )

    c3 = ExaModels.constraint(
        w,
        q[b.f_idx] +
        b.c6 * vm[b.f_bus]^2 +
        b.c4 * (vm[b.f_bus] * vm[b.t_bus] * cos(va[b.f_bus] - va[b.t_bus])) -
        b.c3 * (vm[b.f_bus] * vm[b.t_bus] * sin(va[b.f_bus] - va[b.t_bus])) for
        b in data.branch
    )

    c4 = ExaModels.constraint(
        w,
        p[b.t_idx] - b.c7 * vm[b.t_bus]^2 -
        b.c1 * (vm[b.t_bus] * vm[b.f_bus] * cos(va[b.t_bus] - va[b.f_bus])) -
        b.c2 * (vm[b.t_bus] * vm[b.f_bus] * sin(va[b.t_bus] - va[b.f_bus])) for
        b in data.branch
    )

    c5 = ExaModels.constraint(
        w,
        q[b.t_idx] +
        b.c8 * vm[b.t_bus]^2 +
        b.c2 * (vm[b.t_bus] * vm[b.f_bus] * cos(va[b.t_bus] - va[b.f_bus])) -
        b.c1 * (vm[b.t_bus] * vm[b.f_bus] * sin(va[b.t_bus] - va[b.f_bus])) for
        b in data.branch
    )

    # c6 = ExaModels.constraint(
    #     w,
    #     va[b.f_bus] - va[b.t_bus] for b in data.branch;
    #     lcon = data.angmin,
    #     ucon = data.angmax,
    # )
    c7 = ExaModels.constraint(
        w,
        p[b.f_idx]^2 + q[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )
    c8 = ExaModels.constraint(
        w,
        p[b.t_idx]^2 + q[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    c9 = ExaModels.constraint(w, b.pd + b.gs * vm[b.i]^2 for b in data.bus)

    c10 = ExaModels.constraint(w, b.qd - b.bs * vm[b.i]^2 for b in data.bus)

    c11 = ExaModels.constraint!(w, c9, a.bus => p[a.i] for a in data.arc)
    c12 = ExaModels.constraint!(w, c10, a.bus => q[a.i] for a in data.arc)

    c13 = ExaModels.constraint!(w, c9, g.bus => -pg[g.i] for g in data.gen)
    c14 = ExaModels.constraint!(w, c10, g.bus => -qg[g.i] for g in data.gen)

    return ExaModels.ExaModel(w; kwargs...)

end


function ac_rectangle_power_model(
    filename = "pglib_opf_case3_lmbd.m";
    backend = nothing,
    T = Float64,
    kwargs...,
)

    data = parse_ac_power_data(filename, backend)

    w = ExaModels.ExaCore(T; backend = backend)

    # real parts of voltage 
    vr = ExaModels.variable(w, length(data.bus); start = fill(1, length(data.bus)))

    # imaginary parts of voltage
    vi = ExaModels.variable(w, length(data.bus);)

    # active power generation
    pg = ExaModels.variable(w, length(data.gen); lvar = data.pmin, uvar = data.pmax)

    # reactive power generation
    qg = ExaModels.variable(w, length(data.gen); lvar = data.qmin, uvar = data.qmax)

    # active power in arc
    p = ExaModels.variable(w, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    # reactive power in arc
    q = ExaModels.variable(w, length(data.arc); lvar = -data.rate_a, uvar = data.rate_a)

    ExaModels.objective(
        w,
        g.cost1 * pg[g.i]^2 + g.cost2 * pg[g.i] + g.cost3 for g in data.gen
    )

    # reference bus
    ExaModels.constraint(
        w,
        vi[i] for i in data.ref_buses
    )

    #  vmin <= v[i] <= vmax
    ExaModels.constraint(
        w,
        vr[i]^2 + vi[i]^2 for i in 1:length(data.bus);
        lcon = data.vmin,
        ucon = data.vmax
    )

    # power balance
    p_bal = ExaModels.constraint(w, b.pd + b.gs * (vi[b.i]^2 + vr[b.i]^2) for b in data.bus)

    q_bal = ExaModels.constraint(w, b.qd - b.bs * (vi[b.i]^2 + vr[b.i]^2) for b in data.bus)

    ExaModels.constraint!(w, p_bal, a.bus => p[a.i] for a in data.arc)

    ExaModels.constraint!(w, q_bal, a.bus => q[a.i] for a in data.arc)

    ExaModels.constraint!(w, p_bal, g.bus => -pg[g.i] for g in data.gen)

    ExaModels.constraint!(w, q_bal, g.bus => -qg[g.i] for g in data.gen)

    # ohms_yt_form
    ExaModels.constraint(
        w,
        p[b.f_idx] - b.c5 * (vi[b.f_bus]^2 + vr[b.f_bus]^2) -
        b.c3 * (vi[b.f_bus] * vi[b.t_bus] + vr[b.f_bus] * vr[b.t_bus]) -
        b.c4 * (vi[b.f_bus] * vr[b.t_bus] - vr[b.f_bus] * vi[b.t_bus]) for
        b in data.branch
    )

    ExaModels.constraint(
        w,
        q[b.f_idx] + b.c6 * (vi[b.f_bus]^2 + vr[b.f_bus]^2) +
        b.c4 * (vi[b.f_bus] * vi[b.t_bus] + vr[b.f_bus] * vr[b.t_bus]) -
        b.c3 * (vi[b.f_bus] * vr[b.t_bus] - vr[b.f_bus] * vi[b.t_bus]) for
        b in data.branch
    )

    # ohms_yt_to
    ExaModels.constraint(
        w,
        p[b.t_idx] - b.c7 * (vi[b.t_bus]^2 + vr[b.t_bus]^2) -
        b.c1 * (vi[b.f_bus] * vi[b.t_bus] + vr[b.f_bus] * vr[b.t_bus]) -
        b.c2 * (-(vi[b.f_bus] * vr[b.t_bus] - vr[b.f_bus] * vi[b.t_bus])) for
        b in data.branch
    )

    ExaModels.constraint(
        w,
        q[b.t_idx] + b.c8 * (vi[b.t_bus]^2 + vr[b.t_bus]^2) +
        b.c2 * (vi[b.f_bus] * vi[b.t_bus] + vr[b.f_bus] * vr[b.t_bus]) -
        b.c1 * (-(vi[b.f_bus] * vr[b.t_bus] - vr[b.f_bus] * vi[b.t_bus])) for
        b in data.branch
    )

    ExaModels.constraint(
        w,
        (vi[b.f_bus] * vr[b.t_bus] - vr[b.f_bus] * vi[b.t_bus])/ (vr[b.f_bus] * vr[b.t_bus] + vi[b.f_bus] * vi[b.t_bus]) for b in data.branch;
        lcon = tan.(data.angmin),
        ucon = tan.(data.angmax)
    )

    ExaModels.constraint(
        w,
        p[b.f_idx]^2 + q[b.f_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )
    ExaModels.constraint(
        w,
        p[b.t_idx]^2 + q[b.t_idx]^2 - b.rate_a_sq for b in data.branch;
        lcon = fill!(similar(data.branch, Float64, length(data.branch)), -Inf),
    )

    return ExaModels.ExaModel(w; kwargs...)
end
