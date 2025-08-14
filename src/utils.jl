function project!(l, x, u; marg = 1e-4)
    map!(x, l, x, u) do l, x, u
        max(l + marg, min(u - marg, x))
    end
end

function btime(f, N)
    f()
    GC.gc()
    return minimum(begin
        t = time_ns()
        f()
        (time_ns() - t) / 1e9
    end for i = 1:N)
end

function benchmark_callbacks(m; N = 100)
    nvar = m.meta.nvar
    ncon = m.meta.ncon
    nnzj = m.meta.nnzj
    nnzh = m.meta.nnzh

    x = copy(m.meta.x0)
    y = similar(m.meta.x0, ncon)
    c = similar(m.meta.x0, ncon)
    g = similar(m.meta.x0, nvar)
    jac = similar(m.meta.x0, nnzj)
    hess = similar(m.meta.x0, nnzh)
    jrows = similar(m.meta.x0, Int, nnzj)
    jcols = similar(m.meta.x0, Int, nnzj)
    hrows = similar(m.meta.x0, Int, nnzh)
    hcols = similar(m.meta.x0, Int, nnzh)

    project!(m.meta.lvar, x, m.meta.uvar)


    tobj = btime(() -> ExaModels.obj(m, x), N)

    tcon = btime(() -> ExaModels.cons!(m, x, c), N)

    tgrad = btime(() -> ExaModels.grad!(m, x, g), N)


    tjac = btime(() -> ExaModels.jac_coord!(m, x, jac), N)

    thess = btime(() -> ExaModels.hess_coord!(m, x, y, hess), N)

    tjacs = btime(() -> ExaModels.jac_structure!(m, jrows, jcols), N)

    thesss = btime(() -> ExaModels.hess_structure!(m, hrows, hcols), N)

    return (
        tobj = tobj,
        tcon = tcon,
        tgrad = tgrad,
        tjac = tjac,
        thess = thess,
        tjacs = tjacs,
        thesss = thesss,
    )
end
