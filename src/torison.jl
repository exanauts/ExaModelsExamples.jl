# Torsion problem
# Liz Dolan - Summer 2000
# Version 2.0 - October 2000
# COPS 3.1 - March 2004

function torison_model(nx, ny; T = Float64, backend = nothing, kwargs...)
    c1 = 5.0
    hx = 1.0 / (nx + 1.0)    # grid spacing
    hy = 1.0 / (ny + 1.0)    # grid spacing
    area = 0.5 * hx * hy     # area of triangle

    # Distance to the boundary.
    D = [min(min(i,nx-i+1)*hx, min(j, ny-j+1)*hy) for i in 0:nx+1, j in 0:ny+1]

    c = ExaModels.ExaCore(T; backend = backend)
    v = ExaModels.variable(c, nx+2, ny+2; start = D)

    ExaModels.constraint(
        c,
        v[i,j] for i in 1:nx+2, j in 1:ny+2;
        lcon = -D,
        ucon = D
    )

    ExaModels.objective(
        c,
        area * (((v[i+1,j] - v[i,j])/hx)^2 + ((v[i,j+1] - v[i,j])/hy)^2) / 2  for i in 1:nx+1, j in 1:ny +1
    )

    ExaModels.objective(
        c,
        area * (((v[i,j] - v[i-1,j])/hx)^2 + ((v[i,j] - v[i,j-1])/hy)^2) / 2 for i in 2:nx+2, j in 2:ny+2 
    )

    ExaModels.objective(
        c,
        - area * c1* (v[i+1, j] + v[i, j] + v[i, j+1]) /3 for i in 1:nx+1, j in 1:ny+1
    )

    ExaModels.objective(
        c,
        - area * c1 * (v[i, j] + v[i-1, j] + v[i, j-1]) / 3 for i in 2:nx+2, j in 2:ny+2
    )

    return ExaModels.ExaModel(c; kwargs...)
end

