# Minimal surface with obstacle problem

#  Find the surface with minimal area, given boundary conditions, 
#  and above an obstacle.

#  This is problem 17=the COPS (Version 3) collection of 
#  E. Dolan and J. More'
#  see "Benchmarking Optimization Software with COPS"
#  Argonne National Labs Technical Report ANL/MCS-246 (2004)
#  classification OBR2-AN-V-V

function minsurf_model(nx::Int, ny::Int;T = Float64, backend = nothing, kwargs...)
    x_mesh = LinRange(0, 1, nx + 2) # coordinates of the mesh points x

    v0 = zeros(nx + 2, ny + 2) # Surface matrix initialization
    for i = 1:(nx + 2), j = 1:(ny + 2)
        v0[i, j] = 1 - (2 * x_mesh[i] - 1)^2
    end

    hx = 1 / (nx + 1)
    hy = 1 / (ny + 1)
    area = 1 // 2 * hx * hy

    c = ExaMOdels.ExaCore(T; backend = backend)
    v = ExaModels.variable(c, nx+2, ny+2; start = v0)

    ExaModels.objective(c, area * (1 + ((v[i + 1, j] - v[i, j]) / hx)^2 + ((v[i, j + 1] - v[i, j]) / hy)^2)^(1 / 2) for
                        i = 1:(nx + 1), j = 1:(ny + 1))
    ExaModels.objective(c, area * (1 + ((v[i - 1, j] - v[i, j]) / hx)^2 + ((v[i, j - 1] - v[i, j]) / hy)^2)^(1 / 2) for
                         i = 2:(nx + 2), j = 2:(ny + 2))    
    
    ExaModels.constraint(
        c,
        v[1, j + 1] for j in 0:ny+1
    )

    ExaModels.constraint(
        c,
        v[nx + 2, j + 1] for j in 0:ny+1
    )

    ExaModels.constraint(
        c,
        v[i + 1, 1] - 1 + (2 * i * hx - 1)^2 for i in 0:nx+1
    )

    ExaModels.constraint(
        c,
        v[i + 1, 1] - 1 + (2 * i * hx - 1)^2 for i in 0:nx+1
    )    

    ExaModels.constraint(
        c,
        v[i + 1, j + 1] for j in 0:ny+1, i in 0:nx+1;
        lcon = 0,
        ucon = Inf
    )

    ExaModels.constraint(
        c,
        v[i + 1, j + 1] for i in Int(floor(0.25 / hx)):Int(ceil(0.75 / hx)), j in Int(floor(0.25 / hy)):Int(ceil(0.75 / hy));
        lcon = 1,
        ucon = Inf
    )
    
    return ExaModels.ExaModel(core; kwargs...)
end


