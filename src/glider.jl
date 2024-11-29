# Hang Glider Problem
# Trapezoidal formulation
# David Bortz - Summer 1998
# COPS 2.0 - September 2000
# COPS 3.0 - November 2002
# COPS 3.1 - March 2004

function glider_model(nh; T = Float64, backend = nothing, kwargs...)
    # Design parameters
    x_0 = 0.0
    y_0 = 1000.0
    y_f = 900.0
    vx_0 = 13.23
    vx_f = 13.23
    vy_0 = -1.288
    vy_f = -1.288
    u_c = 2.5
    r_0 = 100.0
    m = 100.0
    g = 9.81
    c0 = 0.034
    c1 = 0.069662
    S = 14.0
    rho = 1.13
    cL_min = 0.0
    cL_max = 1.4

    c = ExaModels.ExaCore(T; backend = backend)
    t_f = ExaModels.variable(c, 1; lvar = 0, start = 1)
    x = ExaModels.variable(c, nh+1; lvar = zeros(nh+1), start = [x_0 + vx_0*(k/nh) for k in 0:nh])
    y = ExaModels.variable(c, nh+1; start = [y_0 + (k/nh)*(y_f - y_0) for k in 0:nh])
    vx = ExaModels.variable(c, nh+1; lvar = zeros(nh+1), start = fill(vx_0, nh+1))
    vy = ExaModels.variable(c, nh+1; start = fill(vy_0, nh+1))
    cL = ExaModels.variable(c, nh+1; lvar = fill(cL_min,nh+1), uvar = fill(cL_max, nh+1), start = fill(cL_max/2, nh+1))

    ExaModels.objective(c, -x[nh+1])

    ExaModels.constraint(
        c,
        x[j] - (x[j-1] + 0.5/nh * (vx[j] + vx[j-1]) * t_f[1]) for j in 2:nh+1
    )

    ExaModels.constraint(
        c,
        y[j] - (y[j-1] + 0.5 * t_f[1]/nh * (vy[j] + vy[j-1])) for j in 2:nh+1
    )

    ExaModels.constraint(
        c,
        vx[i] - (vx[i-1] + 
        0.5 * t_f[1]/nh * ((-(0.5*cL[i]*rho*S*(sqrt(vx[i]^2 + (vy[i] - u_c*(1 - (x[i]/r_0 - 2.5)^2)*exp(-(x[i]/r_0 - 2.5)^2))^2))^2)*((vy[i] - u_c*(1 - (x[i]/r_0 - 2.5)^2)*exp(-(x[i]/r_0 - 2.5)^2))/( sqrt(vx[i]^2 + (vy[i] - u_c*(1 - (x[i]/r_0 - 2.5)^2)*exp(-(x[i]/r_0 - 2.5)^2))^2))) - (0.5*(c0+c1*cL[i]^2)*rho*S*(sqrt(vx[i]^2 + (vy[i] - u_c*(1 - (x[i]/r_0 - 2.5)^2)*exp(-(x[i]/r_0 - 2.5)^2))^2))^2)*(vx[i]/(sqrt(vx[i]^2 + (vy[i] - u_c*(1 - (x[i]/r_0 - 2.5)^2)*exp(-(x[i]/r_0 - 2.5)^2))^2))))/m + 
        (-(0.5*cL[i-1]*rho*S*(sqrt(vx[i-1]^2 + (vy[i-1] - u_c*(1 - (x[i-1]/r_0 - 2.5)^2)*exp(-(x[i-1]/r_0 - 2.5)^2))^2))^2)*((vy[i-1] - u_c*(1 - (x[i-1]/r_0 - 2.5)^2)*exp(-(x[i-1]/r_0 - 2.5)^2))/( sqrt(vx[i-1]^2 + (vy[i-1] - u_c*(1 - (x[i-1]/r_0 - 2.5)^2)*exp(-(x[i-1]/r_0 - 2.5)^2))^2))) - (0.5*(c0+c1*cL[i-1]^2)*rho*S*(sqrt(vx[i-1]^2 + (vy[i-1] - u_c*(1 - (x[i-1]/r_0 - 2.5)^2)*exp(-(x[i-1]/r_0 - 2.5)^2))^2))^2)*(vx[i-1]/(sqrt(vx[i-1]^2 + (vy[i-1] - u_c*(1 - (x[i-1]/r_0 - 2.5)^2)*exp(-(x[i-1]/r_0 - 2.5)^2))^2))))/m)) for i in 2:nh+1
    )
    #The following two lengthy expressions represent vx_dot and vy_dot respectively. This is due to ExaModels does not support @expression.
    # @variables(model, begin
    # 0 <= t_f,                       (start=1.0)
    # 0.0 <= x[k=0:nh],               (start=x_0 + vx_0*(k/nh))
    # y[k=0:nh],                      (start=y_0 + (k/nh)*(y_f - y_0))
    # 0.0 <= vx[k=0:nh],              (start=vx_0)
    # vy[k=0:nh],                     (start=vy_0)
    # cL_min <= cL[k=0:nh] <= cL_max, (start=cL_max/2.0)
    # end)

    ExaModels.constraint(
        c,
        vy[i] - (vy[i-1] + 
        0.5 * t_f[1]/nh * (((0.5*cL[i]*rho*S*(sqrt(vx[i]^2 + (vy[i] - u_c*(1 - (x[i]/r_0 - 2.5)^2)*exp(-(x[i]/r_0 - 2.5)^2))^2))^2)*(vx[i]/(sqrt(vx[i]^2 + (vy[i] - u_c*(1 - (x[i]/r_0 - 2.5)^2)*exp(-(x[i]/r_0 - 2.5)^2))^2))) - (0.5*(c0+c1*cL[i]^2)*rho*S*(sqrt(vx[i]^2 + (vy[i] - u_c*(1 - (x[i]/r_0 - 2.5)^2)*exp(-(x[i]/r_0 - 2.5)^2))^2))^2)*((vy[i] - u_c*(1 - (x[i]/r_0 - 2.5)^2)*exp(-(x[i]/r_0 - 2.5)^2))/(sqrt(vx[i]^2 + (vy[i] - u_c*(1 - (x[i]/r_0 - 2.5)^2)*exp(-(x[i]/r_0 - 2.5)^2))^2))))/m - g + 
        ((0.5*cL[i-1]*rho*S*(sqrt(vx[i-1]^2 + (vy[i-1] - u_c*(1 - (x[i-1]/r_0 - 2.5)^2)*exp(-(x[i-1]/r_0 - 2.5)^2))^2))^2)*(vx[i-1]/(sqrt(vx[i-1]^2 + (vy[i-1] - u_c*(1 - (x[i-1]/r_0 - 2.5)^2)*exp(-(x[i-1]/r_0 - 2.5)^2))^2))) - (0.5*(c0+c1*cL[i-1]^2)*rho*S*(sqrt(vx[i-1]^2 + (vy[i-1] - u_c*(1 - (x[i-1]/r_0 - 2.5)^2)*exp(-(x[i-1]/r_0 - 2.5)^2))^2))^2)*((vy[i-1] - u_c*(1 - (x[i-1]/r_0 - 2.5)^2)*exp(-(x[i-1]/r_0 - 2.5)^2))/(sqrt(vx[i-1]^2 + (vy[i-1] - u_c*(1 - (x[i-1]/r_0 - 2.5)^2)*exp(-(x[i-1]/r_0 - 2.5)^2))^2))))/m - g)) for i in 2:nh+1
    )

    ExaModels.constraint(
        c,
        x[1] - x_0
    )

    ExaModels.constraint(
        c,
        y[1] - y_0
    )

    ExaModels.constraint(
        c,
        y[nh+1] - y_f
    )

    ExaModels.constraint(
        c,
        vx[1] - vx_0
    )

    ExaModels.constraint(
        c,
        vx[nh+1] - vx_f
    )

    ExaModels.constraint(
        c,
        vy[1] - vy_0
    )

    ExaModels.constraint(
        c,
        vy[nh+1] - vy_f
    )

    ExaModels.ExaModel(c,; kwargs...)
end