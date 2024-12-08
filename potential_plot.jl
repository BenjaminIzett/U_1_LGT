include("data_handler.jl")

using Statistics
using LsqFit
using Plots


function static_quark_potential(data, R_values, τ_values)
    (Nm, NRτ) = size(data)
    if NRτ != length(R_values) * length(τ_values)
        display("Incorrect data size for given R and τ values")
    end

    shaped_data = reshape(data, (Nm, length(τ_values), length(R_values)))
    #shaped_data[:,:,R] gives data for R loops
    ratio = shaped_data[:, 2:end, :] ./ shaped_data[:, 1:end-1, :]

    #may need to remove negative ratios
    ln_ratio = log.(ratio)
    # mean_values = mean(ln_ratio, dims=(1,))
    # std_values = std(ln_ratio, dims=(1,))



    Vs = mean(ln_ratio[:, 3:end, :], dims=(1, 2,))
    V_stds = std(mean(ln_ratio[:, 3:end, :], dims=(2,)), dims=(1,))
    # display(-Vs[1, 1, :])
    # display(V_stds[1, 1, :])
    Vplot = scatter(R_values, -Vs[1, 1, :], yerrors=V_stds[1, 1, :], xlims=[0, R_values[end] + 1], ylims=[0.1, -Vs[1, 1, end] + 0.1], markershape=:x, label="Mine")
    xlabel!(Vplot, "R")
    ylabel!(Vplot, "V(R)")

    function model(R, p)
        p[1] .+ p[2] .* R .+ p[3] * log.(R)
    end

    p0 = [0, 0.1, 0.1]

    Vs = mean(ln_ratio[:, 3:end, :], dims=(2,))
    fits = mapslices(slice -> curve_fit(model, collect(R_values), -slice, p0), Vs, dims=(3,))
    shaped_fits = hcat([fit.param for fit in fits[:, 1, 1]]...)
    params = mean(shaped_fits, dims=(2,))
    param_std = std(shaped_fits, dims=(2,))
    # params = fit.param
    display(params)
    display(param_std)
    model_Rs = 1:0.25:10
    plot!(Vplot, model_Rs, model(model_Rs, params),label="Fit")
    hamer_data = [0.1718100890207715, 0.2626112759643917, 0.3373887240356083, 0.40445103857566767, 0.4667655786350148, 0.528486646884273, 0.5919881305637983, 0.6537091988130562, 0.7053412462908013, 0.7658753709198816]

    scatter!(Vplot, 1:10, hamer_data,markershape=:x,label="Hamer")

    Vplot
end

function loop_ratio(data, R_values, τ_values, R)
    (Nm, NRτ) = size(data)
    if NRτ != length(R_values) * length(τ_values)
        display("Incorrect data size for given R and τ values")
    end


    shaped_data = reshape(data, (Nm, length(τ_values), length(R_values)))
    #shaped_data[:,:,R] gives data for R loops
    ratio = shaped_data[:, 2:end, :] ./ shaped_data[:, 1:end-1, :]

    #may need to remove negative ratios
    ln_ratio = log.(ratio)

    Vs = mean(ln_ratio, dims=(1,))
    display(Vs[:, :, R])
    V_stds = std(ln_ratio, dims=(1,))

    loop_plot = scatter(R_values[1:end-1], -Vs[1, :, R], yerror=V_stds[1, :, R], markershape=:x, ylims=(-0.5, 0.5), yticks=-0.4:0.2:0.4)
    hline!(loop_plot, [mean(-Vs[1, 1:end, R])])
    display(mean(-Vs[1, 1:end, R]))

    loop_plot
end
