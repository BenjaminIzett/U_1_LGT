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



    Vs = mean(ln_ratio[:, 4:end, :], dims=(1, 2,))
    V_stds = std(mean(ln_ratio[:, 4:end, :], dims=(2,)), dims=(1,))
    # display(-Vs[1, 1, :])
    # display(V_stds[1, 1, :])
    Vplot = scatter(R_values, -Vs[1, 1, :], yerrors=V_stds[1, 1, :], xlims=[0, R_values[end] + 1], ylims=[0.1, 0.8], markershape=:x, label="Mine")
    xlabel!(Vplot, "R")
    ylabel!(Vplot, "V(R)")

    function model(R, p)
        p[1] .+ p[2] .* R .+ p[3] * log.(R)
    end

    p0 = [0, 0.1, 0.1]

    Vs = mean(ln_ratio[:, 4:end, :], dims=(2,))
    fits = mapslices(slice -> LsqFit.curve_fit(model, collect(R_values), -slice, p0), Vs, dims=(3,))
    shaped_fits = hcat([fit.param for fit in fits[:, 1, 1]]...)
    params = mean(shaped_fits, dims=(2,))
    param_std = std(shaped_fits, dims=(2,))
    # params = fit.param
    display(params)
    display(param_std)
    model_Rs = 1:0.25:10
    # plot!(Vplot, model_Rs, model(model_Rs, params),label="Fit")
    hamer_data = [0.1718100890207715, 0.2626112759643917, 0.3373887240356083, 0.40445103857566767, 0.4667655786350148, 0.528486646884273, 0.5919881305637983, 0.6537091988130562, 0.7053412462908013, 0.7658753709198816]
    parisi_data = [-0.270280803927181, -0.341940855677859, -0.403856582168469, -0.461186508246316, -0.517907758987399, -0.577859171281773, -0.634255662731754, -0.712949807856125]
    parisi_data_up = [-0.269726982119319,-0.340904733659609,-0.402046743646662,-0.458269368513928,-0.512631882878799,-0.567254338464991,-0.612883732076678,-0.675755437848076]
    parisi_data_down = [-0.270837989157468,-0.34298944304376,-0.405703874534129,-0.464200069035835,-0.523471039247711,-0.589446398339986,-0.658832927365397,-0.762140052046897]

    parisi_errors = [0.00366015714714038, 0.00304208243034963, 0.00255081308372604, 0.00273536422256527, 0.00772093397443272, 0.012863480564136, 0.00245097916517631, 0]
    parisi_errors_2 = [0.000557185230286583,0.00104858736590124,0.00184729236565989,0.00301356078951859,0.00556328026031117,0.0115872270582126,0.0245772646336436,0.0491902441907718]
    scatter!(Vplot, 1:10, hamer_data, markershape=:x, label="Hamer")
    scatter!(Vplot, 1.8:8.8, -parisi_data, markershape=:x, label="Parisi")#,yerror=parisi_errors_2)
    # scatter!(Vplot, 2:9, -parisi_data_up, markershape=:x, label="ParisiU")#,yerror=parisi_errors)
    # scatter!(Vplot, 2:9, -parisi_data_down, markershape=:x, label="ParisiD")#,yerror=parisi_errors)

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

    loop_plot = scatter(R_values[1:end-1], -Vs[1, :, R], yerror=V_stds[1, :, R], markershape=:x)
    hline!(loop_plot, [mean(-Vs[1, 4:end, R])])
    display(mean(-Vs[1, 1:end, R]))

    loop_plot
end

# DataHandler.analyse_data("measurements/hamer_loop_check_2.txt","analysis/hamer_loop_check_fa09_2.txt",3:102,(mean for _ in 1:100))
# data=DataHandler.load_data("analysis/hamer_loop_check_2.txt")
# static_quark_potential(data,1:10,1:10)

# DataHandler.analyse_data("measurements/parisi_potential_check_fa099_β2_2.txt", "analysis/parisi_potential_check_fa099_β2_2.txt", 3:83, (mean for _ in 1:81))
data = DataHandler.load_data("analysis/parisi_potential_check_fa099_β2_1.txt")
static_quark_potential(data, 1:10, 1:10)

loop_ratio(data,1:10,1:10,8)
