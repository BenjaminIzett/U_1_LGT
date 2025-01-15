include("../measurement.jl")
include("../measurement_functions.jl")
include("../u_1_lgt.jl")
include("../hmc.jl")
include("../data_handler.jl")
include("../autocorrelation.jl")

using CurveFit
using Statistics
using Plots

filename = "loop_autocorrelation_2.txt"

β = 2

update = HMC.hmc_run
update_args = (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.1, (β,))
lf_index = 3
Δτ_index = 4

# update = HMC.fa_hmc_run_3d
# κ = 0.5
# inverse_FK = HMC.inv_FK_3d_kappa(16, 16, 16, κ)
# update_args = (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.1, inverse_FK, (β,))
# lf_index = 3
# Δτ_index = 4
measurement_functions = (MeasurementFunctions.parisi_loop_range_3d,)
measurement_args = ((β, [(3, 3), (3, 7), (7, 7)]),)
measurement_info = (β, 16)


if isfile("measurements/$filename")
    @info "Found measurement file, skipping measurement step." filename
else

    Measurement.repeat_optimise_measure(filename, U_1_LGT.zero_lattice_3d(16, 16, 16), 250, 1, 5000, update, update_args, measurement_functions, measurement_args, measurement_info, lf_index, Δτ_index, 250, 1000, 10)
end

function autocorrelation_a(a, data)

    N = length(data)

    mean_data = mean(data)
    Γ_0 = var(data)

    Γ_a = (1 / (N - a)) * sum((data[1:N-a] .- mean_data) .* (data[1+a:N] .- mean_data))

    Γ_a / Γ_0

end
function autocorrelation_range(a_range, data)
    map(a -> autocorrelation_a(a, data), a_range)
end
function calc_autocorrelation_index(a_range, index, data)
    mapslices(d -> autocorrelation_range(a_range, d), data[:, :, index], dims=2)
end
function calc_autocorrelation_indices(a_range, indices, data)
    map(index -> calc_autocorrelation_index(a_range, index, data), indices)
end
# function calc_autocorrelation(a_range, indices, data)
#     autocorrelation_range(a_range,data[:, indices, :])
# end

function plot_autocorrelation(a_range, index, data, fit_cutoff)

    autocorrelation = calc_autocorrelation_index(a_range, index, data)
    # plt = plot()
    # for line in autocorrelation_by_index
    plt = plot(a_range, mean(autocorrelation, dims=1)[1, :], ribbon=std(autocorrelation, dims=1)[1, :], xlabel="Lag time", ylabel="Autocorrelation time", label="data")
    # end
    xs = 0:fit_cutoff-1

    # fit = linear_fit(xs, log.(mean(autocorrelation, dims=1)[1, 1:fit_cutoff]))
    # plt = plot!(plt, plot_xs, exp.(fit[1] .+ plot_xs .* fit[2]), label="Exponential fit, exponent: $(round(fit[2],digits=3))")

    fits = mapslices(l -> linear_fit(xs, log.(l[1:fit_cutoff])), autocorrelation, dims=2)
    p1 = map(fit -> fit[1], fits)
    p2 = map(fit -> fit[2], fits)
    mean_fit = (mean(p1), mean(p2))
    # display(mean_fit)
    plt = plot!(plt, xs, exp.(mean_fit[1] .+ xs .* mean_fit[2]), label="Exponential fit, exponent: $(round(mean_fit[2],digits=3)),std: $(round(std(p2),digits=3))")
    plt
end
function plot_autocorrelation!(plt, a_range, index, data, fit_cutoff)

    autocorrelation = calc_autocorrelation_index(a_range, index, data)
    # plt = plot()
    # for line in autocorrelation_by_index
    plt = plot!(plt, a_range, mean(autocorrelation, dims=1)[1, :], ribbon=std(autocorrelation, dims=1)[1, :], xlabel="Lag time", ylabel="Autocorrelation time", label="data")
    # end
    xs = 0:fit_cutoff-1

    # fit = linear_fit(xs, log.(mean(autocorrelation, dims=1)[1, 1:fit_cutoff]))
    # plt = plot!(plt, plot_xs, exp.(fit[1] .+ plot_xs .* fit[2]), label="Exponential fit, exponent: $(round(fit[2],digits=3))")

    fits = mapslices(l -> linear_fit(xs, log.(l[1:fit_cutoff])), autocorrelation, dims=2)
    p1 = map(fit -> fit[1], fits)
    p2 = map(fit -> fit[2], fits)
    mean_fit = (mean(p1), mean(p2))
    # display(mean_fit)
    plt = plot!(plt, xs, exp.(mean_fit[1] .+ xs .* mean_fit[2]), label="Exponential fit, exponent: $(round(mean_fit[2],digits=3)),std: $(round(std(p2),digits=3))")
    plt
end

function plot_log_autocorrelation(a_range, index, data, fit_cutoff)

    autocorrelation = calc_autocorrelation_index(a_range, index, data)
    # plt = plot()
    # for line in autocorrelation_by_index
    plt = plot(a_range, log.(mean(autocorrelation, dims=1)[1, :]), ribbon=std(autocorrelation, dims=1)[1, :], xlabel="Lag time", ylabel="Autocorrelation time", label="data")
    # end
    xs = 1:fit_cutoff
    plot_xs = 0:fit_cutoff
    fit = linear_fit(xs, log.(mean(autocorrelation, dims=1)[1, 1:fit_cutoff]))
    plt = plot!(plt, plot_xs, fit[1] .+ plot_xs .* fit[2], label="Exponential fit, exponent: $(round(fit[2],digits=3))")
    plt
end

data = DataHandler.load_data_reshaped("measurements/$filename")
# display(data)
plt = plot_autocorrelation(0:70, [5], data, 30)

data = DataHandler.load_data_reshaped("measurements/loop_autocorrelation_fa05_2.txt")
# display(data)
plt = plot_autocorrelation!(plt, 0:70, [5], data, 30)
plt |> display
# plot_log_autocorrelation(0:70, [3], data, 50) |> display
# mean_data = mean(data, dims=1)
# DataHandler.save_data("measurements/mean/$filename")


# DataHandler.analyse_data("measurements/$filename", "analysis/$filename", (1, 3, 3), (mean, mean, std))
