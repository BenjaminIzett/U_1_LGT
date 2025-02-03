include("../measurement.jl")
include("../measurement_functions.jl")
include("../u_1_lgt.jl")
include("../hmc.jl")
include("../data_handler.jl")
include("../autocorrelation.jl")

using CurveFit
using Statistics
using Plots


β = 3
trajectory_length = 1

# update = HMC.hmc_run
# update_args = (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.1, (β,))


update = HMC.fa_hmc_run_3d
κ = 0.99
inverse_FK = HMC.inv_FK_3d_kappa(16, 16, 16, κ)
update_args = (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.1, inverse_FK, (β,))
filename = "hamer_loop_autocorrelation_fa099_β3_τ2_1.txt"


lf_index = 3
Δτ_index = 4
measurement_functions = (MeasurementFunctions.hamer_loop_range_3d,)
measurement_args = ((β, 0.7, 10, [(1, 1), (2, 2), (3, 3), (7, 7)]),)
measurement_info = (β, 16)


if isfile("measurements/$filename")
    @info "Found measurement file, skipping measurement step." filename
else

    Measurement.repeat_optimise_measure(filename, U_1_LGT.zero_lattice_3d(16, 16, 16), 250, 1, 5000, update, update_args, measurement_functions, measurement_args, measurement_info, lf_index, Δτ_index, trajectory_length, 250, 1000, 5)
end

# function autocorrelation_a(a, data)

#     N = length(data)

#     mean_data = mean(data)
#     Γ_0 = var(data)

#     Γ_a = (1 / (N - a)) * sum((data[1:N-a] .- mean_data) .* (data[1+a:N] .- mean_data))

#     Γ_a / Γ_0

# end
# function autocorrelation_range(a_range, data)
#     map(a -> autocorrelation_a(a, data), a_range)
# end
# function calc_autocorrelation_index(a_range, index, data)
#     mapslices(d -> autocorrelation_range(a_range, d), data[:, :, index], dims=2)
# end
# function calc_autocorrelation_indices(a_range, indices, data)
#     map(index -> calc_autocorrelation_index(a_range, index, data), indices)
# end
# function calc_autocorrelation(a_range, indices, data)
#     autocorrelation_range(a_range,data[:, indices, :])
# end

function _plot_autocrrelation(a_range, index, data, label, trajectory_length; plt=nothing, color=nothing)
    if plt === nothing
        plt = plot()
    end
    if color === nothing
        color = rand(collect(palette(:auto, 10)))
    end

    autocorrelation = Autocorrelation.calc_autocorrelation_index(a_range, index, data)
    # plt = plot()
    # for line in autocorrelation_by_index
    N_repeats = size(data)[1]

    # values = mean(autocorrelation, dims=1)[1, :]
    # errors = std(autocorrelation, dims=1)[1, :] ./ sqrt(N_repeats)
    # positive_cutoff = all(values .- errors .> 0) ? length(values) : findfirst(x -> x <= 0, values .- errors) - 1

    plt = plot!(plt, a_range / trajectory_length, mean(autocorrelation, dims=1)[1, :], ribbon=std(autocorrelation, dims=1)[1, :] ./ sqrt(N_repeats), xlabel="Lag time", ylabel="Autocorrelation time", label=label, color=color)
    # end
    # xs = 0:fit_cutoff-1

    # fit = linear_fit(xs, log.(mean(autocorrelation, dims=1)[1, 1:fit_cutoff]))
    # plt = plot!(plt, plot_xs, exp.(fit[1] .+ plot_xs .* fit[2]), label="Exponential fit, exponent: $(round(fit[2],digits=3))")

    # fits = mapslices(l -> linear_fit(fit_window, log.(l[fit_window.+1])), autocorrelation, dims=2)

    # p1 = map(fit -> fit[1], fits)
    # p2 = map(fit -> fit[2], fits)
    # mean_fit = (mean(p1), mean(p2))
    # display(mean_fit)
    # plt = plot!(plt, xs, exp.(mean_fit[1] .+ xs .* mean_fit[2]), label="Exponential fit, exponent: $(round(mean_fit[2],digits=3)),std: $(round(std(p2),digits=3))")
    # fit = linear_fit(fit_window, log.(values[fit_window.+1]))

    # fit_values = fit[1] .+ a_range[1:positive_cutoff] .* fit[2]

    # plt = plot!(plt, a_range[1:positive_cutoff], exp.(fit_values), label="Linear fit, slope: $(round(fit[2],digits=3))")
    # plt = scatter!(plt, [fit_window[1], fit_window[end]], exp.(fit[1] .+ [fit_window[1], fit_window[end]] .* fit[2]), marker=:x, color=:orange, label="Fit window")
    plt
end

function plot_autocorrelation(a_range, index, data, label, trajectory_length; color=nothing)
    _plot_autocrrelation(a_range, index, data, label, trajectory_length; color=color)
end
function plot_autocorrelation!(plt, a_range, index, data, label, trajectory_length; color=nothing)
    _plot_autocrrelation(a_range, index, data, label, trajectory_length; plt=plt, color=color)
end




function _plot_log_autocorrelation(a_range, index, data, fit_window, label, trajectory_length; plt=nothing, color=nothing)
    if plt === nothing
        plt = plot()
    end
    if color === nothing
        color = rand(collect(palette(:auto, 10)))
    end

    autocorrelation = Autocorrelation.calc_autocorrelation_index(a_range, index, data)
    # plt = plot()
    # for line in autocorrelation_by_index
    N_repeats = size(data)[1]
    values = mean(autocorrelation, dims=1)[1, :]
    errors = std(autocorrelation, dims=1)[1, :] ./ sqrt(N_repeats)
    positive_cutoff = all(values .- errors .> 0) ? length(values) : findfirst(x -> x <= 0, values .- errors) - 1

    plt = plot!(plt, a_range[1:positive_cutoff] * trajectory_length, values[1:positive_cutoff], ribbon=errors[1:positive_cutoff], yscale=:log10, xlabel="Lag time", ylabel="Autocorrelation time", label=label, ylims=[0.01, 1], color=color)

    fit = linear_fit(fit_window, log.(values[fit_window.+1]))

    fit_values = fit[1] .+ a_range[1:positive_cutoff] .* fit[2]

    plt = plot!(plt, a_range[1:positive_cutoff], exp.(fit_values), yscale=:log10, label="Linear fit, slope: $(round(fit[2],digits=3))", color=color)
    plt = scatter!(plt, [fit_window[1], fit_window[end]], exp.(fit[1] .+ [fit_window[1], fit_window[end]] .* fit[2]), marker=:x, color=color, label="Fit window")
    plt

end
function plot_log_autocorrelation(a_range, index, data, fit_window, label, trajectory_length; color=nothing)
    _plot_log_autocorrelation(a_range, index, data, fit_window, label, trajectory_length; color=color)
end
function plot_log_autocorrelation!(plt, a_range, index, data, fit_window, label, trajectory_length; color=nothing)
    _plot_log_autocorrelation(a_range, index, data, fit_window, label, trajectory_length; plt=plt, color=color)
end

colors = collect(palette(:default))
data = DataHandler.load_data_reshaped("measurements/$filename")
# display(data)
plt = plot_autocorrelation(0:160, [4], data, filename, 2; color=colors[1])


# compare_filename="hamer_loop_autocorrelation_fa09_β22_1.txt"
# data = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt = plot_autocorrelation!(plt, 0:150, [3], data, 1,compare_filename)
# compare_filename="hamer_loop_autocorrelation_fa095_β22_1.txt"
# data = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt = plot_autocorrelation!(plt, 0:150, [3], data, 1,compare_filename)
# compare_filename="hamer_loop_autocorrelation_fa099_β22_1.txt"
# data = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt = plot_autocorrelation!(plt, 0:150, [3], data, 1,compare_filename)
# compare_filename="hamer_loop_autocorrelation_fa0999_β22_1.txt"
# data = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt = plot_autocorrelation!(plt, 0:150, [3], data, 1,compare_filename)
# display(data)


compare_filename = "hamer_loop_autocorrelation_fa099_β3_1.txt"
data_compare = DataHandler.load_data_reshaped("measurements/$compare_filename")
plt = plot_autocorrelation!(plt, 0:80, [4], data_compare, compare_filename, 1; color=colors[2])
plt |> display

plt2 = plot_log_autocorrelation(0:80, [3], data, 2:10, filename,1; color=colors[1])
plt2 = plot_log_autocorrelation!(plt2, 0:80, [3], data_compare, 2:10, compare_filename,1; color=colors[2])
plt2 |> display

# plot_log_autocorrelation(0:70, [3], data, 50) |> display
# mean_data = mean(data, dims=1)
# DataHandler.save_data("measurements/mean/$filename")


# DataHandler.analyse_data("measurements/$filename", "analysis/$filename", (1, 3, 3), (mean, mean, std))


#bigger kappa, use shorter trajectory length
