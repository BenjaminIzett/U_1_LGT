include("../measurement.jl")
include("../measurement_functions.jl")
include("../u_1_lgt.jl")
include("../hmc.jl")
include("../data_handler.jl")
include("../autocorrelation.jl")
include("../aic.jl")

using CurveFit
using Statistics
using Plots


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

    plt = plot!(plt, a_range * trajectory_length, mean(autocorrelation, dims=1)[1, :], ribbon=std(autocorrelation, dims=1)[1, :] ./ sqrt(N_repeats), xlabel="Lag time", ylabel="Autocorrelation time", label=label, color=color)
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
    positive_cutoff = (all(values .- errors .> 0) ? length(values) : findfirst(x -> x <= 0, values .- errors)) - 1

    plt = plot!(plt, a_range[1:positive_cutoff] * trajectory_length, values[1:positive_cutoff], ribbon=errors[1:positive_cutoff], yscale=:log10, xlabel="Lag time", ylabel="Autocorrelation time", label=label, ylims=[0.01, 1], color=color, size=(800, 600))

    lower_bound = findfirst(v -> v >= first(fit_window), a_range * trajectory_length)
    upper_bound = findfirst(v -> v >= last(fit_window), a_range * trajectory_length)
    # fit = linear_fit((a_range[lower_bound:upper_bound] * trajectory_length), log.(values[lower_bound:upper_bound]))

    # fit_values = fit[1] .+ a_range[1:positive_cutoff] * trajectory_length .* fit[2]

    # plt = plot!(plt, a_range[1:positive_cutoff] * trajectory_length, exp.(fit_values), yscale=:log10, label="Linear fit, slope: $(round(fit[2],digits=3))", color=color)
    # plt = scatter!(plt, [a_range[lower_bound] * trajectory_length, a_range[upper_bound] * trajectory_length], exp.(fit[1] .+ [a_range[lower_bound] * trajectory_length, a_range[upper_bound] * trajectory_length] .* fit[2]), marker=:x, color=color, label="Fit window")

    #AIC version

    fit = AIC.fit_window(Polynomial, (1,), a_range[1:positive_cutoff] * trajectory_length, log.(values[1:positive_cutoff]), errors[1:positive_cutoff], 2)

    plt = plot!(plt, a_range[1:positive_cutoff] * trajectory_length, exp.(fit.(a_range[1:positive_cutoff] * trajectory_length)), label="Linear fit, slope: $(round(fit.coeffs[2],digits=3))", color=color)
    # display(fit)
    (plt, fit.coeffs[2])

end
function plot_log_autocorrelation(a_range, index, data, fit_window, label, trajectory_length; color=nothing)
    _plot_log_autocorrelation(a_range, index, data, fit_window, label, trajectory_length; color=color)
end
function plot_log_autocorrelation!(plt, a_range, index, data, fit_window, label, trajectory_length; color=nothing)
    _plot_log_autocorrelation(a_range, index, data, fit_window, label, trajectory_length; plt=plt, color=color)
end

colors = collect(palette(:default))
# filename = "hamer_loop_autocorrelation_fa09_β22_1.txt"
# data = DataHandler.load_data_reshaped("measurements/$filename")
# display(data)
# plt = plot_autocorrelation(0:15, [4], data, filename, 0.85; color=colors[1])


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


# compare_filename = "hamer_loop_autocorrelation_fa09_β3_1.txt"
# data_compare = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt = plot_autocorrelation!(plt, 0:15, [4], data_compare, compare_filename, 1; color=colors[2])

# compare_filename = "hamer_loop_autocorrelation_fa099_β3_1.txt"
# data_compare = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt = plot_autocorrelation!(plt, 0:15, [4], data_compare, compare_filename, 1; color=colors[3])


# compare_filename = "hamer_loop_autocorrelation_fa09_β3_1.txt"
# data_compare = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt = plot_autocorrelation!(plt, 0:15, [4], data_compare, compare_filename, 1; color=colors[2])
# compare_filename = "hamer_loop_autocorrelation_fa09_β3_τ075_1.txt"
# data_compare = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt = plot_autocorrelation!(plt, 0:15, [4], data_compare, compare_filename, 0.75; color=colors[3])
# compare_filename = "hamer_loop_autocorrelation_fa09_β3_τ06_1.txt"
# data_compare = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt = plot_autocorrelation!(plt, 0:15, [4], data_compare, compare_filename, 0.6; color=colors[4])
# compare_filename = "hamer_loop_autocorrelation_fa09_β3_τ08_1.txt"
# data_compare = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt = plot_autocorrelation!(plt, 0:15, [4], data_compare, compare_filename, 0.8; color=colors[5])

# plt |> display


βs = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
trajectory_lengths = [0.6, 0.7, 0.8, 0.9, 1.0, 1.2]
κs = [0.8, 0.9, 0.95, 0.99]
loop = [5]

plt = plot()
layout = @layout [a b; c d; e f]
for β in βs
    plts = []
    coeffs = []
    coeffs_fa = []
    for tl in trajectory_lengths
        filename = "loop_autocorrelation_β$(β)_τ$(tl).txt"
        data = DataHandler.load_data_reshaped("measurements/$filename")

        xmax = 20
        upper = Int(round(xmax / tl))
        local plt = plot_autocorrelation(0:upper, loop, data, filename, tl, color=colors[1])
        # local (plt, coeff) = plot_log_autocorrelation(0:upper, loop, data, (1, 2), filename, tl, color=colors[1])
        # push!(coeffs, (β, tl, coeff))
        Plots.xlims!(plt, (0, xmax))
        Plots.ylims!(plt, (0, 1))
        for (i, κ) in enumerate(κs)
            # display((β, tl, κ))
            filename_fa = "loop_autocorrelation_fa$(κ)_β$(β)_τ$(tl).txt"
            data = DataHandler.load_data_reshaped("measurements/$filename_fa")

            plt = plot_autocorrelation!(plt, 0:upper, loop, data, filename_fa, tl, color=colors[i+1])
            # (plt, coeff_fa) = plot_log_autocorrelation!(plt, 0:upper, loop, data, (1, 2), filename_fa, tl, color=colors[i+1])
            # push!(coeffs_fa, (β, tl, κ, coeff_fa))
        end
        push!(plts, plt)
    end
    # display(coeffs[argmin(getindex.(coeffs, 3))])
    # display(coeffs_fa[argmin(getindex.(coeffs_fa, 4))])
    plot(plts..., layout=layout, size=(1200, 1800)) |> display

end

# compare_filename = "hamer_loop_autocorrelation_fa09_β3_τ085_1.txt"
# data_compare = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt2 = plot_log_autocorrelation(0:15, [4], data_compare, (1, 3.2), compare_filename, 0.85; color=colors[1])
# compare_filename = "hamer_loop_autocorrelation_fa09_β3_1.txt"
# data_compare = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt2 = plot_log_autocorrelation!(plt2, 0:15, [4], data_compare, (3, 7), compare_filename, 1; color=colors[2])
# compare_filename = "hamer_loop_autocorrelation_fa09_β3_τ075_1.txt"
# data_compare = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt2 = plot_log_autocorrelation!(plt2, 0:15, [4], data_compare, (1, 3.2), compare_filename, 0.75; color=colors[3])
# compare_filename = "hamer_loop_autocorrelation_fa09_β3_τ06_1.txt"
# data_compare = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt2 = plot_log_autocorrelation!(plt2, 0:15, [4], data_compare, (1.6, 4), compare_filename, 0.6; color=colors[4])
# compare_filename = "hamer_loop_autocorrelation_fa09_β3_τ08_1.txt"
# data_compare = DataHandler.load_data_reshaped("measurements/$compare_filename")
# plt2 = plot_log_autocorrelation!(plt2, 0:15, [4], data_compare, (1.6, 3.2), compare_filename, 0.8; color=colors[5])

# plt2 = plot_log_autocorrelation(0:80, [3], data, 2:10, filename,1; color=colors[1])
# plt2 = plot_log_autocorrelation!(plt2, 0:80, [3], data_compare, 2:10, compare_filename,1; color=colors[2])
# plt2 |> display

# plot_log_autocorrelation(0:70, [3], data, 50) |> display
# mean_data = mean(data, dims=1)
# DataHandler.save_data("measurements/mean/$filename")


# DataHandler.analyse_data("measurements/$filename", "analysis/$filename", (1, 3, 3), (mean, mean, std))
