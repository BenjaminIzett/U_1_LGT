using Plots
using CurveFit
include("data_handler.jl")


function plot_av_plaquette(filename; plot_kwargs...)
    data = DataHandler.load_data(filename)
    plt = plot(data[:, 1], data[:, 2], yerror=data[:, 3]; plot_kwargs...)
    xaxis!(plt, "β")
    yaxis!(plt, "<P>")
    return plt
end


function check_hamer(filename, show_small)
    hamer_data_P = [0.475, 0.629, 0.656, 0.704, 0.748, 0.790, 0.806, 0.834, 0.854, 0.869, 0.881]
    hamer_data_β = [1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0]
    plt = plot_av_plaquette(filename, label="Mine")
    plot!(plt, hamer_data_β, hamer_data_P, label="Hamer")
    if show_small
        plot!(plt, collect(0.2:0.1:1), collect((0.2:0.1:1)) / 2, label="Small behaviour")
    end
    plt
end

function plot_autocorrelation(filename, fit_cutoff; error_file="")
    data = DataHandler.load_data(filename)

    if error_file != ""
        errors = DataHandler.load_data(error_file)
        plt = plot(data[2:end], ribbon=errors[2:end], xlabel="Lag time", ylabel="Autocorrelation time", label="data")

    else
        plt = plot(data[2:end], xlabel="Lag time", ylabel="Autocorrelation time", label="data")
    end
    xs = 1:fit_cutoff
    plot_xs = 0:40
    fit = linear_fit(xs, log.(data[2:fit_cutoff+1]))
    plt = plot!(plt, plot_xs, exp.(fit[1] .+ plot_xs .* fit[2]), label="Exponential fit, exponent: $(round(fit[2],digits=3))")
end
function plot_autocorrelation_log(filename, linear_cutoff; error_file="")
    data = DataHandler.load_data(filename)[1, :]
    xs = collect(1:linear_cutoff)
    # fit = log_fit(xs, data[2:linear_cutoff])
    # display(log.(data[2:linear_cutoff]))
    fit = linear_fit(xs, log.(data[2:linear_cutoff+1]))
    display(fit)
    positive_cutoff = findfirst(n -> n <= 0, data)
    plt = plot(data[2:positive_cutoff-1], yaxis=:log, xlabel="Lag time", ylabel="Autocorrelation time", label="data")
    # display(xs)
    # display(fit[1] .+ xs .* fit[2])
    plt = plot!(plt, xs, exp.(fit[1] .+ xs .* fit[2]), yaxis=:log, label="Linear fit, slope: $(round(fit[2],digits=3))")
    # display(0 .+ xs * m)
    plt

end

function plot_av_plaquette_v_trajectory(filename)
    data = DataHandler.load_data(filename)

    plt = plot(data[:, 3])
    xaxis!(plt, "Trajectory")
    yaxis!(plt, "<P>")
    plt
end
function plt_hot_cold_plaquette_v_trajectory(filename_hot, filename_cold)
    data_hot = DataHandler.load_data(filename_hot)
    data_cold = DataHandler.load_data(filename_cold)
    plt = plot(data_hot[:, 3], label="Random initial field")
    plt = plot!(plt, data_cold[:, 3], label="Zero initial field")
    xaxis!(plt, "Trajectory")
    yaxis!(plt, "<P>")

    plt
end


check_hamer("analysis/test_a.txt", false) |> display
# check_hamer("analysis/test_full.txt", true) |> display
# plot_autocorrelation("analysis/autocorrelation_1_05_1000.txt") |> display
# plot_autocorrelation("analysis/autocorrelation_1_05_100.txt") |> display
# plot_autocorrelation("analysis/autocorrelation_1_1_100.txt") |> display

# plot_autocorrelation("analysis/autocorrelation_05_05_1000.txt") |> display
# plot_autocorrelation("analysis/autocorrelation_3_05_1000.txt") |> display

# plot_autocorrelation("analysis/autocorrelation_1_05_1000_1.txt") |> display
# plot_autocorrelation("analysis/autocorrelation_1_05_1000_2.txt") |> display

# plot_autocorrelation("analysis/autocorrelation_1_05_10000.txt", 10) |> display
# plot_autocorrelation_log("analysis/autocorrelation_1_05_10000.txt", 10) |> display

# plot_autocorrelation("analysis/mean/autocorrelation_mean.txt", 15; error_file="analysis/std/autocorrelation_mean.txt") |> display
# plot_autocorrelation_log("analysis/mean/autocorrelation_mean.txt", 15; error_file="analysis/std/autocorrelation_mean.txt") |> display


# plt_hot_cold_plaquette_v_trajectory("measurements/plaquette_v_trajectory_hot.txt", "measurements/plaquette_v_trajectory_cold.txt")
