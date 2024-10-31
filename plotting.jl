using Plots
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

function plot_autocorrelation(filename)
    data = DataHandler.load_data(filename)
    plt = plot(data[2:end], xlabel="Lag time", ylabel="Autocorrelation time")

end
function plot_autocorrelation_log(filename)
    data = DataHandler.load_data(filename)
    plt = plot(abs.(data[2:end]), xaxis="Lag time", yaxis=:log, ylabel="Autocorrelation time")

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


# check_hamer("analysis/test_a.txt", false) |> display
# check_hamer("analysis/test_full.txt", true) |> display
# plot_autocorrelation("analysis/autocorrelation_1_05_1000.txt") |> display
# plot_autocorrelation("analysis/autocorrelation_1_05_100.txt") |> display
# plot_autocorrelation("analysis/autocorrelation_1_1_100.txt") |> display

# plot_autocorrelation("analysis/autocorrelation_05_05_1000.txt") |> display
# plot_autocorrelation("analysis/autocorrelation_3_05_1000.txt") |> display

# plot_autocorrelation("analysis/autocorrelation_1_05_1000_1.txt") |> display
# plot_autocorrelation("analysis/autocorrelation_1_05_1000_2.txt") |> display

plot_autocorrelation("analysis/autocorrelation_1_05_10000.txt") |> display
plot_autocorrelation_log("analysis/autocorrelation_1_05_10000.txt") |> display

plt_hot_cold_plaquette_v_trajectory("measurements/plaquette_v_trajectory_hot.txt", "measurements/plaquette_v_trajectory_cold.txt")
