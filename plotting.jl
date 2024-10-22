using Plots
include("data_handler.jl")


function plot_av_plaquette(filename; plot_kwargs...)
    data = DataHandler.load_data(filename)
    plt = plot(data[:, 1], data[:, 2], yerror=data[:, 3]; plot_kwargs...)
    xaxis!(plt, "β")
    yaxis!(plt, "<P>")
    return plt
end


function check_hamer(filename)
    hamer_data_P = [0.475, 0.629, 0.656, 0.704, 0.748, 0.790, 0.806, 0.834, 0.854, 0.869, 0.881]
    hamer_data_β = [1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0]
    plot = plot_av_plaquette(filename, label="Mine")
    plot!(plot, hamer_data_β, hamer_data_P, label="Hamer")
end

check_hamer("analysis/test_a.txt")
