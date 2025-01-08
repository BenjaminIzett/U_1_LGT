include("../measurement.jl")
include("../measurement_functions.jl")
include("../u_1_lgt.jl")
include("../hmc.jl")
include("../data_handler.jl")

using Statistics
using Plots

filename = "hamer_loops_1.txt"


update_args_list = [(U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.05, (β,)) for β in dataβ]
measurement_args_list = [((β,),) for β in dataβ]
measurement_info_list = [(β, 16) for β in dataβ]

if isfile("measurements/$filename")
    @info "Found measurement file, skipping measurement step." filename
else
    Measurement.optimised_measure_range(filename, U_1_LGT.zero_lattice_3d(16, 16, 16), 250, 10, 50, HMC.hmc_run, update_args_list, (MeasurementFunctions.mean_plaquette_3d,), measurement_args_list, measurement_info_list, 3, 4, 250, 1000)
end

DataHandler.analyse_data("measurements/$filename", "analysis/$filename", (1, 3, 3), (mean, mean, std))

function plot_av_plaquette(filename; plot_kwargs...)
    data = DataHandler.load_data(filename)
    plt = scatter(data[:, 1], data[:, 2], yerror=data[:, 3]; plot_kwargs...)
    xaxis!(plt, "β")
    yaxis!(plt, "<P>")
    return plt
end


function check_hamer(filename, show_small)
    hamer_data_P = [0.475, 0.629, 0.656, 0.704, 0.748, 0.790, 0.806, 0.834, 0.854, 0.869, 0.881]
    hamer_data_β = [1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0]
    plt = plot_av_plaquette(filename, label="Mine", markershape=:x)
    scatter!(plt, hamer_data_β, hamer_data_P, label="Hamer", markershape=:x)
    if show_small
        plot!(plt, collect(0.2:0.1:1), collect((0.2:0.1:1)) / 2, label="Small behaviour")
    end
    plt
end

check_hamer("analysis/$filename", true)
