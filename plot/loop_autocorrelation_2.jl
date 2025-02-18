include("../measurement.jl")
include("../measurement_functions.jl")
include("../u_1_lgt.jl")
include("../hmc.jl")
include("../data_handler.jl")
include("../autocorrelation.jl")

using CurveFit
using Statistics
using Plots


function loop_autocorrelation_fa(β, trajectory_length, κ, filename, measurement_functions, measurement_args, measurement_info)

    update = HMC.fa_hmc_run_3d

    inverse_FK = HMC.inv_FK_3d_kappa(32, 32, 32, κ)
    update_args = (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.1, inverse_FK, (β,))

    lf_index = 3
    Δτ_index = 4


    if isfile("measurements/$filename")
        @info "Found measurement file, skipping measurement step." filename
    else
        Measurement.repeat_optimise_measure(filename, U_1_LGT.zero_lattice_3d(32, 32, 32), 250, 1, 10000, update, update_args, measurement_functions, measurement_args, measurement_info, lf_index, Δτ_index, trajectory_length, 500, 1000, 5)
    end

end

function loop_autocorrelation(β, trajectory_length, filename, measurement_functions, measurement_args, measurement_info)


    update = HMC.hmc_run
    update_args = (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.1, (β,))

    lf_index = 3
    Δτ_index = 4


    if isfile("measurements/$filename")
        @info "Found measurement file, skipping measurement step." filename
    else
        Measurement.repeat_optimise_measure(filename, U_1_LGT.zero_lattice_3d(32, 32, 32), 250, 1, 10000, update, update_args, measurement_functions, measurement_args, measurement_info, lf_index, Δτ_index, trajectory_length, 500, 1000, 5)
    end

end

βs = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
trajectory_lengths = [0.6, 0.7, 0.8, 0.9, 1.0, 1.2]
κs = [0.8, 0.9, 0.95, 0.99]

combinations = [(β, tl, κ) for β in βs for tl in trajectory_lengths for κ in κs]
for comb in combinations
    β, trajectory_length, κ = comb#combinations[parse(Int, ARGS[1])]

    filename = "loop_autocorrelation_fa$(κ)_β$(β)_τ$(trajectory_length).txt"
    measurement_functions = (MeasurementFunctions.mean_plaquette_3d, MeasurementFunctions.parisi_loop_range_3d,)
    measurement_args = ((β,), (β, [(2, 2), (4, 4), (8, 8), (16, 16)]),)
    measurement_info = (β, 32, κ, trajectory_length)

    loop_autocorrelation_fa(β, trajectory_length, κ, filename, measurement_functions, measurement_args, measurement_info)
end
