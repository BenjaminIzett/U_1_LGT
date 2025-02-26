include("data_handler.jl")
include("autocorrelation.jl")

using Statistics


# filename = "hamer_loop_autocorrelation_fa09_β3_1.txt"
# filename = "hamer_loop_autocorrelation_β3_1.txt"
function run_autocorrelation(filename, a_range, data_range)
    data = DataHandler.load_data_reshaped("measurements/$filename")
    display(size(data))
    N_repeats, run_length, N_different_measurements = size(data)

    # data_range = 5:24
    autocorrelations = Autocorrelation.calc_autocorrelation_indices(a_range, data_range, data)
    #(data_range,[N_repeats,a_range])
    # display(autocorrelations)
    # display(size(autocorrelations[1]))

    # display(int_autocorr_times)
    for run_to_analyse in 1:N_repeats
        mean_values = mean(data[run_to_analyse, :, data_range], dims=1)[1, :]

        stderr_values = std(data[run_to_analyse, :, data_range], dims=1)[1, :] / sqrt(run_length)
        run_autocorrelations = map(t -> t[run_to_analyse, :], autocorrelations)
        int_autocorr_times = Autocorrelation.int_autocorrelation_time.(run_autocorrelations, run_length)
        τ_ints = first.(int_autocorr_times)
        errors = stderr_values .* sqrt.(2 * τ_ints)

        δτ_ints = last.(int_autocorr_times)
        DataHandler.save_data_custom("autocorrelation_analysis/$filename", [mean_values, errors, stderr_values, τ_ints, δτ_ints], "a", "5")
    end
end

# filename = "teper_fa099_β22_10_50.txt"

# filenames = ["teper_fa099_β22_8_80.txt", "teper_fa099_β22_10_50.txt", "teper_fa099_β22_14_64.txt", "teper_fa099_β22_18_50.txt", "teper_fa099_β22_22_44.txt", "teper_fa099_β22_26_42.txt", "teper_fa099_β22_34_34.txt", "teper_fa099_β22_38_38.txt"]

# for filename in filenames
#     run_autocorrelation(filename)
# end

#WILSON LINES
βs = [2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]
a_range = 0:1000
data_range = [2, 3]
for β in βs

    filename = "wilson_line_abs_fa0.99_β$(β)_τ0.6_2.txt"
    run_autocorrelation(filename, a_range, data_range)
end

#LOOP AUTOCORRELATIONS
# βs = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
# trajectory_lengths = [0.6, 0.7, 0.8, 0.9, 1.0, 1.2]
# κs = [0.8, 0.9, 0.95, 0.99]

# combinations = [(β, tl, κ) for β in βs for tl in trajectory_lengths for κ in κs]
# for comb in combinations
#     β, trajectory_length, κ = comb#combinations[parse(Int, ARGS[1])]

#     filename = "loop_autocorrelation_fa$(κ)_β$(β)_τ$(trajectory_length).txt"
#     run_autocorrelation(filename, 0:10000, [2, 3, 4, 5,6])
# end

# βs = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
# trajectory_lengths = [0.6, 0.7, 0.8, 0.9, 1.0, 1.2]

# combinations = [(β, tl) for β in βs for tl in trajectory_lengths]
# for comb in combinations
#     β, trajectory_length = comb#combinations[parse(Int, ARGS[1])]

#     filename = "loop_autocorrelation_β$(β)_τ$(trajectory_length).txt"
#     run_autocorrelation(filename, 0:10000, [2, 3, 4, 5, 6])
# end
