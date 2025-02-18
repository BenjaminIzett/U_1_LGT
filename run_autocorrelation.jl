include("data_handler.jl")
include("autocorrelation.jl")

using Statistics


# filename = "hamer_loop_autocorrelation_fa09_β3_1.txt"
# filename = "hamer_loop_autocorrelation_β3_1.txt"
function run_autocorrelation(filename)
    data = DataHandler.load_data_reshaped("measurements/$filename")
    display(size(data))
    N_repeats, run_length, N_different_measurements = size(data)

    data_range = 5:24
    autocorrelations = Autocorrelation.calc_autocorrelation_indices(0:100, data_range, data)

    int_autocorr_times = Autocorrelation.int_autocorrelation_time.(autocorrelations, run_length)
    # display(first.(int_autocorr_times))
    for run_to_analyse in 1:N_repeats
        mean_values = mean(data[run_to_analyse, :, data_range], dims=1)[1, :]

        stderr_values = std(data[run_to_analyse, :, data_range], dims=1)[1, :] / sqrt(run_length)

        τ_ints = first.(int_autocorr_times)
        errors = stderr_values .* sqrt.(2 * τ_ints)

        δτ_ints = last.(int_autocorr_times)
        DataHandler.save_data_custom("autocorrelation_analysis/$filename", [mean_values, errors, stderr_values, τ_ints, δτ_ints], "a", "5")
    end
end

# filename = "teper_fa099_β22_10_50.txt"

filenames = ["teper_fa099_β22_8_80.txt", "teper_fa099_β22_10_50.txt", "teper_fa099_β22_14_64.txt", "teper_fa099_β22_18_50.txt", "teper_fa099_β22_22_44.txt", "teper_fa099_β22_26_42.txt", "teper_fa099_β22_34_34.txt", "teper_fa099_β22_38_38.txt"]

for filename in filenames
    run_autocorrelation(filename)
end
