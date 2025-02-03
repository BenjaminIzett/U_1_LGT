include("data_handler.jl")
include("autocorrelation.jl")

using Statistics

filename = "parisi_loop_check_fa099_β22_2.txt"

# filename = "hamer_loop_autocorrelation_fa09_β3_1.txt"
# filename = "hamer_loop_autocorrelation_β3_1.txt"

data = DataHandler.load_data_reshaped("measurements/$filename")
display(size(data))
N_repeats, run_length, N_different_measurements = size(data)

autocorrelations = Autocorrelation.calc_autocorrelation_indices(0:100, 3:90, data)

int_autocorr_times = Autocorrelation.int_autocorrelation_time.(autocorrelations, run_length)
# display(first.(int_autocorr_times))
run_to_analyse = 1
mean_values = mean(data[run_to_analyse, :, 3:end], dims=1)[1, :]

stderr_values = std(data[run_to_analyse, :, 3:end], dims=1)[1, :] / sqrt(run_length)

τ_ints = first.(int_autocorr_times)
errors = stderr_values .* sqrt.(2 * τ_ints)

δτ_ints = last.(int_autocorr_times)
DataHandler.save_data("autocorrelation_analysis/$filename", [mean_values, errors, stderr_values, τ_ints, δτ_ints])
