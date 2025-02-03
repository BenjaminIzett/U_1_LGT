module Autocorrelation

using Statistics

include("data_handler.jl")

function autocorrelation_a(a, data)

    N = length(data)

    mean_data = mean(data)
    Γ_0 = var(data)

    Γ_a = (1 / (N - a)) * sum((data[1:N-a] .- mean_data) .* (data[1+a:N] .- mean_data))

    Γ_a / Γ_0

end
function autocorrelation(a_range, filename)
    DataHandler.analyse_data("measurements/$filename", "analysis/$filename", tuple(1, [5 for _ in a_range]...), tuple(length, [data -> autocorrelation_a(n, data) for n in a_range]...))

end
function autocorrelation(a_range, index, filename)
    DataHandler.analyse_data("measurements/$filename", "analysis/$filename", tuple(1, [index for _ in a_range]...), tuple(length, [data -> autocorrelation_a(n, data) for n in a_range]...))
end

function autocorrelation(a_range, indices, filename)
    DataHandler.analyse_data("measurements/$filename", "analysis/$filename", tuple(1, Iterators.flatten([[index for _ in a_range] for index in indices])...), tuple(length, Iterators.flatten([[data -> autocorrelation_a(n, data) for n in a_range] for _ in indices])...))

end


function autocorrelation_range(a_range, data)
    map(a -> autocorrelation_a(a, data), a_range)
end
function calc_autocorrelation_index(a_range, index, data)
    mapslices(d -> autocorrelation_range(a_range, d), data[:, :, index], dims=2)
end
function calc_autocorrelation_indices(a_range, indices, data)
    map(index -> calc_autocorrelation_index(a_range, index, data), indices)
end


function int_autocorrelation_time(ρ_a, N; print_M=false)
    # The search for M is not optimal but good enough
    M = 1
    τ_int = 1 / 2 + ρ_a[1]
    while M < 4 * τ_int + 1
        M += 1
        if M > length(ρ_a)
            print("ρ_a needs to be calculated for larger values of a")
            break
        else
            τ_int += ρ_a[M]
        end
    end
    δτ_int = sqrt((4 * M + 2) / N) * τ_int
    if print_M
        display(M)
    end
    return τ_int, δτ_int
end


function int_autocorrelation_time_file(filename)
    data = DataHandler.load_data("analysis/$filename")
    N = data[1]
    ρ_a = data[2:end]
    return int_autocorrelation_time(ρ_a, N)
end

end


# Autocorrelation.autocorrelation(0:70, "test_full.txt")
# Autocorrelation.autocorrelation(0:70, "autocorrelation_1_05_1000.txt")
# Autocorrelation.autocorrelation(0:70, "autocorrelation_1_05_100.txt")
# Autocorrelation.autocorrelation(0:70, "autocorrelation_1_1_100.txt")
# Autocorrelation.autocorrelation(0:70, "autocorrelation_05_05_1000.txt")
# Autocorrelation.autocorrelation(0:70, "autocorrelation_1_05_1000_1.txt")
# Autocorrelation.autocorrelation(0:70, "autocorrelation_1_05_1000_2.txt")


# Autocorrelation.autocorrelation(0:70, "autocorrelation_1_05_10000.txt")
# for n in 1:10
#     Autocorrelation.autocorrelation(0:70, "autocorrelation/autocorrelation_1_10000_$(n).txt")
# end

# Autocorrelation.specific_heat_capacity("plaquette_heat_capacity.txt")
