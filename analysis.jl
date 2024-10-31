include("data_handler.jl")

# DataHandler.autocorrelation(0:70, "test_full.txt")
# DataHandler.autocorrelation(0:70, "autocorrelation_1_05_1000.txt")
# DataHandler.autocorrelation(0:70, "autocorrelation_1_05_100.txt")
# DataHandler.autocorrelation(0:70, "autocorrelation_1_1_100.txt")
# DataHandler.autocorrelation(0:70, "autocorrelation_05_05_1000.txt")
# DataHandler.autocorrelation(0:70, "autocorrelation_1_05_1000_1.txt")
# DataHandler.autocorrelation(0:70, "autocorrelation_1_05_1000_2.txt")


# DataHandler.autocorrelation(0:70, "autocorrelation_1_05_10000.txt")


DataHandler.specific_heat_capacity("plaquette_heat_capacity.txt")


function int_autocorrelation_time(ρ_a, N)
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

    return τ_int, δτ_int
end

function int_autocorrelation_time(filename)
    data = DataHandler.load_data("analysis/$filename")
    N = data[1]
    ρ_a = data[2:end]
    int_autocorrelation_time(ρ_a, N)
end

int_autocorrelation_time("autocorrelation_1_05_10000.txt")

# DataHandler.analyse_data("measurements/long_test_full.txt", "analysis/long_test_full.txt", (3, 3,), (mean, std,))
