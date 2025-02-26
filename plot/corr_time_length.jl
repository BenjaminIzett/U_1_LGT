include("../measurement.jl")
include("../measurement_functions.jl")
include("../u_1_lgt.jl")
include("../hmc.jl")
include("../data_handler.jl")
include("../autocorrelation.jl")


βs = [2.0, 2.2, 2.4, 2.6, 2.8]
trajectory_lengths = [0.6]#[0.6, 0.7, 0.8, 0.9, 1.0, 1.2]
κs = [0.99]#[0.8, 0.9, 0.95, 0.99]
loop_index =2

plot_βs = [1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8]
plot_ξs = sqrt.([
    6.01319873068591,
    9.18273645546373,
    10.9716962055815,
    14.9338625002154,
    20.1994697639187,
    27.2429620334371,
    36.0244924904622,
    47.7400506881898,
    64.8793592644238,
    83.3554162669879,
    108.82411667791,
    136.411050364542,
    171.411790704009
])
τ_ints = []
δτ_ints = []
for trajectory_length in trajectory_lengths

    ξs = map(β -> plot_ξs[findfirst(v -> v == β, plot_βs)], βs)
    for β in βs
        filename = "loop_autocorrelation_β$(β)_τ$(trajectory_length).txt"
        data = DataHandler.load_data_reshaped("autocorrelation_analysis/$filename")
        push!(τ_ints, mean(data[:, 4, loop_index]))
        push!(δτ_ints, mean(data[:, 5, loop_index]))
    end
    plot(ξs, τ_ints, yerrors=δτ_ints, xaxis=:log, yaxis=:log) |> display
end

for κ in κs
    for trajectory_length in trajectory_lengths
        fa_τ_ints = []
        fa_δτ_ints = []
        ξs = map(β -> plot_ξs[findfirst(v -> v == β, plot_βs)], βs)
        for β in βs
            filename = "loop_autocorrelation_fa$(κ)_β$(β)_τ$(trajectory_length).txt"
            data = DataHandler.load_data_reshaped("autocorrelation_analysis/$filename")
            push!(fa_τ_ints, mean(data[:, 4, loop_index]))
            push!(fa_δτ_ints, mean(data[:, 5, loop_index]))
        end
        plt = plot(ξs, fa_τ_ints, yerrors=fa_δτ_ints, xaxis=:log, yaxis=:log)
        plt = plot!(plt, ξs, τ_ints, yerrors=δτ_ints, xaxis=:log, yaxis=:log)
        plt |> display
    end
end

plot_βs = [1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8]
plot_ξs = sqrt.([
    6.01319873068591,
    9.18273645546373,
    10.9716962055815,
    14.9338625002154,
    20.1994697639187,
    27.2429620334371,
    36.0244924904622,
    47.7400506881898,
    64.8793592644238,
    83.3554162669879,
    108.82411667791,
    136.411050364542,
    171.411790704009
])
τ_ints = []
δτ_ints = []
for trajectory_length in trajectory_lengths

    ξs = map(β -> plot_ξs[findfirst(v -> v == β, plot_βs)], βs)
    for β in βs
        filename = "loop_autocorrelation_β$(β)_τ$(trajectory_length).txt"
        data = DataHandler.load_data_reshaped("autocorrelation_analysis/$filename")
        push!(τ_ints, mean(data[:, 4, loop_index]))
        push!(δτ_ints, mean(data[:, 5, loop_index]))
    end
    plot(ξs, τ_ints, yerrors=δτ_ints, xaxis=:log, yaxis=:log) |> display
end

for κ in κs
    for trajectory_length in trajectory_lengths
        fa_τ_ints = []
        fa_δτ_ints = []
        ξs = map(β -> plot_ξs[findfirst(v -> v == β, plot_βs)], βs)
        for β in βs
            filename = "loop_autocorrelation_fa$(κ)_β$(β)_τ$(trajectory_length).txt"
            data = DataHandler.load_data_reshaped("autocorrelation_analysis/$filename")
            push!(fa_τ_ints, mean(data[:, 4, loop_index]))
            push!(fa_δτ_ints, mean(data[:, 5, loop_index]))
        end
        plt = plot(ξs, fa_τ_ints, yerrors=fa_δτ_ints, xaxis=:log, yaxis=:log)
        plt = plot!(plt, ξs, τ_ints, yerrors=δτ_ints, xaxis=:log, yaxis=:log)
        plt |> display
    end
end
