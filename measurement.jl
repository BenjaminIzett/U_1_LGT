include("action.jl")
include("hmc.jl")
include("data_handler.jl")

using Statistics
import Plots
using DelimitedFiles
using LinearAlgebra

#initial field
#update field
#steps_to_thermalise
#measurement interval
#number of measurements
#list of measurent functions
#update_field args
#measurement args


function mean_plaquette3d(ϕ, J)
    1 - Action.S_SAshift_3d(ϕ, J) / (length(ϕ) * J)
end

function mean_plaquette2d(ϕ, J)

    1 - Action.S_SAshift_2d(ϕ, J) / (length(ϕ) * J)

end

const unit_vectors_3d = [I[1:3, k] for k in 1:3]
# function loop3d_old(ϕ, μ, ν, rμ, rν, n)
#     if μ == ν
#         return 0
#     end

#     #can use mod1 to wrap indices

#     μp_indices = [n + unit_vectors_3d[μ] * r for r in 0:rμ-1]
#     μn_indices = [n + unit_vectors_3d[ν] * (rν - 1) + unit_vectors_3d[μ] * r for r in (rμ-1):-1:0]
#     νp_indices = [n + unit_vectors_3d[μ] * (rμ - 1) + unit_vectors_3d[ν] * r for r in 0:rν-1]
#     νn_indices = [n + unit_vectors_3d[ν] * r for r in rν-1:-1:0]

#     loop = sum(ϕ[μ, ind...] for ind in μp_indices) +
#            sum(ϕ[ν, ind...] for ind in νp_indices) +
#            sum(-ϕ[μ, ind...] for ind in μn_indices) +
#            sum(-ϕ[ν, ind...] for ind in νn_indices)
#     loop
# end

function loop3d(ϕ, μ, ν, rμ, rν, n)
    if μ == ν
        return 0
    end
    size_μ = size(ϕ)[μ+1]
    size_ν = size(ϕ)[ν+1]

    α = 6 - μ - ν # the index perp to loop

    μp_indices = mod1.(n[μ] .+ (0:rμ-1), size_μ)
    νp_indices = mod1.(n[ν] .+ (0:rν-1), size_ν)
    #can reverse order in last two terms but since we are just adding it doesnt matter
    loop = sum(ϕ[μ, (unit_vectors_3d[α] * n[α] + unit_vectors_3d[μ] * ind + unit_vectors_3d[ν] * νp_indices[1])...] for ind in μp_indices) +
           sum(ϕ[ν, (unit_vectors_3d[α] * n[α] + unit_vectors_3d[μ] * μp_indices[end] + unit_vectors_3d[ν] * ind)...] for ind in νp_indices) +
           sum(-ϕ[μ, (unit_vectors_3d[α] * n[α] + unit_vectors_3d[μ] * ind + unit_vectors_3d[ν] * νp_indices[end])...] for ind in μp_indices) +
           sum(-ϕ[ν, (unit_vectors_3d[α] * n[α] + unit_vectors_3d[μ] * μp_indices[1] + unit_vectors_3d[ν] * ind)...] for ind in νp_indices)

    loop
end

function benchmark_loops(n_runs)
    ϕ = randn((3, 16, 16, 16))
    tests = []
    for _ in n_runs
        μ = rand([1, 2, 3])
        ν = rand([1, 2, 3])
        rμ = rand(1:16)
        rν = rand(1:16)
        n = [rand(1:16), rand(1:16), rand(1:16)]
        n = [1, 1, 1]
        push!(tests, (μ, ν, rμ, rν, n))
    end
    # @time map(t -> loop3d(ϕ, t...), tests)
    # @time map(t -> loop3d(ϕ, t...), tests)
    @time map(t -> loop3d_v2(ϕ, t...), tests)
    @time map(t -> loop3d_v2(ϕ, t...), tests)
end

# function optimal_hmc_step(ϕ_init, S, dSdϕ, S_args)

# end


function acceptance_rate(ϕ, update, update_args, steps, initial_steps)
    accept_count = 0
    for _ in 1:initial_steps
        ϕ = update(ϕ, update_args...)[1]
    end
    for _ in 1:steps
        ϕ, ΔH, accepted = update(ϕ, update_args...)

        accept_count += accepted
    end
    accept_count / steps
end


function measure(filename, mode, ϕ_init, initial_steps, measurement_interval, Nmeasurements, update, update_args, measurement_functions, measurement_args, measurement_info)
    # assumes the first arguement of update is ϕ
    ϕ = ϕ_init
    for _ in 1:initial_steps
        ϕ = update(ϕ, update_args...)[1]
    end

    accept_count = 0
    open("measurements/$filename", mode) do io
        write(io, "# $Nmeasurements\n")
        for ith_measurement in 1:Nmeasurements
            # print(ith_measurement)
            for _ in 1:measurement_interval

                ϕ, ΔH, accepted = update(ϕ, update_args...)
                accept_count += accepted
            end


            measurement = map((f, args) -> f(ϕ, args...), measurement_functions, measurement_args)
            writedlm(io, [measurement_info... measurement...])
        end
    end
    accept_rate = accept_count / (Nmeasurements * measurement_interval)
    # write(file, "final_config", ϕ)
    DataHandler.save_field("final_fields/$filename", ϕ)


    print("acceptance rate: $accept_rate\n")
    # end
    #save field for later con
end





# dataβ = [1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0]
# dataβ_full = [0.2, 0.4, 0.6, 0.8, 1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]

# dataP = [0.475, 0.629, 0.656, 0.704, 0.748, 0.790, 0.806, 0.834, 0.854, 0.869, 0.881]


# for β in dataβ_full
#     measure("test_full.txt", "a", zeros(Float64, (3, 16, 16, 16)), 100, 10, 20, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# end

# run16 = [0.614, 0.732, 0.744, 0.770, 0.793, 0.816, 0.827, 0.846, 0.862, 0.876, 0.887]
# run16 = [0.475, 0.630, 0.654, 0.704, 0.749, 0.792, 0.807, 0.836, 0.854, 0.869, 0.881]
# run16_wider = [-0.601, -0.232, 0.128, 0.440, 0.895, 0.903, 0.909, 0.913, 0.919]
# dataβ_full = [0.2, 0.4, 0.6, 0.8, 1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]

# run16_full = [-0.601, -0.232, 0.128, 0.440, 0.614, 0.732, 0.744, 0.770, 0.793, 0.816, 0.827, 0.846, 0.862, 0.876, 0.887, 0.895, 0.903, 0.909, 0.913, 0.919]
# dataβ_wider = [0.2, 0.4, 0.6, 0.8, 3.2, 3.4, 3.6, 3.8, 4.0]


# half16 = [-0.057,0.253,0.307,0.416,0.493,0.581,0.617,0.670,0.711,0.738,0.762]
# p = Plots.plot(dataβ, dataP, label="Hamer")
# p = Plots.plot!(p, dataβ, run16, label="Mine")
# p = Plots.plot!(dataβ[1:4], dataβ[1:4] / 2)
# p = Plots.xlabel!(p, "J")
# p = Plots.ylabel!(p, "<P>")


# p |> display
# β = 1
# measure("autocorrelation_1_05_1000.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.93

# measure("autocorrelation_1_05_100.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 100, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.91

# measure("autocorrelation_1_1_100.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 100, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.1, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.77

# β = 0.5
# measure("autocorrelation_05_05_1000.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.977

# β = 3
# measure("autocorrelation_3_05_1000.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.737


# β = 1
# measure("autocorrelation_1_05_1000_1.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.946

# β = 1
# measure("autocorrelation_1_05_1000_2.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.939

# β = 1
# measure("autocorrelation_1_05_10000.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 10000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# 0.9343


β = 1
Δτ = 0.115
lf_steps = Int(round(1 / Δτ))
# acceptance_rate(zeros(Float64, (3, 16, 16, 16)), HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), 200, 100)
# measure("plaquette_v_trajectory_cold.txt", "w", zeros(Float64, (3, 16, 16, 16)), 0, 1, 150, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
β = 1
Δτ = 0.05
lf_steps = Int(round(1 / Δτ))
# acceptance_rate((rand(Float64,(3,16,16,16))*2 .-1)*pi, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), 200, 100)
# measure("plaquette_v_trajectory_hot.txt", "w", (rand(Float64, (3, 16, 16, 16)) * 2 .- 1) * pi, 0, 1, 150, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))

β = 1
Δτ = 0.115
lf_steps = Int(round(1 / Δτ))
for _ in 1:10
    measure("plaquette_heat_capacity.txt", "a", (zeros(Float64, (3, 16, 16, 16)) * 2 .- 1) * pi, 250, 10, 10, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
end
