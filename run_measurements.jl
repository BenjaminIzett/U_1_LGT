include("measurement.jl")
include("measurement_functions.jl")
include("u_1_lgt.jl")
include("hmc.jl")


# dataβ = [1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0]
# dataβ_full = [0.2, 0.4, 0.6, 0.8, 1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]

# dataP = [0.475, 0.629, 0.656, 0.704, 0.748, 0.790, 0.806, 0.834, 0.854, 0.869, 0.881]


# for β in dataβ_full
#     measure("test_full.txt", "a", zeros(Float64, (3, 16, 16, 16)), 100, 10, 20, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# end

# run16 = [0.614, 0.732, 0.744, 0.770, 0.793, 0.816, 0.827, 0.846, 0.862, 0.876, 0.887]
# run16 = [0.475, 0.630, 0.654, 0.704, 0.749, 0.792, 0.807, 0.836, 0.854, 0.869, 0.881]
# run16_wider = [-0.601, -0.232, 0.128, 0.440, 0.895, 0.903, 0.909, 0.913, 0.919]
# dataβ_full = [0.2, 0.4, 0.6, 0.8, 1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]

# run16_full = [-0.601, -0.232, 0.128, 0.440, 0.614, 0.732, 0.744, 0.770, 0.793, 0.816, 0.827, 0.846, 0.862, 0.876, 0.887, 0.895, 0.903, 0.909, 0.913, 0.919]
# dataβ_wider = [0.2, 0.4, 0.6, 0.8, 3.2, 3.4, 3.6, 3.8, 4.0]


# half16 = [-0.057, 0.253, 0.307, 0.416, 0.493, 0.581, 0.617, 0.670, 0.711, 0.738, 0.762]
# p = Plots.plot(dataβ, dataP, label="Hamer")
# p = Plots.plot!(p, dataβ, run16, label="Mine")
# p = Plots.plot!(dataβ[1:4], dataβ[1:4] / 2)
# p = Plots.xlabel!(p, "β")
# p = Plots.ylabel!(p, "<P>")


# p |> display
# β = 1
# measure("autocorrelation_1_05_1000.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# 0.93

# measure("autocorrelation_1_05_100.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 100, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# 0.91

# measure("autocorrelation_1_1_100.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 100, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.1, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# 0.77

# β = 0.5
# measure("autocorrelation_05_05_1000.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# 0.977

# β = 3
# measure("autocorrelation_3_05_1000.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# 0.737


# β = 1
# measure("autocorrelation_1_05_1000_1.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# 0.946

# β = 1
# measure("autocorrelation_1_05_1000_2.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# 0.939

# β = 1
# measure("autocorrelation_1_05_10000.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 10000, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# 0.9343
# #

# β = 2
# Δτ = 0.08
# lf_steps = Int(round(1 / Δτ))
# display(acceptance_rate(zeros(Float64, (3, 16, 16, 16)), HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, lf_steps, Δτ, (β,)), 200, 100))
# measure("plaquette_v_trajectory_cold.txt", "w", zeros(Float64, (3, 16, 16, 16)), 0, 1, 150, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))

# β = 1
# Δτ = 0.05
# lf_steps = Int(round(1 / Δτ))
# acceptance_rate((rand(Float64, (3, 16, 16, 16)) * 2 .- 1) * pi, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, lf_steps, Δτ, (β,)), 200, 100)
# measure("plaquette_v_trajectory_hot.txt", "w", (rand(Float64, (3, 16, 16, 16)) * 2 .- 1) * pi, 0, 1, 150, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))

# β = 1
# Δτ = 0.115
# lf_steps = Int(round(1 / Δτ))
# for nth_run in 1:10
#     measure("autocorrelation/autocorrelation_$(β)_10000_$(nth_run).txt", "w", zeros(Float64, (3, 16, 16, 16)), 250, 1, 10000, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# end

# β = 1
# Δτ = 0.115
# lf_steps = Int(round(1 / Δτ))
# measure("wilson_loop1.txt", "w", zeros(Float64, (3, 16, 16, 16)), 250, 10, 2, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, lf_steps, Δτ, (β,)), (range_loop3d,), ((1:8, 1:8,),), (β, 16))

# β = 2
# R = 5
# Δτ = 0.08
# lf_steps = Int(round(1 / Δτ))
# for _ in 1:10
#     measure("wilson_loop_run_B2_R5_ot_1000_10_real.txt", "a", zeros(Float64, (3, 16, 16, 16)), 250, 10, 1000, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, lf_steps, Δτ, (β,)), (rect_wilson_loop_3d, rect_wilson_loop_3d, rect_wilson_loop_3d, rect_wilson_loop_3d, rect_wilson_loop_3d, rect_wilson_loop_3d, rect_wilson_loop_3d, rect_wilson_loop_3d), ((R, 1, β, 0.7, 10), (R, 2, β, 0.7, 10), (R, 3, β, 0.7, 10), (R, 4, β, 0.7, 10), (R, 5, β, 0.7, 10), (R, 6, β, 0.7, 10), (R, 7, β, 0.7, 10), (R, 8, β, 0.7, 10)), (β, 16))
# end

# β = 2
# Δτ = 0.08
# lf_steps = Int(round(1 / Δτ))
# measurement_functions = (rect_wilson_loop_3d for R in 1:10 for τ in 1:10)
# measurement_function_arguments = ((R, τ, β, 0.7, 10) for R in 1:10 for τ in 1:10)
# for _ in 1:10
#     measure("wilson_loop_run_long_thermal_range_ot_110_1500_10_real.txt", "a", zeros(Float64, (3, 16, 16, 16)), 1000, 10, 1500, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, lf_steps, Δτ, (β,)), measurement_functions, measurement_function_arguments, (β, 16))
# end

# β = 2
# Δτ = 0.05
# lf_steps = Int(round(1 / Δτ))
# measurement_functions = (parisi_3d_loops for _ in 1:36)
# measurement_function_arguments = ((R, τ, β) for R in 3:10 for τ in R:10)
# for _ in 1:10
#     measure("wl_check_p_2_1000_1000.txt", "a", zeros(Float64, (3, 32, 32, 32)), 1000, 10, 1000, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, lf_steps, Δτ, (β,)), measurement_functions, measurement_function_arguments, (β, 32))
# end

# β = 2
# Δτ = 0.05
# lf_steps = Int(round(1 / Δτ))
# measurement_functions = (rect_wilson_loop_3d for R in 2:10 for τ in 2:10)
# measurement_function_arguments = ((R, τ, β, 0.7, 10) for R in 2:10 for τ in 2:10)
# for _ in 1:10
#     measure("wilson_loop32_run_range_no_1000_10_real.txt", "a", zeros(Float64, (3, 16, 16, 16)), 250, 10, 1000, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, lf_steps, Δτ, (β,)), measurement_functions, measurement_function_arguments, (β, 16))
# end

# acceptance_rate(zeros(Float64, (3, 32, 32, 32)), HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, lf_steps, Δτ, (β,)), 500, 250)

# β = 3
# Δτ = 0.125
# lf_steps = Int(round(1 / Δτ))
# inverse_FK = HMC.inv_FK_3d_mass(16, 16, 16, 10000)
# acceptance_rate(zeros(Float64, (3, 16, 16, 16)), HMC.fa_hmc_run_3d, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, inverse_FK, lf_steps, Δτ, (β,)), 500, 250)

# dataβ = [1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0]
# dataβ_full = [0.2, 0.4, 0.6, 0.8, 1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]

# dataP = [0.475, 0.629, 0.656, 0.704, 0.748, 0.790, 0.806, 0.834, 0.854, 0.869, 0.881]


# for β in dataβ
#     measure("fa_hamer_p_4.txt", "a", zeros(Float64, (3, 16, 16, 16)), 250, 10, 100, HMC.fa_hmc_run_3d, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, inverse_FK, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# end

# β = 2
# Δτ = 0.1
# lf_steps = Int(round(1 / Δτ))
# inverse_FK = HMC.inv_FK_3d_kappa(16, 16, 16, 0.5)
# acceptance_rate(zeros(Float64, (3, 16, 16, 16)), HMC.fa_hmc_run_3d, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, inverse_FK, lf_steps, Δτ, (β,)), 500, 250)
# opt = optimise_update_args(zeros(Float64, (3, 16, 16, 16)), HMC.fa_hmc_run_3d, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, inverse_FK, lf_steps, Δτ, (β,)), 4, 5, 250, 250)

# measurement_functions = (rect_wilson_loop_3d for R in 1:10 for τ in 1:10)
# measurement_function_arguments = ((R, τ, β, 0.7, 10) for R in 1:10 for τ in 1:10)
# for _ in 1:10
#     measure("fa_wilson_loop_5.txt", "a", zeros(Float64, (3, 16, 16, 16)), 500, 10, 500, HMC.fa_hmc_run_3d, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, inverse_FK, lf_steps, Δτ, (β,)), measurement_functions, measurement_function_arguments, (β, 16))
# end


# @time measure("timetest.txt", "a", zeros(Float64, (3, 16, 16, 16)), 250, 10, 1000, HMC.fa_hmc_run_3d, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, inverse_FK, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# @time measure("timetest.txt", "a", zeros(Float64, (3, 16, 16, 16)), 250, 10, 1000, HMC.fa_hmc_run_3d, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, inverse_FK, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))

# @time measure("timetest.txt", "a", zeros(Float64, (3, 16, 16, 16)), 250, 10, 1000, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# @time measure("timetest.txt", "a", zeros(Float64, (3, 16, 16, 16)), 250, 10, 1000, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))


# β = 1
# Δτ = 0.1
# lf_steps = Int(round(1 / Δτ))
# inverse_FK = HMC.inv_FK_3d(16, 16, 16, 0.01)
# measure("autocorrelation/no_fa_plaquette_test_1.txt", "w", zeros(Float64, (3, 16, 16, 16)), 250, 1, 10000, HMC.hmc_run, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# β = 1
# Δτ = 0.15
# lf_steps = Int(round(1 / Δτ))
# measure("autocorrelation/with_fa_plaquette_test_2.txt", "w", zeros(Float64, (3, 16, 16, 16)), 250, 1, 10000, HMC.fa_hmc_run_3d, (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, inverse_FK, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))



β = 2

update = HMC.hmc_run
update_args = (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.1, (β,))
lf_index = 3
Δτ_index = 4

# update = HMC.fa_hmc_run_3d
# κ = 0.25
# inverse_FK = HMC.inv_FK_3d_kappa(16, 16, 16, κ)
# update_args = (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.1, inverse_FK, (β,))
# lf_index = 3
# Δτ_index = 4
# measurement_functions = (MeasurementFunctions.parisi_loop_range_3d,)
# measurement_args = ((β, [(3, 3), (3, 7), (7, 7)]),)
# measurement_info = (β, 16)
# Measurement.measure("loop_autocorrelation_1.txt", "a", U_1_LGT.zero_lattice_3d(16, 16, 16), 250, 1, 1000, update, update_args, measurement_functions, measurement_args, measurement_info)

# update = HMC.fa_hmc_run_3d
# κ = 0.5
# inverse_FK = HMC.inv_FK_3d_kappa(16, 16, 16, κ)
# update_args = (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.1, inverse_FK, (β,))
# lf_index = 3
# Δτ_index = 4
# measurement_functions = (MeasurementFunctions.hamer_loop_range_3d,)
measurement_functions = (MeasurementFunctions.mean_plaquette_3d,)
measurement_args = ((β,),)
measurement_info = (β, 16)
# Measurement.repeat_optimise_measure("hamer_check1.txt", U_1_LGT.zero_lattice_3d(16, 16, 16), 1000, 10, 1000, update, update_args, measurement_functions, measurement_args, measurement_info, lf_index, Δτ_index, 250, 1000, 10)

# optimised_args = Measurement.optimise_update_args(U_1_LGT.zero_lattice_3d(16, 16, 16), update, update_args, lf_index, Δτ_index, 250, 1000)
# Measurement.measure("timecheck.txt", "w", U_1_LGT.zero_lattice_3d(16, 16, 16), 200, 10, 5, update, update_args, measurement_functions, measurement_args, measurement_info)
# @time Measurement.measure("timecheck.txt", "w", U_1_LGT.zero_lattice_3d(16, 16, 16), 1000, 1000, 5, update, update_args, measurement_functions, measurement_args, measurement_info)

#fa:196s vs 160s 20% slower

ϕ_init = U_1_LGT.rand_lattice_3d(16, 16, 16)

random_update_args = (U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 1, 0.1, (β,))
ϕ = Measurement.initialise_random_start(ϕ_init, update, random_update_args, 1000)
ϕ = Measurement.initialise_random_start(ϕ, update, random_update_args, 1000)
# for _ in 1000
#     ϕ = HMC.hmc_run_noMS(ϕ, U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 1, 0.1, (β,))
# end
Measurement.repeat_optimise_measure("hamer_check_random_start_1.txt", ϕ, 1000, 10, 1000, update, update_args, measurement_functions, measurement_args, measurement_info, lf_index, Δτ_index, 250, 1000, 10)
