module Measurement
include("u_1_lgt.jl")
include("hmc.jl")
include("data_handler.jl")

using Statistics
import Plots
using DelimitedFiles
using LinearAlgebra

import ShiftedArrays as sa
import SpecialFunctions as sf
#initial field
#update field
#steps_to_thermalise
#measurement interval
#number of measurements
#list of measurent functions
#update_field args
#measurement args

function acceptance_rate(ϕ_init, update, update_args, initial_steps, measurement_steps)
    accept_count = 0
    ϕ = ϕ_init
    for _ in 1:initial_steps
        ϕ = update(ϕ, update_args...)[1]
    end
    for _ in 1:measurement_steps
        ϕ, ΔH, accepted = update(ϕ, update_args...)

        accept_count += accepted
    end
    rate = accept_count / measurement_steps
end

function initialise_random_start(ϕ_init, update, update_args, steps)
    accept_count = 0
    ϕ = ϕ_init
    for _ in 1:steps
        ϕ, ΔH, accepted = update(ϕ, update_args...)
        accept_count += accepted
    end

    rate = accept_count / steps
    display("Random start accept rate: $rate")
    return ϕ
end


function optimise_update_args(ϕ, update, update_args, lf_index, Δτ_index, initial_steps, measurement_steps)
    # Determines the optimal number of leapfrog steps (lf_steps) and thus the time step (Δτ)
    # under the simplifying constraint that lf_steps * Δτ = 1. (Not always optimal)

    # Should probably change update_args to a dictionary instead of using indices
    min_rate = 0.625
    optimal_rate = 0.65
    max_rate = 0.75
    max_count = 10

    lf_steps_initial = max(1, update_args[lf_index])
    optimised_update_args = [i for i in update_args]
    optimised_update_args[lf_index] = lf_steps_initial
    optimised_update_args[Δτ_index] = 1 / lf_steps_initial
    rate = acceptance_rate(ϕ, update, optimised_update_args, initial_steps, measurement_steps)

    flag = false
    itr_count = 0
    while !flag && (itr_count < max_count) && ((rate < min_rate) | (rate > max_rate))

        new_lf_steps = Int(round(optimised_update_args[lf_index] * (1 + optimal_rate - rate)))

        if new_lf_steps == optimised_update_args[lf_index]
            if rate < min_rate
                if new_lf_steps > 1
                    new_lf_steps -= 1
                else
                    flag = true
                end
            else
                new_lf_steps += 1
            end
        end
        optimised_update_args[lf_index] = new_lf_steps
        optimised_update_args[Δτ_index] = 1 / new_lf_steps

        rate = acceptance_rate(ϕ, update, optimised_update_args, measurement_steps, initial_steps)
        itr_count += 1
    end
    if flag || (itr_count >= max_count)
        display("Failed to find optimised parameters, rate: $rate, lf_steps: $(optimised_update_args[lf_index]), itr_count: $itr_count.")
    else
        display("Successfully found optimised parameters, rate: $rate, lf_steps: $(optimised_update_args[lf_index]), itr_count: $itr_count.")
    end
    return tuple(optimised_update_args...)
end


function measure(filename, mode, ϕ_init, initial_steps, measurement_interval, N_measurements, update, update_args, measurement_functions, measurement_args, measurement_info)
    # assumes the first arguement of update is ϕ
    ϕ = ϕ_init
    for _ in 1:initial_steps
        ϕ = update(ϕ, update_args...)[1]
    end

    accept_count = 0
    open("measurements/$filename", mode) do io
        write(io, "# $N_measurements\n")
        for ith_measurement in 1:N_measurements
            # print(ith_measurement)
            for _ in 1:measurement_interval

                ϕ, ΔH, accepted = update(ϕ, update_args...)
                accept_count += accepted
            end

            # print(measurement_functions, measurement_args)
            measurement = map((f, args) -> f(ϕ, args...), measurement_functions, measurement_args)
            # display(measurement)
            # writedlm(io, [measurement_info... Iterators.flatten(measurement...)...])
            # writedlm(io, transpose([measurement_info..., measurement...]))
            writedlm(io, transpose([measurement_info..., Iterators.flatten(measurement)...]))
        end
    end
    accept_rate = accept_count / (N_measurements * measurement_interval)
    # write(file, "final_config", ϕ)
    DataHandler.save_field("final_fields/$filename", ϕ)


    print("acceptance rate: $accept_rate\n")
    # end
    #save field for later con
end

function optimised_measure(filename, mode, ϕ_init, initial_steps, measurement_interval, N_measurements, update, update_args, measurement_functions, measurement_args, measurement_info, lf_index, Δτ_index, optimise_initial_steps, optimise_measurement_steps)
    optimised_update_args = optimise_update_args(ϕ_init, update, update_args, lf_index, Δτ_index, optimise_initial_steps, optimise_measurement_steps)
    measure(filename, mode, ϕ_init, initial_steps, measurement_interval, N_measurements, update, optimised_update_args, measurement_functions, measurement_args, measurement_info)
end

function repeat_optimise_measure(filename, ϕ_init, initial_steps, measurement_interval, N_measurements, update, update_args, measurement_functions, measurement_args, measurement_info, lf_index, Δτ_index, optimise_initial_steps, optimise_measurement_steps, N_repeats)
    optimised_update_args = optimise_update_args(ϕ_init, update, update_args, lf_index, Δτ_index, optimise_initial_steps, optimise_measurement_steps)
    for _ in 1:N_repeats
        measure(filename, "a", ϕ_init, initial_steps, measurement_interval, N_measurements, update, optimised_update_args, measurement_functions, measurement_args, measurement_info)
    end
end

function measure_range(filename, ϕ_init, initial_steps, measurement_interval, N_measurements, update, update_args_list, measurement_functions, measurement_args_list, measurement_info_list)
    for (update_args, measurement_args, measurement_info) in zip(update_args_list, measurement_args_list, measurement_info_list)
        measure(filename, "a", ϕ_init, initial_steps, measurement_interval, N_measurements, update, update_args, measurement_functions, measurement_args, measurement_info)
    end
end

function optimised_measure_range(filename, ϕ_init, initial_steps, measurement_interval, N_measurements, update, update_args_list, measurement_functions, measurement_args_list, measurement_info_list, lf_index, Δτ_index, optimise_initial_steps, optimise_measurement_steps)
    # Note if lf_step or Δτ are included in measurement_info or (this should not happen) measurement_args they will not be updated
    for (update_args, measurement_args, measurement_info) in zip(update_args_list, measurement_args_list, measurement_info_list)
        optimised_measure(filename, "a", ϕ_init, initial_steps, measurement_interval, N_measurements, update, update_args, measurement_functions, measurement_args, measurement_info, lf_index, Δτ_index, optimise_initial_steps, optimise_measurement_steps)
    end
end

end
