include("action.jl")
include("hmc.jl")
using Statistics
import Plots


using JLD
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


function measure(file_name, ϕ_init, initial_steps, measurement_interval, Nmeasurements, update, update_args, measurement_functions, measurement_args)
    # temporarily disabled saving to file and just output mean values
    ϕ = ϕ_init
    for _ in 1:initial_steps
        ϕ = update(ϕ, update_args...)[1]
    end
    testres = []
    accept_count = 0
    # jldopen("measurements/$file_name.jdl", "w") do file
    for ith_measurement in 1:Nmeasurements
        # print(ith_measurement)
        for _ in 1:measurement_interval

            ϕ, ΔH, accepted = update(ϕ, update_args...)
            accept_count += accepted
        end


        m1 = mean_plaquette3d(ϕ, measurement_args[1][1])

        push!(testres, m1)

        measurement = map(fargs_pair -> fargs_pair[1](ϕ, fargs_pair[2]...), collect(zip(measurement_functions, measurement_args)))
        # write(file, "measurement_$ith_measurement", measurement)
    end
    accept_rate = accept_count / (Nmeasurements * measurement_interval)
    # write(file, "final_config", ϕ)
    print("mean\n")
    print(mean(testres))
    print("\n")
    print("acceptance rate: $accept_rate\n")
    # end
    #save field for later con
end





dataβ = [1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0]
dataP = [0.475, 0.629, 0.656, 0.704, 0.748, 0.790, 0.806, 0.834, 0.854, 0.869, 0.881]


for β in dataβ
    measure("test3d-$β", zeros(Float64, (3, 16, 16, 16)), 100, 100, 10, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),))
end

run16 = [0.614, 0.732, 0.744, 0.770, 0.793, 0.816, 0.827, 0.846, 0.862, 0.876, 0.887]

p = Plots.plot(dataβ, dataP, label="Hamer")
p = Plots.plot!(p, dataβ, run16, label="Mine")
p = Plots.plot!(dataβ[1:4], dataβ[1:4] / 2)
p = Plots.xlabel!(p, "J")
p = Plots.ylabel!(p, "<P>")


p |> display
