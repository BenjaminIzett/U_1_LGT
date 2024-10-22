include("action.jl")
include("hmc.jl")
include("data_handler.jl")

using Statistics
import Plots
using DelimitedFiles


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


function measure(filename, mode, ϕ_init, initial_steps, measurement_interval, Nmeasurements, update, update_args, measurement_functions, measurement_args, measurement_info)
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





dataβ = [1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0]
dataP = [0.475, 0.629, 0.656, 0.704, 0.748, 0.790, 0.806, 0.834, 0.854, 0.869, 0.881]


for β in dataβ
    measure("test_a.txt", "a", zeros(Float64, (3, 16, 16, 16)), 100, 5, 20, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
end

# run16 = [0.614, 0.732, 0.744, 0.770, 0.793, 0.816, 0.827, 0.846, 0.862, 0.876, 0.887]
run16 = [0.475, 0.630, 0.654, 0.704, 0.749, 0.792, 0.807, 0.836, 0.854, 0.869, 0.881]
run16_wider = [-0.601, -0.232, 0.128, 0.440, 0.895, 0.903, 0.909, 0.913, 0.919]
dataβ_full = [0.2, 0.4, 0.6, 0.8, 1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]

run16_full = [-0.601, -0.232, 0.128, 0.440, 0.614, 0.732, 0.744, 0.770, 0.793, 0.816, 0.827, 0.846, 0.862, 0.876, 0.887, 0.895, 0.903, 0.909, 0.913, 0.919]
dataβ_wider = [0.2, 0.4, 0.6, 0.8, 3.2, 3.4, 3.6, 3.8, 4.0]


# half16 = [-0.057,0.253,0.307,0.416,0.493,0.581,0.617,0.670,0.711,0.738,0.762]
p = Plots.plot(dataβ, dataP, label="Hamer")
p = Plots.plot!(p, dataβ, run16, label="Mine")
p = Plots.plot!(dataβ[1:4], dataβ[1:4] / 2)
p = Plots.xlabel!(p, "J")
p = Plots.ylabel!(p, "<P>")


p |> display
