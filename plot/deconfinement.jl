include("../measurement.jl")
include("../measurement_functions.jl")
include("../u_1_lgt.jl")
include("../hmc.jl")
include("../data_handler.jl")
include("../autocorrelation.jl")
include("../aic.jl")

βs = [2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0]
a_range = 0:100
data_range = [2]
wilson_lines = []
errors = []
for β in βs
    filename = "autocorrelation_analysis/wilson_line_abs_fa0.99_β$(β)_τ0.6_2.txt"
    data = DataHandler.load_data_reshaped(filename)
    # display(data)
    push!(wilson_lines, data[1,1,2])
    push!(errors, data[1,2,2])
end

plot(βs, wilson_lines, yerrors=errors) |> display
