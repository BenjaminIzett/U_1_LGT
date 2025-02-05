include("../data_handler.jl")
using Statistics
using Plots

filename = "teper_fa099_β22.txt"
raw_data = DataHandler.load_data_reshaped("measurements/$filename")

data = raw_data[:, :, 5:end]

mean_data = mean(data, dims=2)

ratio = mean_data[:, 1, 2:end] ./ mean_data[:, 1, 1:end-1]
ln_ratio = -log.(ratio)
# plot(2:8,ln_ratio[1,:])|>display
Vs = mean(ln_ratio, dims=2)[:, 1]
plot(2:8, -log.(ratio[3, :])) |> display
ls = [10, 14, 18] .* sqrt(0.020581)
display(Vs)
plot(ls, Vs ./ sqrt(0.020581)) |> display
