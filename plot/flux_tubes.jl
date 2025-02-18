include("../data_handler.jl")
include("../aic.jl")
using Statistics
using Plots
using CurveFit

filename = "teper_fa099_β22_5.txt"
raw_data = DataHandler.load_data_reshaped("autocorrelation_analysis/$filename")

mean_data = raw_data[:, 1, :]

errors = raw_data[:, 2, :]


ratio = mean_data[:, 2:end] ./ mean_data[:, 1:end-1]

ratios = []
ratio_errors = []
(N, M) = size(mean_data)
for n in 1:N
    endpoint = findfirst(t -> (t <= 0), mean_data[n, :])
    if endpoint == nothing
        endpoint = M
    end
    push!(ratios, -log.(mean_data[n, 2:(endpoint-1)] ./ mean_data[n, 1:(endpoint-2)]))
    push!(ratio_errors, sqrt.((errors[n, 2:(endpoint-1)] ./ mean_data[n, 2:(endpoint-1)]) .^ 2 .+ (errors[n, 1:(endpoint-2)] ./ mean_data[n, 1:(endpoint-2)]) .^ 2))
end

K = 5
# display(ratios[K][15])
testfit = AIC.fit_window(Polynomial, (0,), 1:length(ratios[K]), ratios[K], ratio_errors[K], 1)
plt = plot(1:length(ratios[K]), ratios[K], yerrors=ratio_errors[K])
plot!(plt, 1:length(ratios[K]), [testfit(4) for _ in 1:length(ratios[K])])
plt |> display

Vs = [AIC.fit_window(Polynomial, (0,), 1:length(ratios[n]), ratios[n], [1 for _ in 1:length(ratios[n])], 1)(5) for n in 1:N]

ls = [8, 10, 14, 18, 22, 26, 34,]# 38, 42, 46, 50, 68]
σ = 0.027423
display(Vs)
plt = scatter(ls * sqrt(σ), Vs / sqrt(σ), xlims=(0, 12), ylims=(0, 14))



tevs = [
    1.1460832745236407, 1.3239237826393797, 3.4975299929428365, 4.179251940719832, 5.513055751587861, 6.224417784050812, 6.935779816513761, 7.528581510232886, 8.190543401552576, 11.352152434721242, 2.0846859562455897, 2.7960479887085388]
tels = [
    1.4836363636363632,
    1.6581818181818178,
    3.6436363636363627,
    4.298181818181816,
    5.629090909090906,
    6.283636363636361,
    6.949090909090906,
    7.6145454545454525,
    8.269090909090908,
    11.258181818181818,
    2.3236363636363633,
    2.9890909090909084]
plt = scatter!(plt, tels, tevs, marker=:x)
plt |> display
display(readline())
