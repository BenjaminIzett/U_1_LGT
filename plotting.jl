using Plots
using CurveFit
include("data_handler.jl")




function plot_autocorrelation(filename, fit_cutoff; error_file="")
    data = DataHandler.load_data(filename)

    if error_file != ""
        errors = DataHandler.load_data(error_file)
        plt = plot(data[2:end], ribbon=errors[2:end], xlabel="Lag time", ylabel="Autocorrelation time", label="data")

    else
        plt = plot(data[2:end], xlabel="Lag time", ylabel="Autocorrelation time", label="data")
    end
    xs = 1:fit_cutoff
    plot_xs = 0:40
    fit = linear_fit(xs, log.(data[2:fit_cutoff+1]))
    plt = plot!(plt, plot_xs, exp.(fit[1] .+ plot_xs .* fit[2]), label="Exponential fit, exponent: $(round(fit[2],digits=3))")
end
function plot_autocorrelation_log(filename, linear_cutoff; error_file="")
    data = DataHandler.load_data(filename)[1, :]
    xs = collect(1:linear_cutoff)
    # fit = log_fit(xs, data[2:linear_cutoff])
    # display(log.(data[2:linear_cutoff]))
    fit = linear_fit(xs, log.(data[2:linear_cutoff+1]))
    display(fit)
    positive_cutoff = findfirst(n -> n <= 0, data)
    plt = plot(data[2:positive_cutoff-1], yaxis=:log, xlabel="Lag time", ylabel="Autocorrelation time", label="data")
    # display(xs)
    # display(fit[1] .+ xs .* fit[2])
    plt = plot!(plt, xs, exp.(fit[1] .+ xs .* fit[2]), yaxis=:log, label="Linear fit, slope: $(round(fit[2],digits=3))")
    # display(0 .+ xs * m)
    plt

end

function plot_av_plaquette_v_trajectory(filename)
    data = DataHandler.load_data(filename)

    plt = plot(data[:, 3])
    xaxis!(plt, "Trajectory")
    yaxis!(plt, "<P>")
    plt
end
function plt_hot_cold_plaquette_v_trajectory(filename_hot, filename_cold)
    data_hot = DataHandler.load_data(filename_hot)
    data_cold = DataHandler.load_data(filename_cold)
    plt = plot(data_hot[:, 3], label="Random initial field")
    plt = plot!(plt, data_cold[:, 3], label="Zero initial field")
    xaxis!(plt, "Trajectory")
    yaxis!(plt, "<P>")

    plt
end

function plot_log_ratio_wilson_loops(W_τ)
    ratios = log.(W_τ[2:end] ./ W_τ[1:end-1])
    τs = collect(1:7)
    plot(τs, ratios, ylims=(-0.55, -0.3), marker=:x)
end

check_hamer("analysis/test_a.txt", false) |> display
check_hamer("analysis/fa_hamer_p_4.txt", false) |> display
#
check_hamer("analysis/test_full.txt", true) |> display
plot_autocorrelation("analysis/autocorrelation_1_05_1000.txt") |> display
plot_autocorrelation("analysis/autocorrelation_1_05_100.txt") |> display
plot_autocorrelation("analysis/autocorrelation_1_1_100.txt") |> display

plot_autocorrelation("analysis/autocorrelation_05_05_1000.txt") |> display
plot_autocorrelation("analysis/autocorrelation_3_05_1000.txt") |> display

plot_autocorrelation("analysis/autocorrelation_1_05_1000_1.txt") |> display
plot_autocorrelation("analysis/autocorrelation_1_05_1000_2.txt") |> display

plot_autocorrelation("analysis/autocorrelation_1_05_10000.txt", 10) |> display
plot_autocorrelation_log("analysis/autocorrelation_1_05_10000.txt", 10) |> display

plot_autocorrelation("analysis/mean/autocorrelation_mean.txt", 15; error_file="analysis/std/autocorrelation_mean.txt") |> display
plot_autocorrelation_log("analysis/mean/autocorrelation_mean.txt", 15; error_file="analysis/std/autocorrelation_mean.txt") |> display


plt_hot_cold_plaquette_v_trajectory("measurements/plaquette_v_trajectory_hot.txt", "measurements/plaquette_v_trajectory_cold.txt")
plot_log_ratio_wilson_loops([0.47502213651053316, 0.28598744547295896, 0.18426005836045833, 0.12132245839749849, 0.08053928509021227, 0.05372076844599194, 0.03574225533086559, 0.023890720866362482]) |> display

plot_log_ratio_wilson_loops([0.652689865407466, 0.43460440810953, 0.290419522053951, 0.194719179021082, 0.130881557501511, 0.0882875282215217, 0.0597037533011341, 0.0402920138185977]) |> display

plot(1:7,[-0.413488245951787,-0.409808840125972,-0.40784679677955,-0.405833883761059,-0.405319679215287,-0.406122267220057,-0.409624139571279],yerror=[0.00279108404905521,0.0030849418757994,0.00367971568560996,0.00409957217941803,0.00423293355370217,0.00674138862741733,0.00837743546333262],ylims=(-0.55, -0.3))
#
plot_autocorrelation("analysis/autocorrelation/no_fa_plaquette_test_1.txt",10) |> display
plot_autocorrelation("analysis/autocorrelation/with_fa_plaquette_test_2.txt",10) |> display
