
include("hmc.jl")
include("u_1_lgt.jl")
include("measurement.jl")
include("measurement_functions.jl")

import Plots
using Symbolics
using Statistics

import ShiftedArrays



function test_reversibility(ϕ, dSdϕ, lf_steps, Δτ, S_args)
    p_0 = randn(Float64, size(ϕ))

    ϕ_new, p_new = HMC.leapfrog(ϕ, p_0, dSdϕ, lf_steps, Δτ, S_args)

    ϕ_rev, p_rev = HMC.leapfrog(ϕ_new, -p_new, dSdϕ, lf_steps, Δτ, S_args)

    m_err_p = sum(abs.(p_0 + p_rev)) / length(p_0)
    m_err_ϕ = sum(abs.(ϕ - ϕ_rev)) / length(ϕ)
    @info "Mean Error when reversing leapforg" lf_steps Δτ m_err_ϕ m_err_p
end

function test_derivative_numerical(S, dSdϕ, N, L, shift)
    dims = [L for _ in 1:N]
    ϕ1 = randn((N, dims...))
    S1 = S(ϕ1, 1)
    ns = [rand(1:L, N) for _ in 1:100]
    res = []
    for n in ns
        for μ in 1:N

            δϕ = zeros(size(ϕ1))
            δϕ[μ, n...] = shift
            ϕ2 = ϕ1 + δϕ

            S2 = S(ϕ2, 1)

            push!(res, abs((S2 - S1) / shift - dSdϕ(ϕ1, 1, μ, n)))
        end
    end
    return sum(res) / length(res)
end
function test_derivative_numericalXY(S, dSdϕ, N, L, shift)
    dims = [L for _ in 1:N]
    ϕ1 = randn(dims...)
    S1 = S(ϕ1, 1)
    ns = [rand(1:L, N) for _ in 1:100]
    res = []
    for n in ns


        δϕ = zeros(size(ϕ1))
        δϕ[n...] = shift
        ϕ2 = ϕ1 + δϕ

        S2 = S(ϕ2, 1)
        # print(dSdϕ(ϕ1, 1)[n...])
        # print('\n')
        push!(res, abs((S2 - S1) / shift - dSdϕ(ϕ1, 1)[n...]))

    end
    return sum(res) / length(res)
end


function test_XY_derivative_numerical(S, dSdϕ, N, L, shift)
    dims = [L for _ in 1:N]
    ϕ1 = randn((dims...))
    S1 = S(ϕ1, 1)
    ns = [rand(1:L, N) for _ in 1:100]
    res = []
    for n in ns
        δϕ = zeros(size(ϕ1))
        δϕ[n...] = shift
        ϕ2 = ϕ1 + δϕ

        S2 = S(ϕ2, 1)

        push!(res, abs((S2 - S1) / shift - dSdϕ(ϕ1, 1, n)))

    end
    return sum(res) / length(res)
end
function test_derivative(dSdϕ_s, dSdϕ_t, N, L)
    β = 1

    dims = [L for _ in 1:N]

    ϕ = randn((N, dims...))
    ns = [rand(1:L, N) for _ in 1:100]
    # display([dSdϕ_s(ϕ, β, μ, n) for μ in 1:N for n in ns])
    sum([abs(dSdϕ_s(ϕ, β, μ, n) - dSdϕ_t(ϕ, β)[μ, n...]) for μ in 1:N for n in ns])
end


function test_p_dist(S, dSdϕ, N, L, runs, lf_steps, Δτ, S_args, initial_steps)
    function KE(p)
        sum(p .^ 2) / (2 * length(p)) # <p^2>/2
    end
    ϕ = zeros(Float64, (N, [L for _ in 1:N]...))
    for _ in 1:initial_steps
        ϕ = HMC.hmc_run(ϕ, S, dSdϕ, lf_steps, Δτ, S_args)[1]
    end
    # ϕ_0 = rand(N, [L for _ in 1:N]...)*2π
    initial_KEs = zeros(runs)
    final_KEs = zeros(runs)
    accepted_count = 0
    for run in 1:runs
        p_0 = randn(Float64, (N, [L for _ in 1:N]...))
        initial_KEs[run] = KE(p_0)
        H_init = S(ϕ, S_args...) + 1 / 2 * sum(p_0 .^ 2)

        ϕ_new, p_new = HMC.leapfrog(ϕ, p_0, dSdϕ, lf_steps, Δτ, S_args)
        # display(maximum(abs.(ϕ_new)))
        H_final = S(ϕ_new, S_args...) + 1 / 2 * sum(p_new .^ 2)

        ΔH = H_final - H_init
        r = rand(Float64)

        if r < exp(-ΔH)
            # return ϕ_new, ΔH, true, p_new
            final_KEs[run] = KE(p_new)
            accepted_count += 1
        else
            # return ϕ, ΔH, false, p
            final_KEs[run] = KE(p_0)
        end

        # p_final = HMC.hmc_run(copy(ϕ), S, dSdϕ, lf_steps, Δτ, S_args)[2]
        # final_KEs[run] = KE(p_final)
    end
    accepted_rate = accepted_count / runs
    @info "Final KE average" mean(final_KEs) mean(initial_KEs) runs accepted_rate
end

function test_p_dist_XY(dSdϕ, N, L, runs, lf_steps, Δτ, S_args)
    function KE(p)
        sum(p .^ 2) / (2 * length(p)) # <p^2>/2
    end
    ϕ_0 = zeros(Float64, [L for _ in 1:N]...)
    initial_KEs = zeros(runs)
    final_KEs = zeros(runs)
    for run in 1:runs
        p_0 = randn(Float64, [L for _ in 1:N]...)
        initial_KEs[run] = KE(p_0)
        p_final = HMC.leapfrog(ϕ_0, p_0, dSdϕ, lf_steps, Δτ, S_args)[2]
        final_KEs[run] = KE(p_final)
    end
    @info "Final KE average" mean(final_KEs) mean(initial_KEs) runs
end

function test_XY_p_dist(dSdϕ, N)
    function KE(p)
        sum(p .^ 2) / (2 * length(p)) # <p^2>/2
    end
    p1ds = [randn(Float64, ([L for _ in 1:N]...)) for L in [10, 30, 50, 100, 200, 500, 1000]]

    display([(length(p), KE(p), KE(HMC.leapfrog(zeros(Float64, size(p)), p, dSdϕ, 50, 0.01, (1,))[2])) for p in p1ds])


    # p2ds = [randn(Float64, (2, N, N)) for N in [10, 30, 50, 100, 200, 500, 1000]]

    # [(length(p), KE(p), KE(HMC.leapfrog(zeros(Float64, size(p)), p, Action.XY_dSdϕ_2d, 100, 0.01, (1,))[2])) for p in p2ds]

end

function test_expval(ϕ, S, dSdϕ, lf_steps, Δτ, S_args, HMC_steps, steps_skipped)
    ΔHs = []
    accepted = []
    for k in 2:HMC_steps+1
        ϕ, ΔH_k, accepted_k = HMC.hmc_run(ϕ, S, dSdϕ, lf_steps, Δτ, S_args)
        push!(ΔHs, ΔH_k)
        push!(accepted, accepted_k)
        if k % 100 == 0
            display(k)
        end
    end

    Plots.plot(ΔHs) |> display
    histogram(ΔHs[steps_skipped:end], bins=50, title="Histogram of ΔH", xlabel="ΔH", ylabel="Frequency") |> display
    meanΔHs = mean(ΔHs[steps_skipped:end])

    exp_ΔH = mean(exp.(ΔHs[steps_skipped:end] * -1))
    std_exp_exp_ΔH = std(exp.(ΔHs[steps_skipped:end] * -1))
    acceptance_rate = sum(accepted) / length(accepted)
    @info "Checking the mean value of the exponental of the ΔH" exp_ΔH std_exp_exp_ΔH lf_steps Δτ HMC_steps steps_skipped acceptance_rate meanΔHs
end

function fa_test_expval3d(ϕ, S, dSdϕ, lf_steps, Δτ, S_args, fa_κ, HMC_steps, steps_skipped)
    (D, Lx, Ly, Lz) = size(ϕ)
    inverse_FK = HMC.inv_FK_3d_mass(Lx, Ly, Lz, fa_κ)
    ΔHs = []
    accepted = []
    for k in 2:HMC_steps+1
        ϕ, ΔH_k, accepted_k = HMC.fa_hmc_run_3d(ϕ, S, dSdϕ, lf_steps, Δτ, inverse_FK, S_args)
        push!(ΔHs, ΔH_k)
        push!(accepted, accepted_k)
        if k % 100 == 0
            display(k)
        end
    end

    Plots.plot(ΔHs) |> display
    meanΔHs = mean(ΔHs[steps_skipped:end])
    histogram(ΔHs[steps_skipped:end], bins=50, title="Histogram of ΔH", xlabel="ΔH", ylabel="Frequency") |> display
    exp_ΔH = mean(exp.(ΔHs[steps_skipped:end] * -1))
    std_exp_exp_ΔH = std(exp.(ΔHs[steps_skipped:end] * -1))
    acceptance_rate = sum(accepted) / length(accepted)
    @info "Checking the exponental of the mean value of ΔH" exp_ΔH std_exp_exp_ΔH lf_steps Δτ HMC_steps steps_skipped acceptance_rate meanΔHs
end

function symbolic2_test()
    @variables x111 x112 x121 x122 x211 x212 x221 x222
    ϕ_symb = reshape([x111 x211 x121 x221 x112 x212 x122 x222], (2, 2, 2))
    display(U_1_LGT.S_2d(ϕ_symb, 1))
    display(U_1_LGT.S_2d(ϕ_symb, 1))
    #
    display(U_1_LGT.dSdϕ_2d(ϕ_symb, 1)[2, 2, 1])
    display(U_1_LGT.dSdϕ_2d(ϕ_symb, 1, 2, [2, 1]))
    display(U_1_LGT.dSdϕ_2d(ϕ_symb, 1)[1, 2, 1])
    display(U_1_LGT.dSdϕ_2d(ϕ_symb, 1, 1, [2, 1]))
    display(U_1_LGT.dSdϕ_2d(ϕ_symb, 1)[2, 2, 2])
    display(U_1_LGT.dSdϕ_2d(ϕ_symb, 1, 2, [2, 2]))
    #

end

function symbolic3_test()

    @variables y[1:3, 1:2, 1:1, 1:1]
    yy = zeros(Num, (3, 2, 1, 1))
    yy[1, 1] = y[1, 1, 1, 1]
    yy[1, 2] = y[1, 2, 1, 1]
    yy[2, 1] = y[2, 1, 1, 1]
    yy[2, 2] = y[2, 2, 1, 1]
    yy[3, 1] = y[3, 1, 1, 1]
    yy[3, 2] = y[3, 2, 1, 1]
    display(yy[1, 1, :, :])
    ϕ3_symb = reshape(yy, (3, 2, 1, 1))
    display(S_3d(ϕ3_symb, 1))
    display(S_v2_3d(ϕ3_symb, 1))

end


function benchmarkS2d()
    ϕs = randn((1000, 2, 100, 100))
    benchmark2d(U_1_LGT.S_2d, ϕs, true)
    t1 = benchmark2d(U_1_LGT.S_2d, ϕs, true)
    benchmark2d(U_1_LGT.S_v2_2d, ϕs, true)
    t2 = benchmark2d(U_1_LGT.S_v2_2d, ϕs, true)
    benchmark2d(U_1_LGT.S_2d, ϕs, true)
    t3 = benchmark2d(U_1_LGT.S_2d, ϕs, true)
    display(sum(abs.(t1 - t3)))
    display(sum(abs.(t2 - t3)))

end
function benchmarkdSdϕ2d()
    ϕs = randn((1000, 2, 100, 100))
    benchmark2d(U_1_LGT.dSdϕ_2d, ϕs, false)
    benchmark2d(U_1_LGT.dSdϕ_2d, ϕs, false)

end
function benchmarkS3d()
    ϕs = randn((100, 3, 20, 20, 20))
    benchmark3d(U_1_LGT.S_3d, ϕs, true)
    t1 = benchmark3d(U_1_LGT.S_3d, ϕs, true)
    benchmark3d(U_1_LGT.S_v2_3d, ϕs, true)
    t2 = benchmark3d(U_1_LGT.S_v2_3d, ϕs, true)
    benchmark3d(U_1_LGT.S_3d, ϕs, true)
    t3 = benchmark3d(U_1_LGT.S_3d, ϕs, true)
    display(sum(abs.(t1 - t3)))
    display(sum(abs.(t2 - t3)))
end
function benchmarkdSdϕ3d()
    ϕs = randn((100, 3, 20, 20, 20))
    benchmark3d(U_1_LGT.dSdϕ_3d, ϕs, false)
    benchmark3d(U_1_LGT.dSdϕ_3d, ϕs, false)

end


function benchmark2d(S_func, ϕs, store)

    function test_store(S_func, ϕs)
        N = size(ϕs)[1]
        return [S_func(ϕs[n, :, :, :], 1) for n in 1:N]

    end
    function test(S_func, ϕs)
        N = size(ϕs)[1]
        for n in 1:N
            S_func(ϕs[n, :, :, :], 1)
        end

    end
    if store
        @time test_store(S_func, ϕs)
    else
        @time test(S_func, ϕs)
    end
end

function benchmark3d(S_func, ϕs, store)

    function test_store(S_func, ϕs)
        N = size(ϕs)[1]
        return [S_func(ϕs[n, :, :, :, :], 1) for n in 1:N]
    end
    function test(S_func, ϕs)
        N = size(ϕs)[1]
        for n in 1:N
            S_func(ϕs[n, :, :, :, :], 1)
        end
    end

    if store
        @time test_store(S_func, ϕs)
    else
        @time test(S_func, ϕs)
    end
end



function ΔHΔτ(ϕ, S, dSdϕ, S_args, hmc_runs, initial_steps)

    for _ in 1:initial_steps
        ϕ = HMC.hmc_run(ϕ, S, dSdϕ, 10, 0.0001, S_args)[1]
    end

    Δτs = [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005]
    mean_ΔHs = []
    for Δτ in Δτs
        lf_steps = Int(round(1 / (100 * Δτ)))
        print("$Δτ\n")
        ΔHs = []
        for run in 1:hmc_runs
            ϕ_f, ΔH, accepted = HMC.hmc_run(ϕ, S, dSdϕ, lf_steps, Δτ, S_args)
            push!(ΔHs, ΔH)
        end
        push!(mean_ΔHs, mean(ΔHs))
    end
    plot = Plots.plot(Δτs, abs.(mean_ΔHs), xaxis=:log, yaxis=:log, markershape=:cross, linealpha=0, legend=:none)
    plot = Plots.xaxis!(plot, "log Δτ")
    plot = Plots.yaxis!(plot, "log mean ΔH")
    plot |> display
end

function random_gauge_transform_action(N, ϕ, S, S_args)
    # assumes first dimension of ϕ is direction: e1,e2,e3
    if N == 3
        χ = (rand(Float64, size(ϕ)[2:end]) * 2 .- 1) * pi
        Δχ1 = ShiftedArrays.circshift(χ, -[1, 0, 0]) - χ
        Δχ2 = ShiftedArrays.circshift(χ, -[0, 1, 0]) - χ
        Δχ3 = ShiftedArrays.circshift(χ, -[0, 0, 1]) - χ
        ϕ_n = copy(ϕ)
        ϕ_n[1, :, :, :] += Δχ1
        ϕ_n[2, :, :, :] += Δχ2
        ϕ_n[3, :, :, :] += Δχ3

        @info "Random Gauge Transformation" abs(S(ϕ, S_args...) - S(ϕ_n, S_args...))
    else
        display("Not implemented for N ≠ 3")
    end
end

function random_gauge_transform_wilson_loop(N, ϕ, loop, loop_args)
    # assumes first dimension of ϕ is direction: e1,e2,e3
    if N == 3
        χ = (rand(Float64, size(ϕ)[2:end]) * 2 .- 1) * pi
        Δχ1 = ShiftedArrays.circshift(χ, -[1, 0, 0]) - χ
        Δχ2 = ShiftedArrays.circshift(χ, -[0, 1, 0]) - χ
        Δχ3 = ShiftedArrays.circshift(χ, -[0, 0, 1]) - χ
        ϕ_n = copy(ϕ)
        ϕ_n[1, :, :, :] += Δχ1
        ϕ_n[2, :, :, :] += Δχ2
        ϕ_n[3, :, :, :] += Δχ3
        @info "Random Gauge Transformation" abs(loop(ϕ, loop_args...) - loop(ϕ_n, loop_args...))
    else
        display("Not implemented for N ≠ 3")
    end
end

function main()

    # test_p_dist(U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 3, 16, 100, 10, 0.05, (1,), 300)

    # test_reversibility(zeros((2, 100, 100)), U_1_LGT.dSdϕ_2d, 1000, 0.05, (1,))

    # test_derivative(U_1_LGT.dSdϕ_2d, U_1_LGT.dSdϕ_2d, 2,100)

    # dSdϕ(ϕ, β, μ, n) = U_1_LGT.dSdϕ_3d(ϕ, β)[μ, n...]
    # test_derivative_numerical(U_1_LGT.S_3d, dSdϕ, 3, 30, 0.0001)


    # test_expval(zeros((2, 10, 10)), U_1_LGT.S_2d, U_1_LGT.dSdϕ_2d, 100, 0.05, (1,), 10000, 300)
    # test_expval(zeros((3, 16, 16, 16)), U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 10, 0.1, (1,), 1000, 100)

    # benchmarkS2d()
    # benchmarkS3d()

    # benchmarkdSdϕ2d()
    # benchmarkdSdϕ3d()
    # symbolic2_test()

    # ΔHΔτ(zeros(3, 16, 16, 16), U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, (1,), 100, 100)

    # random_gauge_transform_action(3, rand(Float64, (3, 32, 32, 32)), U_1_LGT.S_3d, (1,))
    # random_gauge_transform_wilson_loop(3, rand(Float64, (3, 32, 32, 32)), site_average_loop3d, (1, 3, 8, 3))
    # random_gauge_transform_wilson_loop(3, rand(Float64, (3, 32, 32, 32)), rect_wilson_loop_3d, (4, 7, 1, 0.7, 10))

    # fa_test_expval3d(zeros((3, 16, 16, 16)), U_1_LGT.S_3d, U_1_LGT.dSdϕ_3d, 3, 1/3, (1,), 0.1, 2000, 200)

    function parisi_loop(ϕ, size)
        MeasurementFunctions.parisi_loop_range_3d(ϕ, 1, [size])[1]
    end
    function flux_tube(ϕ, τ)
        MeasurementFunctions.flux_tubes_x_3d(ϕ, 0.7, 10, [τ])[1]
    end
    # random_gauge_transform_wilson_loop(3, U_1_LGT.rand_lattice_3d(16, 16, 16), parisi_loop, ((10, 2),))
    random_gauge_transform_wilson_loop(3, U_1_LGT.rand_lattice_3d(16, 16, 16), flux_tube, [5])
end

main()
