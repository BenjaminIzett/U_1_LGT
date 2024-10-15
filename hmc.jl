module HMC

export leapfrog
export hmc_run

function hmc_run(ϕ, S, dSdϕ, lf_steps, Δτ, S_args)
    # Do I need to change the standard deviation of the momenta p? No there is a scale between computer time and these momenta

    p_0 = randn(Float64, size(ϕ))

    H_init = S(ϕ, S_args...) + 1 / 2 * sum(p_0 .^ 2)

    ϕ_new, p_new = leapfrog(copy(ϕ), p_0, dSdϕ, lf_steps, Δτ, S_args)

    H_final = S(ϕ_new, S_args...) + 1 / 2 * sum(p_new .^ 2)

    ΔH = H_final - H_init
    r = rand(Float64)

    if r < exp(-ΔH)
        return ϕ_new, ΔH, true
    else
        return ϕ, ΔH, false
    end
end

function leapfrog(ϕ_0, p_0, dSdϕ, lf_steps, Δτ, S_args)
    ϕ_half = ϕ_0 + p_0 * Δτ / 2
    for _ in 1:lf_steps-1
        p_0 = p_0 - dSdϕ(ϕ_half, S_args...) * Δτ
        ϕ_half = ϕ_half + p_0 * Δτ
    end

    p_final = p_0 - dSdϕ(ϕ_half, S_args...) * Δτ
    ϕ_final = ϕ_half + p_final * Δτ / 2
    return ϕ_final, p_final
end


function leapfrog_pfirst(ϕ_0, p_0, dSdϕ, N, Δτ, S_args)
    for _ in 1:N
        p_0 = p_0 - Δτ / 2 .* dSdϕ(ϕ_0, S_args...)
        ϕ_0 = ϕ_0 + Δτ .* p_0
        p_0 = p_0 - Δτ / 2 .* dSdϕ(ϕ_0, S_args...)
    end
    return ϕ_0, p_0
end

end
