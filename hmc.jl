module HMC

import FFTW
import ShiftedArrays as sa

export leapfrog
export hmc_run

export fa_hmc_run_3d
export fa_leapfrog



function hmc_run(ϕ, S, dSdϕ, lf_steps, Δτ, S_args)
    # Do I need to change the standard deviation of the momenta p? No there is a scale between computer time and these momenta

    p_0 = randn(Float64, size(ϕ))

    H_init = S(ϕ, S_args...) + 1 / 2 * sum(p_0 .^ 2)

    ϕ_new, p_new = leapfrog(ϕ, p_0, dSdϕ, lf_steps, Δτ, S_args)

    H_final = S(ϕ_new, S_args...) + 1 / 2 * sum(p_new .^ 2)

    ΔH = H_final - H_init
    r = rand(Float64)

    if r < exp(-ΔH)
        return ϕ_new, ΔH, true
    else
        return ϕ, ΔH, false
    end
end
function hmc_run_noMS(ϕ, S, dSdϕ, lf_steps, Δτ, S_args)
    # Useful when initialising from random
    p_0 = randn(Float64, size(ϕ))

    ϕ_new, p_new = leapfrog(ϕ, p_0, dSdϕ, lf_steps, Δτ, S_args)

    return ϕ_new
end

function leapfrog(ϕ, p, dSdϕ, lf_steps, Δτ, S_args)
    ϕ_half = ϕ + p * Δτ / 2
    for _ in 1:lf_steps-1
        p = p - dSdϕ(ϕ_half, S_args...) * Δτ
        ϕ_half = ϕ_half + p * Δτ
    end

    p_final = p - dSdϕ(ϕ_half, S_args...) * Δτ
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

function inv_FK_3d_mass(Lx, Ly, Lz, M)
    # K = M^2 - ∂^2
    inverse_FK_along_x = reshape([1 / (4 * sin(π * k1 / Lx)^2 + 4 * sin(π * k2 / Ly)^2 + 4 * sin(π * k3 / Lz)^2 + M^2) for k3 in 0:Lz-1 for k2 in 0:Ly-1 for k1 in 0:Lx-1], (Lx, Ly, Lz))
    inverse_FK = zeros((3, Lx, Ly, Lz))
    inverse_FK[1, :, :, :] = inverse_FK_along_x
    inverse_FK[2, :, :, :] = inverse_FK_along_x
    inverse_FK[3, :, :, :] = inverse_FK_along_x
    return inverse_FK
end

function inv_FK_3d_kappa(Lx, Ly, Lz, κ)
    # K = (1-κ)-κ∂^2
    inverse_FK_along_x = reshape([1 / (κ * (4 * sin(π * k1 / Lx)^2 + 4 * sin(π * k2 / Ly)^2 + 4 * sin(π * k3 / Lz)^2 - 1) + 1) for k3 in 0:Lz-1 for k2 in 0:Ly-1 for k1 in 0:Lx-1], (Lx, Ly, Lz))
    inverse_FK = zeros((3, Lx, Ly, Lz))
    inverse_FK[1, :, :, :] = inverse_FK_along_x
    inverse_FK[2, :, :, :] = inverse_FK_along_x
    inverse_FK[3, :, :, :] = inverse_FK_along_x
    return inverse_FK
end



function sample_Fp_3d(inverse_FK)
    inverse_FK_size = size(inverse_FK)
    (D, Lx, Ly, Lz) = inverse_FK_size
    Π_1 = map(Ak -> randn() * sqrt((Lx * Ly * Lz) / Ak), inverse_FK[1, :, :, :])
    Π_2 = map(Ak -> randn() * sqrt((Lx * Ly * Lz) / Ak), inverse_FK[2, :, :, :])
    Π_3 = map(Ak -> randn() * sqrt((Lx * Ly * Lz) / Ak), inverse_FK[3, :, :, :])

    Fp = zeros(ComplexF64, inverse_FK_size)
    Fp[1, :, :, :] = fourier_space_momenta_3d(Π_1, Lx, Ly, Lz)
    Fp[2, :, :, :] = fourier_space_momenta_3d(Π_2, Lx, Ly, Lz)
    Fp[3, :, :, :] = fourier_space_momenta_3d(Π_3, Lx, Ly, Lz)
    return Fp
end

function Fke(inverse_FK, Fp, Nsites)
    real(sum((conj(Fp) .* inverse_FK .* Fp)) / (2 * Nsites))
end

function fa_hmc_run_3d(ϕ, S, dSdϕ, lf_steps, Δτ, inverse_FK, S_args)

    # Momenta generation and fft_dims do not generalise
    # The are the dimentions to apply transforms to
    # Can be changed to a parameter or generalised to (2,...,D+1)
    (D, Lx, Ly, Lz) = size(inverse_FK)
    Nsites = Lx * Ly * Lz
    fft_dims = (2, 3, 4)
    Fp_0 = sample_Fp_3d(inverse_FK)
    p_0 = real(FFTW.ifft(Fp_0, fft_dims))
    H_init = S(ϕ, S_args...) + Fke(inverse_FK, Fp_0, Nsites)

    ϕ_new, p_new, Fp_new = fa_leapfrog(ϕ, p_0, dSdϕ, inverse_FK, fft_dims, lf_steps, Δτ, S_args)

    H_final = S(ϕ_new, S_args...) + Fke(inverse_FK, Fp_new, Nsites)

    ΔH = H_final - H_init

    r = rand(Float64)

    if r < exp(-ΔH)
        return ϕ_new, ΔH, true
    else
        return ϕ, ΔH, false
    end
end

function fourier_space_momenta_3d(Π, Lx, Ly, Lz)
    # Not optimal
    # Assumes Lx,Ly,Lz are all even

    # Can do to remove Lx,Ly,Lz from arguments
    # (Lx, Ly, Lz) = size(Π)

    i_real_indices = [0, Int(Lx / 2)]
    j_real_indices = [0, Int(Ly / 2)]
    k_real_indices = [0, Int(Lx / 2)]
    function conj_index(i, j, k)
        (mod1(Lx - i + 1, Lx), mod1(Ly - j + 1, Ly), mod1(Lz - k + 1, Lz))
    end

    Fp = zeros(ComplexF64, (Lx, Ly, Lz))
    for i in 0:Lx-1
        for j in 0:Ly-1
            for k in 0:Int(Lz / 2)

                if i in i_real_indices && j in j_real_indices && k in k_real_indices
                    Fp[i+1, j+1, k+1] = Π[i+1, j+1, k+1]
                else
                    Fp[i+1, j+1, k+1] = complex(Π[i+1, j+1, k+1], Π[conj_index(i, j, k)...]) / sqrt(2)
                    Fp[conj_index(i, j, k)...] = conj(Fp[i+1, j+1, k+1])
                end
            end
        end
    end

    return Fp
end

function fa_leapfrog(ϕ, p, dSdϕ, inverse_FK, fft_dims, lf_steps, Δτ, S_args)

    ϕ_half = ϕ + real(FFTW.ifft(inverse_FK .* FFTW.fft(p, fft_dims), fft_dims)) * Δτ / 2

    for _ in 1:lf_steps-1
        p = p - dSdϕ(ϕ_half, S_args...) * Δτ
        ϕ_half = ϕ_half + real(FFTW.ifft(inverse_FK .* FFTW.fft(p, fft_dims), fft_dims)) * Δτ
    end
    p_final = p - dSdϕ(ϕ_half, S_args...) * Δτ
    Fp_final = FFTW.fft(p_final, fft_dims)
    ϕ_final = ϕ_half + real(FFTW.ifft(inverse_FK .* Fp_final, fft_dims)) * Δτ / 2
    return ϕ_final, p_final, Fp_final
end

end
