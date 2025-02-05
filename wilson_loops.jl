module WilsonLoops
using Statistics
using LinearAlgebra

import ShiftedArrays as sa
import SpecialFunctions as sf

export *

const unit_vectors_3d = [I[1:3, k] for k in 1:3]

function loop_3d(ϕ, μ, ν, rμ, rν, n)
    if μ == ν
        return 0
    end
    size_μ = size(ϕ)[μ+1]
    size_ν = size(ϕ)[ν+1]

    α = 6 - μ - ν # the index perp to loop

    μp_indices = mod1.(n[μ] .+ (0:rμ), size_μ)
    νp_indices = mod1.(n[ν] .+ (0:rν), size_ν)


    #can reverse order in last two terms but since we are just adding it doesnt matter

    loop = sum(ϕ[μ, (unit_vectors_3d[α] * n[α] + unit_vectors_3d[μ] * ind + unit_vectors_3d[ν] * νp_indices[1])...] for ind in μp_indices[1:end-1]) +
           sum(ϕ[ν, (unit_vectors_3d[α] * n[α] + unit_vectors_3d[μ] * μp_indices[end] + unit_vectors_3d[ν] * ind)...] for ind in νp_indices[1:end-1]) +
           sum(-ϕ[μ, (unit_vectors_3d[α] * n[α] + unit_vectors_3d[μ] * ind + unit_vectors_3d[ν] * νp_indices[end])...] for ind in μp_indices[1:end-1]) +
           sum(-ϕ[ν, (unit_vectors_3d[α] * n[α] + unit_vectors_3d[μ] * μp_indices[1] + unit_vectors_3d[ν] * ind)...] for ind in νp_indices[1:end-1])
    exp(complex(0, loop))
end

function site_average_loop3d(ϕ, μ, ν, rμ, rν)
    mean(loop_3d(ϕ, μ, ν, rμ, rν, [nx, ny, nz]) for nx in 1:size(ϕ)[2] for ny in 1:size(ϕ)[3] for nz in 1:size(ϕ)[4])
end

function Rτ_loop_3d(ϕ, R, τ)
    mean([site_average_loop3d(ϕ, 1, 3, R, τ), site_average_loop3d(ϕ, 2, 3, R, τ)])
end

function general_hamer_loop_3d(ϕ, ϕ_ape_smeared, ϕ_t, spatial_path, τ)
    # Assumes that the temporal direction is in the e_3 direction ie [0,0,1] basis vector
    # τ should be <= half the temporal lattice width
    # a path in the spatial plane is defined as a list of tuples (sign, direction, x, y)

    (D, Lx, Ly, Lz) = size(ϕ)
    if D != 3
        @warn "Only supports 3D lattices" D
        return 0
    end

    plane_basis = [(1, 0), (0, 1)]
    #could move this calculation elsewhere and make it a parameter
    endpoint =
        if spatial_path[end][1] == 1
            spatial_path[end][3:4] .+ spatial_path[end][1] .* plane_basis[spatial_path[end][2]]
        else
            spatial_path[end][3:4]
        end

    loops = zeros(ComplexF64, Lz)
    for k1 in 1:Lz
        k2 = mod1(k1 + τ, Lz)

        #top and bottom planes where the loops move along spatially

        p1 = ϕ_ape_smeared[:, :, :, k1]

        p2 = ϕ_ape_smeared[:, :, :, k2]
        # temporal_inds =
        #     if mod1(k1 + 1, Lz) < mod1(k2 - 2, Lz)
        #         mod1(k1 + 1, Lz):mod1(k2 - 2, Lz)
        #     else
        #         mod1(k1 + 1, Lz):-1:mod1(k2 - 2, Lz)
        #     end

        temporal_inds = mod1.(k1+1:k1+τ-2, Lz)

        # temporal_inds = mod1(k1 + 1, Lz):mod1(k2 - 2, Lz)

        # a 2x2x1 array such that temporal_lines[x,y] gives the sum of links along the z from z=i1 to z=i2-1
        # only works for τ>2
        temporal_lines =
            if τ == 0
                ones(ComplexF64, size(ϕ_t[:, :, 1]))
            elseif τ == 1
                exp.(complex.(0, ϕ[3, :, :, k1]))
            elseif τ == 2
                exp.(complex.(0, ϕ[3, :, :, k1])) .* exp.(complex.(0, ϕ[3, :, :, mod1(k2 - 1, Lz)]))
            else
                exp.(complex.(0, ϕ[3, :, :, k1])) .* prod(ϕ_t[:, :, temporal_inds], dims=3) .* exp.(complex.(0, ϕ[3, :, :, mod1(k2 - 1, Lz)]))

            end

        # temporal_lines = exp.(complex.(0, ϕ[3, :, :, k1])) .* prod(ϕ_t[:, :, temporal_inds], dims=3) .* exp.(complex.(0, ϕ[3, :, :, mod1(k2 - 1, Lz)]))

        # p1_along_path = prod(r[1] * sa.circshift(p1[r[2], :, :], ((1, 1) .- r[3:4])) for r in spatial_path)
        p1_paths = [r[1] == 1 ? sa.circshift(p1[r[2], :, :], ((1, 1) .- r[3:4])) : conj(sa.circshift(p1[r[2], :, :], ((1, 1) .- r[3:4]))) for r in spatial_path]
        p1_along_path = reduce((F, t) -> F .* t, p1_paths)

        # subtract because we are going backwards along the path
        # p2_along_path = prod(r[1] * sa.circshift(p2[r[2], :, :], ((1, 1) .- r[3:4])) for r in spatial_path)
        p2_paths = [r[1] == -1 ? sa.circshift(p2[r[2], :, :], ((1, 1) .- r[3:4])) : conj(sa.circshift(p2[r[2], :, :], ((1, 1) .- r[3:4]))) for r in spatial_path]
        p2_along_reversed_path = reduce((F, t) -> F .* t, p2_paths)

        loops_k1 = p1_along_path .* circshift(temporal_lines, ((1, 1) .- endpoint)) .* p2_along_reversed_path .* conj(temporal_lines)

        loops[k1] = mean(loops_k1)
    end
    mean(loops)
end

function hamer_loop_3d(ϕ, ϕ_ape_smeared, U_t, R, τ, β, α, N_smears)
    # α = 0.7
    # ϕ_ape_smeared = N_ape_smearing_3d(exp.(complex.(0, ϕ[1:2, :, :, :])), α, N_smears)
    # U_t = thermally_average_3d(ϕ, β)
    # ϕ_ape_smeared = exp.(im * ϕ[1:2, :, :, :])
    # ϕ_t = exp.(im * ϕ[3, :, :, :])
    x_path = [(1, 1, i, 1) for i in 1:R]
    y_path = [(1, 2, 1, i) for i in 1:R]
    x_loops = general_hamer_loop_3d(ϕ, ϕ_ape_smeared, U_t, x_path, τ)
    y_loops = general_hamer_loop_3d(ϕ, ϕ_ape_smeared, U_t, y_path, τ)
    return real(mean([x_loops, y_loops]))
end


function ape_smearing_3d(exp_iϕ, α)
    # maybe change to specify a plane instead of smearing all of them
    # Assumes the temporal direction is in e3 = [0, 0, 1]
    ϕ_smeared = zeros(ComplexF64, (2, size(exp_iϕ)[2:end]...))

    ϕ_smeared[1, :, :, :] = α * exp_iϕ[1, :, :, :] + ((1 - α) / 2) * (
        exp_iϕ[2, :, :, :] .*
        sa.circshift(exp_iϕ[1, :, :, :], -[0, 1, 0]) .*
        conj(sa.circshift(exp_iϕ[2, :, :, :], -[1, 0, 0])) +
        conj(sa.circshift(exp_iϕ[2, :, :, :], -[0, -1, 0])) .*
        sa.circshift(exp_iϕ[1, :, :, :], -[0, -1, 0]) .*
        sa.circshift(exp_iϕ[2, :, :, :], -[1, -1, 0])
    )
    ϕ_smeared[2, :, :, :] = α * exp_iϕ[2, :, :, :] + ((1 - α) / 2) * (
        exp_iϕ[1, :, :, :] .*
        sa.circshift(exp_iϕ[2, :, :, :], -[1, 0, 0]) .*
        conj(sa.circshift(exp_iϕ[1, :, :, :], -[0, 1, 0])) +
        conj(sa.circshift(exp_iϕ[1, :, :, :], -[-1, 0, 0])) .*
        sa.circshift(exp_iϕ[2, :, :, :], -[-1, 0, 0]) .*
        sa.circshift(exp_iϕ[1, :, :, :], -[-1, 1, 0])
    )

    # Project onto group element
    # Maybe check that elements are non zero
    ϕ_smeared = ϕ_smeared ./ abs.(ϕ_smeared)
    return ϕ_smeared
end

function N_ape_smearing_3d(exp_iϕ, α, N)
    for _ in 1:N
        exp_iϕ = ape_smearing_3d(exp_iϕ, α)
    end
    return exp_iϕ
end

function cos_2_sum_as_cos(A_1, A_2, θ_1, θ_2)
    # Expresses A_1*cos(x-θ_1)+A_2*cos(x-θ_2) as R*cos(x-θ)
    # Returns (R,θ)
    R = sqrt(A_1^2 + A_2^2 + 2 * A_1 * A_2 * cos(θ_1 - θ_2))
    θ = atan(A_1 * sin(θ_1) + A_2 * sin(θ_2), A_1 * cos(θ_1) + A_2 * cos(θ_2))
    return (R, θ)
end
function cos_4_sum_as_cos(A_1, A_2, A_3, A_4, θ_1, θ_2, θ_3, θ_4)
    # Expresses A_1*cos(x-θ_1)+A_2*cos(x-θ_2)+A_3*cos(x-θ_3)+A_4*cos(x-θ_4) as R*cos(x-θ)
    # Returns (R,θ)
    (R_12, θ_12) = cos_2_sum_as_cos(A_1, A_2, θ_1, θ_2)
    (R_34, θ_34) = cos_2_sum_as_cos(A_3, A_4, θ_3, θ_4)
    (R, θ) = cos_2_sum_as_cos(R_12, R_34, θ_12, θ_34)
    return (R, θ)
end

function parisi_dir_3d(ϕ, β, direction)
    # Thermally averages in given direction
    if direction == 1
        θ_1 = sa.circshift(ϕ[2, :, :, :], -[1, 0, 0]) -
              sa.circshift(ϕ[1, :, :, :], -[0, 1, 0]) -
              ϕ[2, :, :, :]
        θ_2 = sa.circshift(ϕ[3, :, :, :], -[1, 0, 0]) -
              sa.circshift(ϕ[1, :, :, :], -[0, 0, 1]) -
              ϕ[3, :, :, :]
        θ_3 = sa.circshift(ϕ[2, :, :, :], -[0, -1, 0]) -
              sa.circshift(ϕ[1, :, :, :], -[0, -1, 0]) -
              sa.circshift(ϕ[2, :, :, :], -[1, -1, 0])
        θ_4 = sa.circshift(ϕ[3, :, :, :], -[0, 0, -1]) -
              sa.circshift(ϕ[1, :, :, :], -[0, 0, -1]) -
              sa.circshift(ϕ[3, :, :, :], -[1, 0, -1])

        X_l = exp.(complex.(0, θ_1)) + exp.(complex.(0, θ_2)) + exp.(complex.(0, θ_3)) + exp.(complex.(0, θ_4))
        d = abs.(X_l)

        U_t = (conj(X_l) .* sf.besseli.(1, β * d)) ./ (d .* sf.besseli.(0, β * d))
        return U_t
    elseif direction == 2
        θ_1 = sa.circshift(ϕ[1, :, :, :], -[0, 1, 0]) -
              sa.circshift(ϕ[2, :, :, :], -[1, 0, 0]) -
              ϕ[1, :, :, :]
        θ_2 = sa.circshift(ϕ[3, :, :, :], -[0, 1, 0]) -
              sa.circshift(ϕ[2, :, :, :], -[0, 0, 1]) -
              ϕ[3, :, :, :]
        θ_3 = sa.circshift(ϕ[1, :, :, :], -[-1, 0, 0]) -
              sa.circshift(ϕ[2, :, :, :], -[-1, 0, 0]) -
              sa.circshift(ϕ[1, :, :, :], -[-1, 1, 0])
        θ_4 = sa.circshift(ϕ[3, :, :, :], -[0, 0, -1]) -
              sa.circshift(ϕ[2, :, :, :], -[0, 0, -1]) -
              sa.circshift(ϕ[3, :, :, :], -[0, 1, -1])

        X_l = exp.(complex.(0, θ_1)) + exp.(complex.(0, θ_2)) + exp.(complex.(0, θ_3)) + exp.(complex.(0, θ_4))
        d = abs.(X_l)

        U_t = (conj(X_l) .* sf.besseli.(1, β * d)) ./ (d .* sf.besseli.(0, β * d))
        return U_t
    elseif direction == 3
        θ_1 = sa.circshift(ϕ[1, :, :, :], -[0, 0, 1]) -
              sa.circshift(ϕ[3, :, :, :], -[1, 0, 0]) -
              ϕ[1, :, :, :]
        θ_2 = sa.circshift(ϕ[2, :, :, :], -[0, 0, 1]) -
              sa.circshift(ϕ[3, :, :, :], -[0, 1, 0]) -
              ϕ[2, :, :, :]
        θ_3 = sa.circshift(ϕ[1, :, :, :], -[-1, 0, 0]) -
              sa.circshift(ϕ[3, :, :, :], -[-1, 0, 0]) -
              sa.circshift(ϕ[1, :, :, :], -[-1, 0, 1])
        θ_4 = sa.circshift(ϕ[2, :, :, :], -[0, -1, 0]) -
              sa.circshift(ϕ[3, :, :, :], -[0, -1, 0]) -
              sa.circshift(ϕ[2, :, :, :], -[0, -1, 1])

        X_l = exp.(complex.(0, θ_1)) + exp.(complex.(0, θ_2)) + exp.(complex.(0, θ_3)) + exp.(complex.(0, θ_4))
        d = abs.(X_l)

        U_t = (conj(X_l) .* sf.besseli.(1, β * d)) ./ (d .* sf.besseli.(0, β * d))
        # R_θ = cos_4_sum_as_cos.(1, 1, 1, 1, -θ_1, -θ_2, -θ_3, -θ_4)

        # ϕ_t = map(R_θ_i -> sf.besseli(1, β * R_θ_i[1]) * exp(complex(0, R_θ_i[2])) / sf.besseli(0, β * R_θ_i[1]), R_θ)

        return U_t
    else
        @warn "Direction must be 1, 2, 3." direction
        return 0
    end
end

function parisi_3d(ϕ, β)
    U_t = zeros(ComplexF64, size(ϕ))
    U_t[1, :, :, :] = parisi_dir_3d(ϕ, β, 1)
    U_t[2, :, :, :] = parisi_dir_3d(ϕ, β, 2)
    U_t[3, :, :, :] = parisi_dir_3d(ϕ, β, 3)
    return U_t
end

function parisi_3d_loops_old(U, U_t, R, τ, β)
    (D, Lx, Ly, Lz) = size(U)
    if D != 3
        @warn "Only supports 3D lattices" D
        return 0
    end
    if τ < 3 || R < 3
        @warn "Only supports loops of size R > 2, τ > 2." R τ
        return 0
    end

    # U_t = parisi_3d(ϕ, β)
    # U = exp.(complex.(0, ϕ))
    # U_t = copy(U)
    Ux = U[1, :, :, :]
    Uy = U[2, :, :, :]
    Uz = U[3, :, :, :]

    Ux_t = U_t[1, :, :, :]
    Uy_t = U_t[2, :, :, :]
    Uz_t = U_t[3, :, :, :]
    # return Ux, Uy, Uz, Ux_t, Uy_t, Uz_t
    x_lines = Ux .* reduce((F, t) -> F .* t, [sa.circshift(Ux_t, -[r, 0, 0]) for r in 1:R-2]) .* sa.circshift(Ux, -[R - 1, 0, 0])
    y_lines = Uy .* reduce((F, t) -> F .* t, [sa.circshift(Uy_t, -[0, r, 0]) for r in 1:R-2]) .* sa.circshift(Uy, -[0, R - 1, 0])
    τ_lines = Uz .* reduce((F, t) -> F .* t, [sa.circshift(Uz_t, -[0, 0, T]) for T in 1:τ-2]) .* sa.circshift(Uz, -[0, 0, τ - 1])

    x_loops = x_lines .* sa.circshift(τ_lines, -[R, 0, 0]) .* conj.(sa.circshift(x_lines, -[0, 0, τ])) .* conj.(τ_lines)

    y_loops = y_lines .* sa.circshift(τ_lines, -[0, R, 0]) .* conj.(sa.circshift(y_lines, -[0, 0, τ])) .* conj.(τ_lines)

    loops = mean(x_loops .+ y_loops) / 2
    return real(loops)
end

function parisi_3d_loops(U, U_t, R, τ, β)
    (D, Lx, Ly, Lz) = size(U)
    if D != 3
        @warn "Only supports 3D lattices" D
        return 0
    end
    if τ < 2 || R < 2
        @warn "Only supports loops of size R > 1, τ > 1." R τ
        return 0
    end

    # U_t = parisi_3d(ϕ, β)
    # U = exp.(complex.(0, ϕ))
    # U_t = copy(U)
    Ux = U[1, :, :, :]
    Uy = U[2, :, :, :]
    Uz = U[3, :, :, :]

    Ux_t = U_t[1, :, :, :]
    Uy_t = U_t[2, :, :, :]
    Uz_t = U_t[3, :, :, :]
    # return Ux, Uy, Uz, Ux_t, Uy_t, Uz_t
    x_lines = ComplexF64(0)
    y_lines = ComplexF64(0)
    if R == 2
        x_lines = Ux .* sa.circshift(Ux, -[R - 1, 0, 0])
        y_lines = Uy .* sa.circshift(Uy, -[0, R - 1, 0])
    else
        x_lines = Ux .* reduce((F, t) -> F .* t, [sa.circshift(Ux_t, -[r, 0, 0]) for r in 1:R-2]) .* sa.circshift(Ux, -[R - 1, 0, 0])
        y_lines = Uy .* reduce((F, t) -> F .* t, [sa.circshift(Uy_t, -[0, r, 0]) for r in 1:R-2]) .* sa.circshift(Uy, -[0, R - 1, 0])
    end
    τ_lines = reduce((F, t) -> F .* t, [sa.circshift(Uz_t, -[0, 0, T]) for T in 0:τ-1])

    x_loops = x_lines .* sa.circshift(τ_lines, -[R, 0, 0]) .* conj.(sa.circshift(x_lines, -[0, 0, τ])) .* conj.(τ_lines)

    y_loops = y_lines .* sa.circshift(τ_lines, -[0, R, 0]) .* conj.(sa.circshift(y_lines, -[0, 0, τ])) .* conj.(τ_lines)

    loops = (mean(x_loops) + mean(y_loops)) / 2
    return real(loops)
end

function flux_tubes_x(ϕ, α, N, τs)
    # Assumes the x direction is 1.

    # summed_ϕ = sum(ϕ[1, :, :, :], dims=1)
    # [mean(cos.(summed_ϕ .- sa.circshift(summed_ϕ, -[0, 0, τ]))) for τ in τs]

    U_smeared = N_ape_smearing_3d(exp.(complex.(0, ϕ)), α, N)
    walls = prod(U_smeared[1, :, :, :], dims=1)
    [sum(real(walls .* sa.circshift(conj(walls), -[0, 0, τ]))) for τ in τs]
end
function flux_tubes_x_test(ϕ, α, N)
    # Assumes the x direction is 1.

    # summed_ϕ = sum(ϕ[1, :, :, :], dims=1)
    # [mean(cos.(summed_ϕ .- sa.circshift(summed_ϕ, -[0, 0, τ]))) for τ in τs]

    U_smeared = N_ape_smearing_3d(exp.(complex.(0, ϕ)), α, N)
    walls = prod(U_smeared[1, :, :, :], dims=1)
    sum(real(walls))
end


end
