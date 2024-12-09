include("action.jl")
include("hmc.jl")
include("data_handler.jl")

using Statistics
import Plots
using DelimitedFiles
using LinearAlgebra

import ShiftedArrays as sa
import SpecialFunctions as sf
#initial field
#update field
#steps_to_thermalise
#measurement interval
#number of measurements
#list of measurent functions
#update_field args
#measurement args

function optimal_lf_steps(update, update_args)
    trajectory_length = 1
    initial_steps = 300
    measurement_steps = 300
    rate_min = 0.65
    rate_max = 0.7

    accept_rate = 0
    while accept_rate <
          for _ in 1:initial_steps
        ϕ = update(ϕ, update_args...)[1]
    end
        for _ in 1:steps
            ϕ, ΔH, accepted = update(ϕ, update_args...)

            accept_count += accepted
        end
        accept_count / steps
    end
end

function mean_plaquette3d(ϕ, J)
    1 - Action.S_SAshift_3d(ϕ, J) / (length(ϕ) * J)
end

function mean_plaquette2d(ϕ, J)

    1 - Action.S_SAshift_2d(ϕ, J) / (length(ϕ) * J)

end

const unit_vectors_3d = [I[1:3, k] for k in 1:3]
function loop3d_old(ϕ, μ, ν, rμ, rν, n)
    if μ == ν
        return 0
    end

    #can use mod1 to wrap indices

    μp_indices = [n + unit_vectors_3d[μ] * r for r in 0:rμ-1]
    μn_indices = [n + unit_vectors_3d[ν] * (rν - 1) + unit_vectors_3d[μ] * r for r in (rμ-1):-1:0]
    νp_indices = [n + unit_vectors_3d[μ] * (rμ - 1) + unit_vectors_3d[ν] * r for r in 0:rν-1]
    νn_indices = [n + unit_vectors_3d[ν] * r for r in rν-1:-1:0]

    loop = sum(ϕ[μ, ind...] for ind in μp_indices) +
           sum(ϕ[ν, ind...] for ind in νp_indices) +
           sum(-ϕ[μ, ind...] for ind in μn_indices) +
           sum(-ϕ[ν, ind...] for ind in νn_indices)
    loop
end

function loop3d(ϕ, μ, ν, rμ, rν, n)
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
    mean(loop3d(ϕ, μ, ν, rμ, rν, [nx, ny, nz]) for nx in 1:size(ϕ)[2] for ny in 1:size(ϕ)[3] for nz in 1:size(ϕ)[4])
end

function Rτ_loop_3d(ϕ, R, τ)
    mean([site_average_loop3d(ϕ, 1, 3, R, τ), site_average_loop3d(ϕ, 2, 3, R, τ)])
end
# function site_dir_average_loop3d(ϕ, rμ, rν)
#     mean(site_average_loop3d(ϕ, μ, ν, rμ, rν) for (μ, ν) in [(1, 2), (1, 3), (2, 3)])
# end
# function range_loop3d(ϕ, rμ_range, rν_range)
#     collect((rμ, rν, site_dir_average_loop3d(ϕ, rμ, rν)) for rμ in rμ_range for rν in rν_range)
# end

function wilson_loop3d(ϕ, ϕ_ape_smeared, ϕ_t, spatial_path, τ)
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

function rect_wilson_loop_3d(ϕ, R, τ, J, α, N_smears)
    # α = 0.7
    ϕ_ape_smeared = N_ape_smearing_3d(exp.(complex.(0, ϕ[1:2, :, :, :])), α, N_smears)
    # ϕ_t = thermally_average_3d(ϕ, J)
    # ϕ_ape_smeared = exp.(im * ϕ[1:2, :, :, :])
    ϕ_t = exp.(im * ϕ[3, :, :, :])
    x_path = [(1, 1, i, 1) for i in 1:R]
    y_path = [(1, 2, 1, i) for i in 1:R]
    x_loops = wilson_loop3d(ϕ, ϕ_ape_smeared, ϕ_t, x_path, τ)
    y_loops = wilson_loop3d(ϕ, ϕ_ape_smeared, ϕ_t, y_path, τ)
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
    (R_12, θ_12) = cos_2_sum_as_cos(A_1, A_2, θ_1, θ_2)
    (R_34, θ_34) = cos_2_sum_as_cos(A_3, A_4, θ_3, θ_4)
    (R, θ) = cos_2_sum_as_cos(R_12, R_34, θ_12, θ_34)
    return (R, θ)
end
function thermally_average_3d(ϕ, β)
    # Thermally averages the temporal direction (assumed e3) of ϕ

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

    ϕ_t = (conj(X_l) .* sf.besseli.(1, β * d)) ./ (d .* sf.besseli.(0, β * d))
    # R_θ = cos_4_sum_as_cos.(1, 1, 1, 1, -θ_1, -θ_2, -θ_3, -θ_4)

    # ϕ_t = map(R_θ_i -> sf.besseli(1, β * R_θ_i[1]) * exp(complex(0, R_θ_i[2])) / sf.besseli(0, β * R_θ_i[1]), R_θ)

    return ϕ_t
end

function parisi_dir(ϕ, β, direction)
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

        ϕ_t = (conj(X_l) .* sf.besseli.(1, β * d)) ./ (d .* sf.besseli.(0, β * d))
        return ϕ_t
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

        ϕ_t = (conj(X_l) .* sf.besseli.(1, β * d)) ./ (d .* sf.besseli.(0, β * d))
        return ϕ_t
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

        ϕ_t = (conj(X_l) .* sf.besseli.(1, β * d)) ./ (d .* sf.besseli.(0, β * d))
        # R_θ = cos_4_sum_as_cos.(1, 1, 1, 1, -θ_1, -θ_2, -θ_3, -θ_4)

        # ϕ_t = map(R_θ_i -> sf.besseli(1, β * R_θ_i[1]) * exp(complex(0, R_θ_i[2])) / sf.besseli(0, β * R_θ_i[1]), R_θ)

        return ϕ_t
    else
        @warn "Direction must be 1, 2, 3." direction
        return 0
    end
end

function parisi(ϕ, β)
    U_t = zeros(ComplexF64, size(ϕ))
    U_t[1, :, :, :] = parisi_dir(ϕ, β, 1)
    U_t[2, :, :, :] = parisi_dir(ϕ, β, 2)
    U_t[3, :, :, :] = parisi_dir(ϕ, β, 3)
    return U_t
end

function parisi_3d_loops(ϕ, R, τ, β)
    (D, Lx, Ly, Lz) = size(ϕ)
    if D != 3
        @warn "Only supports 3D lattices" D
        return 0
    end
    if τ < 3 || R < 3
        @warn "Only supports loops of size R > 2, τ > 2." R τ
        return 0
    end

    U_t = parisi(ϕ, β)
    U = exp.(complex.(0, ϕ))
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

# function n_ape_smearing_3d(ϕ, α, n)
#     ϕ_smeared = copy(ϕ)
#     for _ in 1:n
#         ϕ_smeared = ape_smearing_3d(ϕ_smeared, α)
#     end
#     return ϕ_smeared
# end

function sum_lattice3d(ϕ, s1, s2, s3)
    #assumes s1, s2, s3 are multiples of lattice size
    size_ϕ = size(ϕ)
    Dx = Int(size_ϕ[2] / s1)
    Dy = Int(size_ϕ[3] / s2)
    Dz = Int(size_ϕ[4] / s3)

    ϕ_sum = zeros(Float64, (3, Dx, Dy, Dz))
    for i in 1:Dx
        for j in 1:Dy
            for k in 1:Dz
                ϕ_sum[:, i, j, k] = sum(ϕ[:, 2*i-1:2*i, 2*j-1:2*j, 2*k-1:2*k], dims=(2, 3, 4))
            end
        end
    end
    ϕ_sum
end

function benchmark_loops(n_runs)
    ϕ = randn((3, 16, 16, 16))
    tests = []
    for _ in n_runs
        μ = rand([1, 2, 3])
        ν = rand([1, 2, 3])
        rμ = rand(1:16)
        rν = rand(1:16)
        n = [rand(1:16), rand(1:16), rand(1:16)]
        n = [1, 1, 1]
        push!(tests, (μ, ν, rμ, rν, n))
    end
    # @time map(t -> loop3d(ϕ, t...), tests)
    # @time map(t -> loop3d(ϕ, t...), tests)
    @time map(t -> loop3d_v2(ϕ, t...), tests)
    @time map(t -> loop3d_v2(ϕ, t...), tests)
end

# function optimal_hmc_step(ϕ_init, S, dSdϕ, S_args)

# end


function acceptance_rate(ϕ, update, update_args, steps, initial_steps)
    accept_count = 0
    for _ in 1:initial_steps
        ϕ = update(ϕ, update_args...)[1]
    end
    for _ in 1:steps
        ϕ, ΔH, accepted = update(ϕ, update_args...)

        accept_count += accepted
    end
    display(accept_count / steps)
end


function measure(filename, mode, ϕ_init, initial_steps, measurement_interval, Nmeasurements, update, update_args, measurement_functions, measurement_args, measurement_info)
    # assumes the first arguement of update is ϕ
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
            # display(measurement)
            # writedlm(io, [measurement_info... Iterators.flatten(measurement...)...])
            writedlm(io, transpose([measurement_info..., measurement...]))
        end
    end
    accept_rate = accept_count / (Nmeasurements * measurement_interval)
    # write(file, "final_config", ϕ)
    DataHandler.save_field("final_fields/$filename", ϕ)


    print("acceptance rate: $accept_rate\n")
    # end
    #save field for later con
end



# dataβ = [1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0]
# dataβ_full = [0.2, 0.4, 0.6, 0.8, 1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]

# dataP = [0.475, 0.629, 0.656, 0.704, 0.748, 0.790, 0.806, 0.834, 0.854, 0.869, 0.881]


# for β in dataβ_full
#     measure("test_full.txt", "a", zeros(Float64, (3, 16, 16, 16)), 100, 10, 20, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# end

# run16 = [0.614, 0.732, 0.744, 0.770, 0.793, 0.816, 0.827, 0.846, 0.862, 0.876, 0.887]
# run16 = [0.475, 0.630, 0.654, 0.704, 0.749, 0.792, 0.807, 0.836, 0.854, 0.869, 0.881]
# run16_wider = [-0.601, -0.232, 0.128, 0.440, 0.895, 0.903, 0.909, 0.913, 0.919]
# dataβ_full = [0.2, 0.4, 0.6, 0.8, 1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]

# run16_full = [-0.601, -0.232, 0.128, 0.440, 0.614, 0.732, 0.744, 0.770, 0.793, 0.816, 0.827, 0.846, 0.862, 0.876, 0.887, 0.895, 0.903, 0.909, 0.913, 0.919]
# dataβ_wider = [0.2, 0.4, 0.6, 0.8, 3.2, 3.4, 3.6, 3.8, 4.0]


# half16 = [-0.057,0.253,0.307,0.416,0.493,0.581,0.617,0.670,0.711,0.738,0.762]
# p = Plots.plot(dataβ, dataP, label="Hamer")
# p = Plots.plot!(p, dataβ, run16, label="Mine")
# p = Plots.plot!(dataβ[1:4], dataβ[1:4] / 2)
# p = Plots.xlabel!(p, "J")
# p = Plots.ylabel!(p, "<P>")


# p |> display
# β = 1
# measure("autocorrelation_1_05_1000.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.93

# measure("autocorrelation_1_05_100.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 100, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.91

# measure("autocorrelation_1_1_100.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 100, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.1, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.77

# β = 0.5
# measure("autocorrelation_05_05_1000.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.977

# β = 3
# measure("autocorrelation_3_05_1000.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.737


# β = 1
# measure("autocorrelation_1_05_1000_1.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.946

# β = 1
# measure("autocorrelation_1_05_1000_2.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 1000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
#0.939

# β = 1
# measure("autocorrelation_1_05_10000.txt", "w", zeros(Float64, (3, 16, 16, 16)), 100, 1, 10000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# 0.9343


β = 2
Δτ = 0.08
lf_steps = Int(round(1 / Δτ))
# display(acceptance_rate(zeros(Float64, (3, 16, 16, 16)), HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), 200, 100))
# measure("plaquette_v_trajectory_cold.txt", "w", zeros(Float64, (3, 16, 16, 16)), 0, 1, 150, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
β = 1
Δτ = 0.05
lf_steps = Int(round(1 / Δτ))
# acceptance_rate((rand(Float64,(3,16,16,16))*2 .-1)*pi, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), 200, 100)
# measure("plaquette_v_trajectory_hot.txt", "w", (rand(Float64, (3, 16, 16, 16)) * 2 .- 1) * pi, 0, 1, 150, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))

# β = 1
# Δτ = 0.115
# lf_steps = Int(round(1 / Δτ))
# for nth_run in 1:10
#     measure("autocorrelation/autocorrelation_$(β)_10000_$(nth_run).txt", "w", zeros(Float64, (3, 16, 16, 16)), 250, 1, 10000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, 10, 0.05, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# end

β = 1
Δτ = 0.115
lf_steps = Int(round(1 / Δτ))
# measure("wilson_loop1.txt", "w", zeros(Float64, (3, 16, 16, 16)), 250, 10, 2, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), (range_loop3d,), ((1:8, 1:8,),), (β, 16))

β = 2
R = 5
Δτ = 0.08
lf_steps = Int(round(1 / Δτ))
# for _ in 1:10
#     measure("wilson_loop_run_B2_R5_ot_1000_10_real.txt", "a", zeros(Float64, (3, 16, 16, 16)), 250, 10, 1000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), (rect_wilson_loop_3d, rect_wilson_loop_3d, rect_wilson_loop_3d, rect_wilson_loop_3d, rect_wilson_loop_3d, rect_wilson_loop_3d, rect_wilson_loop_3d, rect_wilson_loop_3d), ((R, 1, β, 0.7, 10), (R, 2, β, 0.7, 10), (R, 3, β, 0.7, 10), (R, 4, β, 0.7, 10), (R, 5, β, 0.7, 10), (R, 6, β, 0.7, 10), (R, 7, β, 0.7, 10), (R, 8, β, 0.7, 10)), (β, 16))
# end

β = 2
Δτ = 0.08
lf_steps = Int(round(1 / Δτ))
measurement_functions = (rect_wilson_loop_3d for R in 1:10 for τ in 1:10)
measurement_function_arguments = ((R, τ, β, 0.7, 10) for R in 1:10 for τ in 1:10)
# for _ in 1:10
#     measure("wilson_loop_run_long_thermal_range_ot_110_1500_10_real.txt", "a", zeros(Float64, (3, 16, 16, 16)), 1000, 10, 1500, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), measurement_functions, measurement_function_arguments, (β, 16))
# end

β = 2
Δτ = 0.05
lf_steps = Int(round(1 / Δτ))
measurement_functions = (parisi_3d_loops for _ in 1:36)
measurement_function_arguments = ((R, τ, β) for R in 3:10 for τ in R:10)
for _ in 1:10
    measure("wl_check_p_2_1000_1000.txt", "a", zeros(Float64, (3, 32, 32, 32)), 1000, 10, 1000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), measurement_functions, measurement_function_arguments, (β, 32))
end

# β = 2
# Δτ = 0.05
# lf_steps = Int(round(1 / Δτ))
# measurement_functions = (rect_wilson_loop_3d for R in 2:10 for τ in 2:10)
# measurement_function_arguments = ((R, τ, β, 0.7, 10) for R in 2:10 for τ in 2:10)
# for _ in 1:10
#     measure("wilson_loop32_run_range_no_1000_10_real.txt", "a", zeros(Float64, (3, 16, 16, 16)), 250, 10, 1000, HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), measurement_functions, measurement_function_arguments, (β, 16))
# end

# acceptance_rate(zeros(Float64,(3,32,32,32)), HMC.hmc_run, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, lf_steps, Δτ, (β,)), 500, 250)

# β = 3
Δτ = 0.125
lf_steps = Int(round(1 / Δτ))
inverse_FK = HMC.inv_FK_3d(16, 16, 16, 0.1)
# acceptance_rate(zeros(Float64, (3, 16, 16, 16)), HMC.fa_hmc_run_3d, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, inverse_FK, lf_steps, Δτ, (β,)), 500, 250)

dataβ = [1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0]
# dataβ_full = [0.2, 0.4, 0.6, 0.8, 1.0, 1.35, 1.41, 1.55, 1.70, 1.90, 2.0, 2.25, 2.5, 2.75, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]

# dataP = [0.475, 0.629, 0.656, 0.704, 0.748, 0.790, 0.806, 0.834, 0.854, 0.869, 0.881]


# for β in dataβ
#     measure("fa_hamer_p_2.txt", "a", zeros(Float64, (3, 16, 16, 16)), 250, 10, 100, HMC.fa_hmc_run_3d, (Action.S_SAshift_3d, Action.dSdϕ_SAshift_3d, inverse_FK, lf_steps, Δτ, (β,)), (mean_plaquette3d,), ((β,),), (β, 16))
# end
