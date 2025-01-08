module U_1_LGT
using LinearAlgebra
import ShiftedArrays as sa

export *

const unit_vectors_2d = [I[1:2, k] for k in 1:2]
const unit_vectors_3d = [I[1:3, k] for k in 1:3]

#2D
function zero_lattice_2d(Lx, Ly)
    return zeros(Float64, (2, Lx, Ly))
end

function rand_lattice_2d(Lx, Ly)
    return π * (2 * rand(Float64, (2, Lx, Ly)) .- 1)
end


function curl_shifted_2d(ϕ, μ, ν, k)
    if μ == ν
        return 0
    else

        return sa.circshift(ϕ[ν, :, :], -(k + unit_vectors_2d[μ])) - sa.circshift(ϕ[ν, :, :], -k) - sa.circshift(ϕ[μ, :, :], -(k + unit_vectors_2d[ν])) + sa.circshift(ϕ[μ, :, :], -k)
    end
end

function S_2d(ϕ, β)
    β * sum(1 .- cos.(curl_shifted_2d(ϕ, 1, 2, [0, 0])))
end

function dSdϕ_2d(ϕ, β)
    dSdϕ = similar(ϕ)
    dSdϕ[1, :, :] = β * (sin.(curl_shifted_2d(ϕ, 2, 1, -[0, 1])) - sin.(curl_shifted_2d(ϕ, 2, 1, [0, 0])))
    dSdϕ[2, :, :] = β * (sin.(curl_shifted_2d(ϕ, 1, 2, -[1, 0])) - sin.(curl_shifted_2d(ϕ, 1, 2, [0, 0])))

    return dSdϕ
end

#3D
function zero_lattice_3d(Lx, Ly, Lz)
    return zeros(Float64, (3, Lx, Ly, Lz))
end

function rand_lattice_3d(Lx, Ly, Lz)
    return π * (2 * rand(Float64, (3, Lx, Ly, Lz)) .- 1)
end

function curl_shifted_3d(ϕ, μ, ν, k)
    # returns an array A such that A[n] = θ_μν(n+k) where n is a lattice site
    if μ == ν
        return 0
    else
        return sa.circshift(ϕ[ν, :, :, :], -(k + unit_vectors_3d[μ])) - sa.circshift(ϕ[ν, :, :, :], -k) - sa.circshift(ϕ[μ, :, :, :], -(k + unit_vectors_3d[ν])) + sa.circshift(ϕ[μ, :, :, :], -k)
    end
end

function S_3d(ϕ, β)
    β * sum(3 .- (cos.(curl_shifted_3d(ϕ, 2, 1, [0, 0, 0])) + cos.(curl_shifted_3d(ϕ, 3, 1, [0, 0, 0])) + cos.(curl_shifted_3d(ϕ, 3, 2, [0, 0, 0]))))
end

function dSdϕ_3d(ϕ, β)

    dSdϕ = similar(ϕ)
    dSdϕ[1, :, :, :] = β * (sin.(curl_shifted_3d(ϕ, 2, 1, -[0, 1, 0])) - sin.(curl_shifted_3d(ϕ, 2, 1, [0, 0, 0])) + sin.(curl_shifted_3d(ϕ, 3, 1, -[0, 0, 1])) - sin.(curl_shifted_3d(ϕ, 3, 1, [0, 0, 0])))
    dSdϕ[2, :, :, :] = β * (sin.(curl_shifted_3d(ϕ, 1, 2, -[1, 0, 0])) - sin.(curl_shifted_3d(ϕ, 1, 2, [0, 0, 0])) + sin.(curl_shifted_3d(ϕ, 3, 2, -[0, 0, 1])) - sin.(curl_shifted_3d(ϕ, 3, 2, [0, 0, 0])))
    dSdϕ[3, :, :, :] = β * (sin.(curl_shifted_3d(ϕ, 1, 3, -[1, 0, 0])) - sin.(curl_shifted_3d(ϕ, 1, 3, [0, 0, 0])) + sin.(curl_shifted_3d(ϕ, 2, 3, -[0, 1, 0])) - sin.(curl_shifted_3d(ϕ, 2, 3, [0, 0, 0])))
    return dSdϕ
end
end
