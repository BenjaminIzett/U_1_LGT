module Action
using LinearAlgebra
import ArrayPadding
import ShiftedArrays
export *

const unit_vectors_2d = [I[1:2, k] for k in 1:2]
const unit_vectors_3d = [I[1:3, k] for k in 1:3]

function padded2d(ϕ)
    # dims = size(ϕ)
    # ϕ_padded = zeros(typeof(ϕ[1]), (dims[1], dims[2] + 2, dims[3] + 2))
    # ϕ_padded[:, 2:end-1, 2:end-1] = ϕ
    # ϕ_padded[:, 1, 2:end-1] = ϕ[:, end, :]
    # ϕ_padded[:, end, 2:end-1] = ϕ[:, 1, :]
    # ϕ_padded[:, 2:end-1, 1] = ϕ[:, :, end]
    # ϕ_padded[:, 2:end-1, end] = ϕ[:, :, 1]

    # ϕ_padded[:, 1, 1] = ϕ[:, end, end]
    # ϕ_padded[:, 1, end] = ϕ[:, end, 1]
    # ϕ_padded[:, end, 1] = ϕ[:, 1, end]
    # ϕ_padded[:, end, end] = ϕ[:, 1, 1]
    # return ϕ_padded
    ArrayPadding.pad(ϕ, :periodic, (0, 1, 1), (0, 1, 1))
end

function padded3d(ϕ)
    ArrayPadding.pad(ϕ, :periodic, (0, 1, 1, 1), (0, 1, 1, 1))
end

function curl2(ϕ, μ, ν, n)
    if μ == ν
        return 0
    else

        return ϕ[ν, (n + unit_vectors_2d[μ])...] - ϕ[ν, n...] - ϕ[μ, (n + unit_vectors_2d[ν])...] + ϕ[μ, n...]
    end
end


function curl_SAshift_2(ϕ, μ, ν, n)
    if μ == ν
        return 0
    else

        return ShiftedArrays.circshift(ϕ[ν, :, :], -(n + unit_vectors_2d[μ])) - ShiftedArrays.circshift(ϕ[ν, :, :], -n) - ShiftedArrays.circshift(ϕ[μ, :, :], -(n + unit_vectors_2d[ν])) + ShiftedArrays.circshift(ϕ[μ, :, :], -n)
    end
end

#3D
function curl3(ϕ, μ, ν, n)
    if μ == ν
        return 0
    else

        return ϕ[ν, (n + unit_vectors_3d[μ])...] - ϕ[ν, n...] - ϕ[μ, (n + unit_vectors_3d[ν])...] + ϕ[μ, n...]
    end
end


function curl_SAshift_3(ϕ, μ, ν, n)
    if μ == ν
        return 0
    else

        return ShiftedArrays.circshift(ϕ[ν, :, :, :], -(n + unit_vectors_3d[μ])) - ShiftedArrays.circshift(ϕ[ν, :, :, :], -n) - ShiftedArrays.circshift(ϕ[μ, :, :, :], -(n + unit_vectors_3d[ν])) + ShiftedArrays.circshift(ϕ[μ, :, :, :], -n)
    end
end


function S_2d(ϕ, J)
    (Nx, Ny) = size(ϕ)[2:3]
    ϕ_padded = padded2d(ϕ)
    J * sum([1 - cos(curl2(ϕ_padded, μ, ν, [nx, ny])) for μ in 1:2 for ν in μ+1:2 for nx in 2:Nx+1 for ny in 2:Ny+1])
end

function S_v2_2d(ϕ, J)
    (Nx, Ny) = size(ϕ)[2:3]
    ϕ_padded = padded2d(ϕ)
    J * sum([1 - cos(curl2(ϕ_padded, 1, 2, [nx, ny])) for nx in 2:Nx+1 for ny in 2:Ny+1])
end

function S_SAshift_2d(ϕ, J)
    J * sum(1 .- cos.(curl_SAshift_2(ϕ, 1, 2, [0, 0])))
end


function dSdϕ_2d(ϕ, J, μ, n)
    ϕ_padded = padded2d(ϕ)
    μ_b = 3 - μ
    J * (sin(curl2(ϕ_padded, μ_b, μ, 1 .+ (n - unit_vectors_2d[μ_b]))) - sin(curl2(ϕ_padded, μ_b, μ, 1 .+ n)))
end

function dSdϕ_SAshift_2d(ϕ, J)
    dSdϕ = similar(ϕ)
    dSdϕ[1, :, :] = J * (sin.(curl_SAshift_2(ϕ, 2, 1, -[0, 1])) - sin.(curl_SAshift_2(ϕ, 2, 1, [0, 0])))
    dSdϕ[2, :, :] = J * (sin.(curl_SAshift_2(ϕ, 1, 2, -[1, 0])) - sin.(curl_SAshift_2(ϕ, 1, 2, [0, 0])))

    return dSdϕ
end

#3D
function S_3d(ϕ, J)
    (Nx, Ny, Nz) = size(ϕ)[2:4]
    ϕ_padded = padded3d(ϕ)
    J * sum([1 - cos(curl3(ϕ_padded, μ, ν, [nx, ny, nz])) for μ in 1:3 for ν in μ+1:3 for nx in 2:Nx+1 for ny in 2:Ny+1 for nz in 2:Nz+1])
end
function S_v2_3d(ϕ, J)
    (Nx, Ny, Nz) = size(ϕ)[2:4]
    ϕ_padded = padded3d(ϕ)
    J * sum([3 - (cos(curl3(ϕ_padded, 2, 1, [nx, ny, nz])) + cos(curl3(ϕ_padded, 3, 1, [nx, ny, nz])) + cos(curl3(ϕ_padded, 3, 2, [nx, ny, nz]))) for nx in 2:Nx+1 for ny in 2:Ny+1 for nz in 2:Nz+1])
end

# these are the functions that are actually used
function S_SAshift_3d(ϕ, J)
    J * sum(3 .- (cos.(curl_SAshift_3(ϕ, 2, 1, [0, 0, 0])) + cos.(curl_SAshift_3(ϕ, 3, 1, [0, 0, 0])) + cos.(curl_SAshift_3(ϕ, 3, 2, [0, 0, 0]))))
end


function dSdϕ_SAshift_3d(ϕ, J)

    dSdϕ = similar(ϕ)
    dSdϕ[1, :, :, :] = J * (sin.(curl_SAshift_3(ϕ, 2, 1, -[0, 1, 0])) - sin.(curl_SAshift_3(ϕ, 2, 1, [0, 0, 0])) + sin.(curl_SAshift_3(ϕ, 3, 1, -[0, 0, 1])) - sin.(curl_SAshift_3(ϕ, 3, 1, [0, 0, 0])))
    dSdϕ[2, :, :, :] = J * (sin.(curl_SAshift_3(ϕ, 1, 2, -[1, 0, 0])) - sin.(curl_SAshift_3(ϕ, 1, 2, [0, 0, 0])) + sin.(curl_SAshift_3(ϕ, 3, 2, -[0, 0, 1])) - sin.(curl_SAshift_3(ϕ, 3, 2, [0, 0, 0])))
    dSdϕ[3, :, :, :] = J * (sin.(curl_SAshift_3(ϕ, 1, 3, -[1, 0, 0])) - sin.(curl_SAshift_3(ϕ, 1, 3, [0, 0, 0])) + sin.(curl_SAshift_3(ϕ, 2, 3, -[0, 1, 0])) - sin.(curl_SAshift_3(ϕ, 2, 3, [0, 0, 0])))
    return dSdϕ
end
end
