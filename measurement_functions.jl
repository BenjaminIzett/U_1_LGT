module MeasurementFunctions
include("u_1_lgt.jl")
include("wilson_loops.jl")

export *

function mean_plaquette_3d(ϕ, β)
    1 - U_1_LGT.S_3d(ϕ, β) / (length(ϕ) * β)
end

function mean_plaquette_2d(ϕ, β)

    1 - U_1_LGT.S_2d(ϕ, β) / (length(ϕ) * β)

end

function hamer_loop_range_3d(ϕ, β, α, N_smears, R_τ_pairs)
    ϕ_ape_smeared = WilsonLoops.N_ape_smearing_3d(exp.(complex.(0, ϕ[1:2, :, :, :])), α, N_smears)
    U_t = WilsonLoops.parisi_dir_3d(ϕ, β, 3)


    [WilsonLoops.hamer_loop_3d(ϕ, ϕ_ape_smeared, U_t, R, τ, β, α, N_smears) for (R, τ) in R_τ_pairs]
end

function parisi_loop_range_3d(ϕ, β, R_τ_pairs)

    U = exp.(complex.(0, ϕ))
    U_t = WilsonLoops.parisi_3d(ϕ, β)

    [WilsonLoops.parisi_3d_loops(U, U_t, R, τ, β) for (R, τ) in R_τ_pairs]
end

function flux_tubes_x_3d(ϕ, τs)
    return WilsonLoops.flux_tubes_x(ϕ, τs)
end

end
