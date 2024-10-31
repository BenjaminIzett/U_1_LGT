module DataHandler
using DelimitedFiles

using Statistics
# export save_field
# export load_field
# export analyse_data
# export load_data
export *

function save_field(filename, ϕ)
    ϕ_size = size(ϕ)
    open(filename, "w") do io
        writedlm(io, [ϕ_size], ' ')
        writedlm(io, collect(Iterators.flatten(ϕ)))
    end
end


function load_field(filename)
    open(filename, "r") do io
        shape_line = readline(io)
        shape = parse.(Int, split(shape_line, ' '))
        reshape(readdlm(io), (shape...))
    end
end

function analyse_data(data_filename, output_filename, indices, functions)
    n_measurements, contents = open(data_filename, "r",) do io
        n_measurements_line = readline(io)
        n_measurements = parse(Int, split(n_measurements_line, ' ')[2])
        contents = readdlm(io, comments=true, comment_char='#')
        n_measurements, contents
    end


    n_points = Int(size(contents)[1] / n_measurements)
    reshaped_contents = permutedims(reshape(contents, (n_measurements, n_points, size(contents)[2])), (2, 1, 3))
    data = map((index, f) -> mapslices(f, reshaped_contents[:, :, index], dims=2), indices, functions)

    open(output_filename, "w") do io
        writedlm(io, hcat(data...))
    end

end
function load_data(filename)
    open(filename, "r") do io
        readdlm(io, comments=true)
    end
end


function autocorrelation_a(a, data)
    N = length(data)
    mean_data = mean(data)
    Γ_0 = var(data)
    Γ_a = (1 / (N - a)) * sum((data[1:N-a] .- mean_data) .* (data[1+a:N] .- mean_data))
    Γ_a / Γ_0
end
function autocorrelation(range, filename)
    DataHandler.analyse_data("measurements/$filename", "analysis/$filename", tuple(1, [3 for _ in range]...), tuple(length, [data -> autocorrelation_a(n, data) for n in range]...))

end
function calc_specific_heat_capacity(energies)
    return mean(energies)^2 - var(energies)
end


function specific_heat_capacity(filename)
    DataHandler.analyse_data("measurements/$filename", "specific_heat_capacity/$filename", (3,), (calc_specific_heat_capacity,))
end
end
