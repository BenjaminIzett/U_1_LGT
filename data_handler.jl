module DataHandler
using DelimitedFiles
using Statistics

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
function save_data(filename, data)
    open(filename, "w") do io
        writedlm(io, data)
    end
end
function save_data_custom(filename, data, mode="w", comment="")
    open(filename, mode) do io
        if (comment != "")
            write(io, "# $comment\n")
        end
        writedlm(io, data)
    end
end


function load_data(filename)
    open(filename, "r") do io
        readdlm(io, comments=true)
    end
end

function load_data_reshaped(filename)
    open(filename, "r") do io
        N_line = readline(io)
        N_measurements = parse(Int, split(N_line, ' ')[2])
        data = readdlm(io, comments=true)
        shape_2d = size(data)
        shape_3d = (N_measurements, Int(shape_2d[1] / N_measurements), shape_2d[2])
        shaped_data = permutedims(reshape(data, shape_3d), (2, 1, 3))
    end
end

end

# function calc_specific_heat_capacity(energies)
#     return mean(energies)^2 - var(energies)
# end

# function specific_heat_capacity(filename)
#     DataHandler.analyse_data("measurements/$filename", "specific_heat_capacity/$filename", (3,), (calc_specific_heat_capacity,))
# end
