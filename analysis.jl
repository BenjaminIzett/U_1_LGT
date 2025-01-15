include("data_handler.jl")

function file_mean(files, output_filename)
    data = []
    for file in files
        data_i = DataHandler.load_data("analysis/$file")
        push!(data, data_i)
    end
    stacked_data = cat(data..., dims=3)
    mean_data = mean(stacked_data, dims=3)
    std_data = std(stacked_data, dims=3)

    #handle saving
    mean_data[:, :, 1], std_data[:, :, 1]

    open("analysis/mean/$output_filename", "w") do io
        writedlm(io, hcat(mean_data...))
    end
    open("analysis/std/$output_filename", "w") do io
        writedlm(io, hcat(std_data...))
    end
end


#
# int_autocorrelation_time("autocorrelation_1_05_10000.txt")

# file_mean(["autocorrelation/autocorrelation_1_10000_$(nth_run).txt" for nth_run in 1:10], "autocorrelation_mean.txt")
#
