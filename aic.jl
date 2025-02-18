module AIC
export *

using CurveFit

function fit_window(fit_type, fit_args, x, y, y_errors, n_params)
    fits = []
    for lower in 1:length(x)
        for upper in (lower+n_params):length(x)
            window = lower:upper
            n = length(x) - length(window)
            fit = curve_fit(fit_type, x[window], y[window], fit_args...)
            y_fit = fit.(x[window])
            χ2 = sum(((y[window] .- y_fit) .^ 2) ./ (y_errors[window] .^ 2))
            prob = exp(-0.5 * χ2 - n_params - n)
            # push!(fits, (prob, χ2, fit, y_fit, window))
            push!(fits, (prob, fit))
        end
    end
    # display(first.(fits))
    N = sum(first.(fits))
    # norm_fits = map(fit -> (fit[1] / N, fit[2:end]...), fits)
    return sum([fit[1] * fit[2] / N for fit in fits])
    # return fits[argmax(first.(fits))][1:end]
end

end
