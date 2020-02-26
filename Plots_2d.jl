#Functions for 2d plotting output

using Plots
using ArgCheck

function plot_stuff(arr, x_label_name, y_label_name, save_bool, file_name="dummy", scaling=1.0, plot_line=false)
    tmp = round.(arr, digits=5)
    if plot_line == false
        Plots.plot(Array(0:length(arr)-1) * scaling, arr, seriestype=:scatter,
                xlabel=x_label_name, ylabel=y_label_name, label=y_label_name)
    else
        Plots.plot(Array(0:length(arr)-1) * scaling, arr,
                xlabel =x_label_name, ylabel=y_label_name, label=y_label_name)
    end
    if save_bool == true
        #@argcheck file_name!="dummy"
        println(file_name)
        Plots.savefig(file_name)
    end
    nothing
end
