#Functions for 2d plotting output

using Plots
using ArgCheck
using PlotlyJS

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

function plot_pipe_x_debug(array, pipe_num, time_frame)
    Plots.plot(array[(1 + (pipe_num - 1) * I_) : (pipe_num * I_), time_frame])
end

function plot_pipe_t_debug(array, pipe_num, start=true)
    if start == true
        Plots.plot(array[1 + (pipe_num - 1) * I_, :])
    else
        Plots.plot(array[ (pipe_num ) * I_, :])
    end
end
