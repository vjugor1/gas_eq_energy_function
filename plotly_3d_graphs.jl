#Functions for 3d plotting of the solution
using PlotlyJS



function topo_surface(data, title_name_)
    z = data
    trace = PlotlyJS.surface(z=z)
    layout = Layout(title=title_name_, autosize=false, width=500,
                    height=500, margin=attr(l=65, r=50, b=65, t=90))
    PlotlyJS.plot(trace, layout)
end
