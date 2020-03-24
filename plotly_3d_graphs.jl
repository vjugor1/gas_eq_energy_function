#Functions for 3d plotting of the solution
using PlotlyJS


function topo_surface(data, title_name_, epsilon_x=1e-1)
    z = data
    I__, M_ = size(data)
    x = collect(1:I__) * epsilon_x
    t = collect(1:M_) * epsilon_t
    trace = PlotlyJS.surface(x=x, y=t, z=z)
    layout = Layout(title=title_name_, autosize=true, width=500,
                    height=500, margin=attr(l=65, r=50, b=65, t=90))
    PlotlyJS.plot(trace, layout)
end


function foo_3d(n)
     topo_surface(Q_prev[(1 + (n - 1) * I_) : (n * I_), :], "Q(x,t)")
end
#{z: z1, type: 'surface'}
