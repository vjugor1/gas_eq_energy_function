# Write results into csv file
using DelimitedFiles
function write_to_csv(I_, M, p, Q, e_dict_p_Q)
    num_of_pipes = length(e_dict_p_Q)
    list_of_boys = []
    for item in e_dict_p_Q
        push!(list_of_boys, item[2])
    end
    for i=1:num_of_pipes
        arr_p = p[(1 + (i - 1) * (I_)) : (i * I_), :]
        arr_Q = Q[(1 + (i - 1) * (I_)) : (i * I_), :]
        pipe_name = string(list_of_boys[i][1],string("-", list_of_boys[i][2]))
        writedlm( string(pipe_name, "_p.csv"),  arr_p, ',')
        writedlm( string(pipe_name, "_Q.csv"),  arr_Q, ',')
    end
end
function write_lengths_to_csv(e_dict_p_Q)
    list_of_boys = []
    for item in e_dict_p_Q
        push!(list_of_boys, item[2])
    end
    list_of_ls = []
    for fella in list_of_boys
        push!(list_of_ls, round(fella[3], digits=3))
    end
    writedlm( "lengths.csv",  list_of_ls, ',')
end
