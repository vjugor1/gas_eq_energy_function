#gen simulated large-scale data

function get_e_dict(depth = 3)
    j_list = []
    n_nodes = sum(2^(collect(1:depth)[i] - 1) for i=1:depth)
    for i=1:sum(n_nodes)
        push!(j_list, string(i))
    end
    e_dict = Dict()
    lay_list = []
    for (inv_lay, lay) in enumerate(collect(1:depth)[end:-1:1])
        pair_cnt = 1
        println("lay=", lay)
        println("inv_lay=", inv_lay)
        curr_node_list = []
        for i=(1:(2^(lay-1)))
            #println("i=",i)

            if inv_lay == 1
                println("i=",i + (inv_lay-1 ) * 2^(lay) )
                push!(curr_node_list, i + (inv_lay-1 ) * 2^(lay) )
            else
                println("i=",i + sum(2^(huh) for huh in collect((lay):(depth))[1:end-1]) )
                push!(curr_node_list, i + sum(2^(huh) for huh in collect((lay):(depth))[1:end-1]))
            end

        end
        push!(lay_list, curr_node_list)

    end
    for i=1:(length(lay_list)-1)
        pair_cnt = 1
        for j=1:length(lay_list[i])
            if (lay_list[i][j] % 2) == 1
                e_dict[string("lay:", string(lay_list[i][j]))] = [string(lay_list[i][j]), string(lay_list[i+1][pair_cnt]), 1, i, i+1, -sqrt(1/8.57)]
                e_dict[string("lay:", string(lay_list[i][j+1]))] = [string(lay_list[i][j+1]), string(lay_list[i+1][pair_cnt]), 1, i, i+1, -sqrt(1/8.57)]
                pair_cnt = pair_cnt + 1
            end
        end
    end
    return j_list, e_dict
end
