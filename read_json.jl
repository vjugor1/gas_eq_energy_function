using JSON
path_belgian = "gas_eq_energy_function/belgian_lanl_ansi.json"
function read_json!(out, file_name)
    dict = Dict()
    open(file_name, "r") do f
        global dict = JSON.parse(f)  # parse and transform data
    println(dict.keys)
    global out = dict
    end
    nothing
end


function get_inc_matrix(data_dict)
    junction_list = []
    for (k, v) in data_dict["junction"]
       push!(junction_list, k)
    end
    edge_dict = Dict()
    for (k, v) in data_dict["connection"]
        edge_dict[k] = (string(v["f_junction"]), string(v["t_junction"]))
    end
    A = zeros(length(junction_list), length(edge_dict))
    for (idx_j, j) in enumerate(junction_list)
        for (idx_e, (k, v)) in enumerate(edge_dict)
            if j == v[1]
                A[idx_j, idx_e] = -1
            end
            if j == v[2]
                A[idx_j, idx_e] = 1
            end
        end
    end
    return junction_list, edge_dict, A
end
