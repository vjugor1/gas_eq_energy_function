using JSON
path_belgian = "gas_eq_energy_function/belgian_lanl_ansi.json"
function read_json!(out, file_name)
    dict = Dict()
    open(file_name, "r") do f
        global dict = JSON.parse(f)  # parse and transform data
    global out = dict
    end
    nothing
end

function get_j_list_e_dict(data_dict)
    junction_list = []
    for (k, v) in data_dict["junction"]
       push!(junction_list, k)
    end
    edge_dict = Dict()
    for (k, v) in data_dict["pipe"]
        edge_dict[k] = [string(v["f_junction"]), string(v["t_junction"]), v["length"] / 100000.0]
    end
    return junction_list, edge_dict
end


function clean_belgian_leave_pipes_only(j_list, e_dict)
    j_to_clean_idx = [1, 24, 22, 6, 19, 18, 12, 10]
    j_list_clean = j_list[setdiff(1:end, j_to_clean_idx)]
    j_list_wiped = setdiff(j_list, j_list_clean)
    e_to_delete = []
    for (k, v) in e_dict
        if ((v[1] in j_list_wiped) | (v[2] in j_list_wiped))
            push!(e_to_delete, k)
        end
    end
    for item in e_to_delete
        delete!(e_dict, item)
    end
    for (k, v) in e_dict
        for (k_, v_) in e_dict
            if ((k != k_) & (v == v_))
                delete!(e_dict, k_)
            end
        end
    end
    return j_list_clean, e_dict
end

function get_clean_belgian_j_e()
    out = Dict()
    read_json!(out, path_belgian)
    println(out)
    j_list, e_dict = get_j_list_e_dict(dict)
    j_list, e_dict = clean_belgian_leave_pipes_only(j_list, e_dict)
    e_tmp = find_edge(e_dict, "11", "17")
    e_dict[e_tmp][3] = 1.35
    return j_list, e_dict
end
