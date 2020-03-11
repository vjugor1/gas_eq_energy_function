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
        edge_dict[k] = (string(v["f_junction"]), string(v["t_junction"]), v["length"] / 100000.0)
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
    return j_list_clean, e_dict
end


function get_inc_matrix(junction_list, edge_list)
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
    return A
end

function get_adj_matrix(junction_list, edge_dict)
    A = zeros(length(junction_list), length(junction_list))
    for (idx_e, (k, v)) in enumerate(edge_dict)
        #v[1], v[2]
        idx_1 = findall(x->x==v[1], junction_list)[1]
        idx_2 = findall(x->x==v[2], junction_list)[1]
        A[idx_1, idx_2] = 1
        A[idx_2, idx_1] = -1
    end
    return A
end

function write_matrix_to_file(file_name, A)
    open(file_name, "w") do io
        i_, j_ = size(A)
        for i=1:i_
            for j=1:j_
                if j == j_
                    write(io, string(string(Int(A[i, j])), "\n") )
                else
                    write(io, string(string(Int(A[i, j])), ",") )
                end
            end
        end
    end
end
