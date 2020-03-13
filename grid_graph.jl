#Functions for grid/plotly_3d_graphs
using LightGraphs
using SimpleWeightedGraphs
include("read_json.jl")
function get_inc_matrix(junction_list, edge_dict)
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

function get_dist_matrix(junction_list, edge_dict)
    A = zeros(length(junction_list), length(junction_list))
    for (idx_e, (k, v)) in enumerate(edge_dict)
        #v[1], v[2]
        idx_1 = findall(x->x==v[1], junction_list)[1]
        idx_2 = findall(x->x==v[2], junction_list)[1]
        A[idx_1, idx_2] = v[3]
        A[idx_2, idx_1] = v[3]
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

function get_src_dest_weights_lists(j_list, e_dict)
    srcs = []
    dest = []
    wght = []
    for (k, v) in e_dict
        src_idx = findall(x->x==v[1], j_list)[1]
        dst_idx = findall(x->x==v[2], j_list)[1]
        push!(srcs, src_idx)
        push!(dest, dst_idx)
        push!(wght, v[3])
    end
    return Array{Int, 1}(srcs), Array{Int, 1}(dest), Array{Float64, 1}(wght)
end

function find_edge(edges_dict, f_junction, t_junction)
    for (k, v) in edges_dict
        if (v[1] == f_junction) & (v[2] == t_junction)
            return k
        end
    end
end

function get_junc_num(junc_name, j_list)
    return findall(x->x==junc_name, j_list)[1]
end

function compute_throughput(p_i_2, p_j_2, alpha_, L_ij)
    return sqrt(abs(p_i_2 - p_j_2) / (alpha_ * L_ij))
end
function compute_right_pressure_2(p_i_2, Q_ij, alpha_, L_ij)
    return p_i_2 - alpha_ * Q_ij^2 * L_ij
end

function get_weighted_graph(j_list, e_dict)
    srcs, dsts, wghts = get_src_dest_weights_lists(j_list, e_dict)
    g = SimpleWeightedDiGraph(srcs, dsts, wghts)
    return g
end

function get_path(g, j_list, e_dict, from, to)
    idx_f = get_junc_num(from, j_list)
    idx_t = get_junc_num(to, j_list)
    path_f_t = enumerate_paths(dijkstra_shortest_paths(g, idx_f), idx_t)
    return path_f_t
end

function convert_path_idx_to_names(path, j_list)
    return j_list[path]
end

function get_path_length(path_f_t)
    length_f_t = 0.0
    #println(path_f_t)
    for i=1:(length(path_f_t)-1)
        edge = find_edge(e_dict, j_list[path_f_t[i]], j_list[path_f_t[i+1]])
        length_f_t = round(length_f_t + e_dict[edge][3], digits=3)
    end
    return length_f_t
end

function set_p_Q_on_path!(e_dict_p_Q, Q, start_p_2, path_names, alpha_)
    for i=1:(length(path_names)-1)
        edge = find_edge(e_dict_p_Q, path_names[i], path_names[i+1])
        if i==1
            push!(e_dict_p_Q[edge], start_p_2)
            p_2_end = compute_right_pressure_2(start_p_2, Q, alpha_, e_dict_p_Q[edge][3])
            push!(e_dict_p_Q[edge], p_2_end)
            push!(e_dict_p_Q[edge], Q)
        else
            edge_prev = find_edge(e_dict, path_names[i-1], path_names[i])
            push!(e_dict_p_Q[edge], e_dict_p_Q[edge_prev][5])
            p_2_end = compute_right_pressure_2(e_dict_p_Q[edge_prev][5], Q, alpha_, e_dict_p_Q[edge][3])
            push!(e_dict_p_Q[edge], p_2_end)
            push!(e_dict_p_Q[edge], Q)
        end
    end
    nothing
end

function set_stat_p_Q(j_list, e_dict, alpha_)
    g = get_weighted_graph(j_list, e_dict)
    e_dict_p_Q = deepcopy(e_dict)
    l_list = []
    Q_list = []
    #"1"-"4" line
    path_1_4 = get_path(g, j_list, e_dict, "1", "4")
    l_1_4 = get_path_length(path_1_4)
    push!(l_list, l_1_4)
    Q_1_4 = compute_throughput(0.96^2, 0.925^2, alpha_, l_1_4)
    path_1_4_names = convert_path_idx_to_names(path_1_4, j_list)
    set_p_Q_on_path!(e_dict_p_Q, Q_1_4, 0.96^2, path_1_4_names, alpha_)
    push!(Q_list, Q_1_4)
    #"5"-"4" line
    path_5_4 = get_path(g, j_list, e_dict, "5", "4")
    path_5_4_names = convert_path_idx_to_names(path_5_4, j_list)
    l_5_4 = get_path_length(path_5_4)
    push!(l_list, l_5_4)
    Q_5_4 = compute_throughput(0.96^2, 0.925^2, alpha_, l_5_4)
    set_p_Q_on_path!(e_dict_p_Q, Q_5_4, 0.96^2, path_5_4_names, alpha_)
    push!(Q_list, Q_5_4)
    #"81"-"11" line
    path_81_11 = get_path(g, j_list, e_dict, "81", "11")
    path_81_11_names = convert_path_idx_to_names(path_81_11, j_list)
    l_81_11 = get_path_length(path_81_11)
    push!(l_list, l_81_11)
    Q_81_11 = compute_throughput(0.827^2, 0.805^2, alpha_, l_81_11)
    set_p_Q_on_path!(e_dict_p_Q, Q_81_11, 0.827^2, path_81_11_names, alpha_)
    push!(Q_list, Q_81_11)
    #"4"-"14" line
    path_4_14 = get_path(g, j_list, e_dict, "4", "14")
    path_4_14_names = convert_path_idx_to_names(path_4_14, j_list)
    l_4_14 = get_path_length(path_4_14)
    push!(l_list, l_4_14)
    Q_4_14 = Q_1_4 + Q_5_4
    set_p_Q_on_path!(e_dict_p_Q, Q_4_14, 0.925^2, path_4_14_names, alpha_)
    push!(Q_list, Q_4_14)
    edge_4_14 = find_edge(e_dict_p_Q, "4", "14")
    #"11"-"14" line
    path_11_14 = get_path(g, j_list, e_dict, "11", "14")
    path_11_14_names = convert_path_idx_to_names(path_11_14, j_list)
    l_11_14 = get_path_length(path_11_14)
    push!(l_list, l_11_14)
    p_14_2 = e_dict_p_Q[edge_4_14][5]
    Q_11_14 = compute_throughput(0.805^2, p_14_2, alpha_, l_11_14)
    set_p_Q_on_path!(e_dict_p_Q, Q_11_14, 0.805^2, path_11_14_names, alpha_)
    push!(Q_list, Q_11_14)
    #"11"-"17" line
    path_11_17 = get_path(g, j_list, e_dict, "11", "17")
    path_11_17_names = convert_path_idx_to_names(path_11_17, j_list)
    l_11_17 = get_path_length(path_11_17)
    push!(l_list, l_11_17)
    Q_11_17 = Q_81_11 - Q_11_14
    #println(Q_11_17)
    set_p_Q_on_path!(e_dict_p_Q, Q_11_17, 0.805^2, path_11_17_names, alpha_)
    push!(Q_list, Q_11_17)
    #"14"-"16" line
    path_14_16 = get_path(g, j_list, e_dict, "14", "16")
    path_14_16_names = convert_path_idx_to_names(path_14_16, j_list)
    l_14_16 = get_path_length(path_14_16)
    push!(l_list, l_14_16)
    Q_14_16 = Q_4_14 + Q_11_14
    set_p_Q_on_path!(e_dict_p_Q, Q_14_16, p_14_2, path_14_16_names, alpha_)
    push!(Q_list, Q_14_16)
    return e_dict_p_Q

end
