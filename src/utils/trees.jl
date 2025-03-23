using PhyloNetworks

const Tree = HybridNetwork

function load_trees(mcmc_path::String)::Vector{Tree}
    trees = readnexus_treeblock(mcmc_path)
    if 0 < length(trees)
        return trees
    end

    open(mcmc_path) do io
        return map(readnewick, io |> eachline)
    end
end

load_trees_from_newick(newick::String)::Vector{Tree} = [readnewick(newick)]
load_trees_from_newick(newicks::Vector{String})::Vector{Tree} = [readnewick(newick) for newick in newicks]


get_tip_names(tree::Tree)::Vector{String} = tiplabels(tree) |> sort
get_leaf_index_mapping(tree::Tree) = Dict(leaf => i for (i, leaf) in get_tip_names(tree) |> enumerate)
get_leaf_index(tree::Tree, leaf::String) = findfirst(x -> x == leaf, sort(tiplabels(tree)))


"""
A patch for the existing method changing the regex of taxon names to also allow periods.
"""
function PhyloNetworks.readnexus_translatetable(io)
    rx_translate = r"^\s*translate"i
    rx_emptyline = r"^\s*$"
    line = readline(io)
    translate = false
    id2name = Dict{Int,String}()
    while true
        if occursin(rx_translate, line)
            translate = true
            break
        elseif occursin(rx_emptyline, line)
            line = readline(io)
        else
            translate = false
            break
        end
    end
    if translate # then read the table
        rx_end = r"^\s*;"
        rx_idname = r"\s*(\d+)\s+([\.\w]+)[\s]*([,;]?)"
        while true
            line = readline(io)
            occursin(rx_emptyline, line) && continue
            if occursin(rx_end, line)
                line = readline(io)
                break
            end
            m = match(rx_idname, line)
            if isnothing(m)
                @warn "problem reading the translate table at line $line.\nnumbers won't be translated to names"
                translate = false
                break
            end
            push!(id2name, parse(Int, m.captures[1]) => String(m.captures[2]))
            if m.captures[3] == ";"
                line = readline(io)
                break
            end
        end
    end
    return line, translate, id2name
end
