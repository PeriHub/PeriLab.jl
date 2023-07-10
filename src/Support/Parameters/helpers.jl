function dof_index(index, inds)
    dest = Vector{eltype(index)}(undef, inds * length(index))
    i = 1
    @inbounds for a in index, b in 1:inds
        dest[i+=1] = a * inds + b - inds
    end
    dest
end

function find_indices(vector, what)
    return findfirst(item -> item == what, vector)
end


function get_blockNodes(blockID)
    maxBlock = maximum(blockID)
    blockNodes = distribution = [collect(1:maxBlock)]
    for i in 1:maxBlock
        blockNodes[i] = find_indices(blockID, i)
    end
    return blockNodes
end