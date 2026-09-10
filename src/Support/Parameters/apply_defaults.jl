# SPDX-License-Identifier: BSD-3-Clause
# Call this once, right after validate_params succeeds, so accessors can
# assume every SField with a default is actually present.

function apply_defaults!(node::SObject, dict::AbstractDict)
    for (key, child) in node.fields
        if !haskey(dict, key)
            if child isa SField && child.default !== nothing
                dict[key] = child.default
            end
            continue
        end
        if child isa SObject || child isa SAny
            apply_defaults!(child, dict[key])
        end
    end
end

function apply_defaults!(node::SAny, dict::AbstractDict)
    for (_, entry) in dict
        apply_defaults!(node.value, entry)
    end
end

apply_defaults!(::SchemaNode, ::Any) = nothing  # SField/SOneOf/SArray: no-op here
