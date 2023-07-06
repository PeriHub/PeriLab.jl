module Correspondence
()
export shapeTensor



function shapeTensor(data)

    neigborlist = data["Neighborlist"]
    bondDamage = data["Bond Damage"]
    bondGeometry = data["Bond Geometry"] #dx,dy,dz,len
    invK = data["Inverse Shape Tensor"]
    nstatus = data["Node Status"]
    nnodes = data.nnodes
    dof = data.dof

    bondCount = 0
    for iID in 1:nnodes
        shapeTensor = zeros(Float32, dof, dof)
        nneighbors = length(neigborlist[iID])
        if !nstatus
            for jID in 1:nneighbors
                for i in 1:dof
                    for f in 1:dof
                        shapeTensor[i, j] += bondGeometry[(bondCount-1)*dof+i] * bondGeometry[(bondCount-1)*dof+j] * bondDamage[bondCount]
                    end
                end
            end
        end
    end
    invK[:, :, i] = inv(shapeTensor)
end

return data
end