using Tensorial
module Correspondence

export shapeTensor

export defGrad

function shapeTensor(data)

    neigborlist = data["Neighborlist"]
    bondDamage = data["Bond Damage"]
    bondGeometry = data["Bond Geometry"] #dx,dy,dz,len
    invK = data["Inverse Shape Tensor"]
    nstatus = data["Node Status"]
    nnodes = data.nnodes
    dof = data.dof

    bondCount = 0
    for iID in nnodes
        @Tensor shapeTensor = zeros(Float32, dof, dof)
        if !nstatus
            for jID in neigborlist[iID]
                bondCount += 1

                shapeTensor += bondGeometry[(bondCount-1)*dof+i] * bondGeometry[(bondCount-1)*dof+j] * bondDamage[bondCount]

            end
        end
    end
    invK[:, :, i] = inv(shapeTensor)
end

return data
end

function defGrad(data)
    coor = data["Deformation", "NP1"]
    neigborlist = data["Neighborlist"]
    bondDamage = data["Bond Damage", "NP1"]
    bondGeometry = data["Bond Geometry"] #dx,dy,dz,len
    invK = data["Inverse Shape Tensor"]
    defGradN = data["Deformation Gradient", "N"]
    defGradNP1 = data["Deformation Gradient", "NP1"]
    deltaDefGrad = data["Delta Deformation Gradient"]
    nstatus = data["Node Status"]
    nnodes = data.nnodes
    dof = data.dof

    bondCount = 0
    for iID in nnodes
        @Tensor defTensor = zeros(Float32, dof, dof)
        @Tensor invKTensor = invK[:, :, iID]
        nneighbors = length(neigborlist[iID])
        if !nstatus
            for jID in 1:nneighbors
                bondCount += 1
                nID = neigborlist[iID, jID]
                for i in 1:dof
                    for j in 1:dof
                        defTensor[i, j] += (coor[dof*(iID-1)+i] - coor[dof*(nID-1)+i]) * bondGeometry[(bondCount-1)*dof+j] * bondDamage[bondCount]
                    end
                end
                defGradNP1[:, :, iID] = defTensor * invKTensor
                deltaDefGrad[:, :, iID] = defGradNP1[:, :, iID] - defGradN[:, :, iID]
            end
        end
    end
end