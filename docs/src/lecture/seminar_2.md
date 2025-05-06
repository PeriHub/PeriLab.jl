## Lecture 2: Fracture

Multiple ways to deal with


## Critical stretch
The critical stretch model is widely used in literature [BobaruF2016](@cite), [MadenciE2014](@cite). It defines the critical length change, or stretch $s$ as a criterion for a damage.

$$s_{crit}\leq s \frac{| \underline{\mathbf{Y}} |}{| \underline{\mathbf{X}} |}$$

$$s_{crit} = \sqrt{\frac{G_{0C}}{[3G+(\frac{3}{4})^4(K-\frac{5G}{3})]\delta}}$$

## Critical energy

The critical energy model introduced by [FosterJT2011](@cite) is valid for state-based peridynamic analysis.
The bond energy is defined as:

$$w_{bond} = \int_{\boldsymbol{\eta}_{final}} (\mathbf{\underline{T}}[x,t]\langle x'-x\rangle - \mathbf{\underline{T}}[x',t]\langle x-x'\rangle)d\boldsymbol{\eta}$$

with the relative displacement vector as:

$$\boldsymbol{\eta}=\mathbf{\underline{u}}[x',t]-\mathbf{\underline{u}}[x,t]$$

If the bond energy is bigger than or equal to the critical energy value, then the bond is considered to be broken:

$$w_{crit} \leq w_{bond}$$

The critical bond energy can be defined as:

$$w_{crit} = \frac{4G_{0C}}{\pi\delta^4}$$

---

## Differences
- discussion in the seminar
    - direction
    - experimental

---
# Numerical integration
**Question**

How many bonds are allowed to break per time step?

Size of the horizon?

**Integration of damage model**

```Julia
## Code concept

#Run material
#Evaluate Damage
#Run Material

omega = ones(length(nodes), nneighbor)


for idt in time

    K = create_K(omega)

    u = inv(K)*F_V # or from time integration

    for iID in nodes
        for jID in nlist[iID]
            s = (u[iID]+x[iID]-(u[jID]+x[jID]))/(x[jID]-x[iID])
            if  s>s_crit
                omega[i,j]=0
            end
        end
    end

    K = create_K(omega)

    u = inv(K)*F_V # or from time integration

    ...
end
```
