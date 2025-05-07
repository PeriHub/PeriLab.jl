# Damage Models Theory

## Critical stretch
The critical stretch model is widely used in literature [BobaruF2016](@cite), [MadenciE2014](@cite). It defines the critical length change, or stretch $s$ as a criterion for a damage.

$$s_{crit}\leq s \frac{| \underline{\mathbf{Y}} |}{| \underline{\mathbf{X}} |}$$

The advantage of this criterion is, that the implementation is rather simple. Also the result is purely geometrical and therefore not influenced by the origin of the bond, because it is neighborhood independent. However, for complex load cases it is to simple and won't work well.

Some literature describes the possibility to compute the critical stretch based on the energy release rate.

$$s_{crit} = \sqrt{\frac{G_{0C}}{[3G+(\frac{3}{4})^4(K-\frac{5G}{3})]\delta}}$$

## Critical energy

The critical energy model introduced by [FosterJT2011](@cite) is valid for state-based peridynamic analysis.
The bond energy is defined as:

$$w_{bond} = \int_{\boldsymbol{\eta}_{final}} (\mathbf{\underline{T}}[x,t]\langle x'-x\rangle - \mathbf{\underline{T}}[x',t]\langle x-x'\rangle)d\boldsymbol{\eta}$$

with the relative displacement vector as:

$$\boldsymbol{\eta}=\mathbf{\underline{u}}[x',t]-\mathbf{\underline{u}}[x,t]$$

If the bond energy is greater than or equal to the critical energy value, then the bond is considered to be broken:

$$w_{crit} \leq w_{bond}$$

Following [FosterJT2011](@cite) the critical bond energy can be defined as:

$$w_{crit} = \frac{4G_{0C}}{\pi\delta^4}$$

!!! info "Direction"
    Because the bond energy potential can be different in $ij$ compared to $ji$, the stiffness matrix can become unsymmetric. This is a difference to the critical stretch model. However, the critical stretch varies this might occur as well.
