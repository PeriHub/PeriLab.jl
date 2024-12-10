# Damage Models Theory

## Critical stretch
The critical stretch model is widely used in literature [BobaruF2016](@cite), [MadenciE2014](@cite). It defines the critical length change, or stretch $s$ as a criterion for a damage.

$$s_{crit}\leq s \frac{| \underline{\mathbf{Y}} |}{| \underline{\mathbf{X}} |}$$

The advantage of this criterion is, that the implementation is rather simple. Also the result is purely geometrical and therefore not influenced by the origin of the bond, because it is neighborhood independed. However, for complex load cases it is to simple and won't work well.

Some literature describes the possibility to compute the critical stretch based on the energy release rate.
