# Damage Models

| Damage Model           | Critical Stretch | Critical Energy |
|------------------------|:----------------:|:---------------:|
| Critical Value         | ✔️| ✔️|
| Interblock Damage      | (✔️)| (✔️)|
| Anisotropic Damage     | (✔️)| (✔️)|

## Critical Stretch
The critical value correspondends to the critical stretch for this model, defined in the [theory manual](@ref "Damage Models Theory").

## Critical Energy

## Interblock Damage

Interlaminar behaviour between different material blocks can be defined using the `Interblock Critical Value` parameter.
If a Bond is crossing a block interface, a user defined critical damage value is applied to the bonds.

As bonds are bidirectional, the critical damage value can be defined for both orientations, for example:

- **Interblock Damage**:
    - Interblock Critical Value 1_2: 0.1
    - Interblock Critical Value 2_1: 0.1

Bonds that aren't crossing a block interface are not affected by the interblock damage model.

![InterBlockDamage](../../assets/InterBlockDamage.svg)
