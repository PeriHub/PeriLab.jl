# Influence Function

The influence function `ω` weights the contribution of each bond within the
horizon. It can be defined either by selecting a **predefined name** or by
providing a **custom mathematical expression**.

## Predefined Functions
The following predefined options are available:
- **1/xi^2**: Classical inverse-square influence function based on the bond
  length `xi`.
- **1**: Constant influence function (all bonds weighted equally).

## Custom Expression
Instead of a predefined name, an arbitrary mathematical expression can be
provided as a string. The following variables are available inside the
expression:
- **xi**: Bond length (Float64).
- **xiX/xiY/xiZ**: Bond direction components (Float64). `xiZ` is only
  meaningful in 3D and defaults to `0.0` in 2D.

Any valid Julia expression using these variables (e.g. `+`, `-`, `*`, `/`,
`^`, `exp`, `sqrt`, ...) can be used.

The yaml definition looks like this
```yaml
Discretization:
  Influence Function: "1/xi^2"
```
or, using a custom expression
```yaml
Discretization:
  Influence Function: "exp(-xi/3)"
```
or
```yaml
Discretization:
  Influence Function: "xiX^2 + xiY^2"
```
