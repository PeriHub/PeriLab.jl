# Degradation Models

Degradation models change the material over time, e.g. by removing nodes. A block uses one
with `Degradation Model`; the solver evaluates them with `Degradation Models: true`.

## Thermal Decomposition

A node decomposes once its temperature reaches the `Decomposition Temperature`: it is
deactivated (`Active` is set to false) and all its bonds, in both directions, are broken. A
decomposed node stays inactive. The temperature comes from the thermal models, so they have to
be active in the same step.

```@eval
using PeriLab, Markdown
Markdown.parse(PeriLab.parameter_table("Thermal Decomposition"))
```

## Bond-based Corrosion

!!! warning "Not implemented"
    The model is registered, but does not compute anything yet.
