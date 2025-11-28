# Plastic models

In PeriLab J2 plasticity model (von Mises plasticity) is implemented for correspondence-based and ordinary state-based Peridynamics. The plastic model was taken from [Peridigm](https://github.com/peridigm/peridigm/blob/master/src/materials/elastic_plastic_correspondence.cxx) and translated into julia. The theoretical basis can be found here [LammiC2014](@cite).

!!! info "Isotropic model"
This model works only for isotropic yield stresses.

## Fundamental Theory

### J2 Plasticity (von Mises Plasticity)

The J2 plasticity theory is based on the von Mises yield criterion, which states that yielding occurs when the von Mises stress reaches a critical value (yield stress). This theory is particularly suitable for ductile materials like metals.

### Stress Decomposition

The total stress tensor **σ** is decomposed into:

$$\boldsymbol{\sigma} = \boldsymbol{s} + p\mathbf{I}$$

where:

- **s** is the deviatoric stress tensor
- _p_ is the spherical (hydrostatic) stress
- **I** is the identity tensor

The spherical stress is calculated as:

$$p = \frac{1}{3}\text{tr}(\boldsymbol{\sigma}) = \frac{1}{3}(\sigma_{11} + \sigma_{22} + \sigma_{33})$$

The deviatoric stress tensor is:

$$\boldsymbol{s} = \boldsymbol{\sigma} - p\mathbf{I}$$

### von Mises Stress

The von Mises equivalent stress is defined as:

$$\sigma_{vM} = \sqrt{\frac{3}{2}\boldsymbol{s}:\boldsymbol{s}} = \sqrt{\frac{3}{2}s_{ij}s_{ij}}$$

This can also be written as:

$$\sigma_{vM} = \sqrt{3J_2}$$

where _J₂_ is the second invariant of the deviatoric stress tensor.

### Yield Criterion

The yield condition is:

$$f = \sigma_{vM} - \sigma_y \leq 0$$

where:

- _f_ is the yield function
- σᵧ is the yield stress (possibly reduced by flaw functions)

**Elastic regime**: If _f_ < 0, the material remains elastic and stresses are unchanged.

**Plastic regime**: If _f_ ≥ 0, plastic deformation occurs and stresses must be returned to the yield surface.

## Return Mapping Algorithm

When the trial stress exceeds the yield surface, a radial return mapping is applied to the deviatoric stress:

### Step 1: Calculate Trial Deviatoric Stress Magnitude

$$\|\boldsymbol{s}^{trial}\| = \frac{\sigma_{vM}^{trial}}{\sqrt{\frac{2}{3}}}$$

### Step 2: Scale Deviatoric Stress to Yield Surface

$$\boldsymbol{s}^{t+\Delta t} = \boldsymbol{s}^{trial} \cdot \frac{\sqrt{\frac{2}{3}} \cdot \sigma_y}{\|\boldsymbol{s}^{trial}\|}$$

### Step 3: Reconstruct Total Stress

$$\boldsymbol{\sigma}^{t+\Delta t} = \boldsymbol{s}^{t+\Delta t} + p^{t+\Delta t}\mathbf{I}$$

Note: The spherical stress remains unchanged during plastic return, as J2 plasticity assumes plastic incompressibility (volume preservation).

## Equivalent Plastic Strain Update

The equivalent plastic strain is updated using an incremental approach that is independent of the specific return mapping algorithm:

### Deviatoric Strain Increment

The deviatoric part of the strain increment is:

$$\boldsymbol{e}^{dev} = \boldsymbol{\varepsilon}^{inc} - \frac{1}{3}\text{tr}(\boldsymbol{\varepsilon}^{inc})\mathbf{I}$$

### Plastic Strain Increment Calculation

The plastic strain increment is computed by projecting onto the flow direction:

$$\mathbf{A} = \boldsymbol{e}^{dev} - \frac{\boldsymbol{s}^{t+\Delta t} - \boldsymbol{s}^t}{2G}$$

$$\mathbf{B} = \frac{1}{2}\left(\frac{\boldsymbol{s}^{t+\Delta t}}{\|\boldsymbol{s}^{t+\Delta t}\|} + \frac{\boldsymbol{s}^t}{\|\boldsymbol{s}^t\|}\right)$$

$$\Delta \varepsilon_p^{eq} = \max\left(0, \sqrt{\frac{2}{3}} \, \mathbf{A}:\mathbf{B}\right)$$

where:

- _G_ is the shear modulus
- The contraction **A**:**B** represents the tensor inner product

### Total Equivalent Plastic Strain

$$\varepsilon_p^{eq,t+\Delta t} = \varepsilon_p^{eq,t} + \Delta \varepsilon_p^{eq}$$

## Material Parameters

The model requires the following material parameters:

1. **Shear Modulus** (_G_): Controls elastic shear response and plastic strain increment calculation
2. **Yield Stress** (σᵧ): Critical stress at which plastic deformation begins

## Flaw Function

The yield stress can be spatially reduced using a flaw function:

$$\sigma_y^{reduced} = f_{flaw}(\mathbf{x}) \cdot \sigma_y$$

This allows modeling of material defects, damage, or heterogeneity.
