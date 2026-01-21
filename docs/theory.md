---
parent: CaptureFOAM v1.0 User Guide
title: Theory
nav_order: 2
usemathjax: true
---

# Theory

---

## Governing Equations

CaptureFOAM relies on native OpenFOAM models for weakly compressible fluids to represent systems containing a multiphase flow with a gaseous bulk phase (the OpenFOAM `fluid` module) and a liquid film phase (the CaptureFOAM `multicomponentFilm` module, derived from the OpenFOAM `film` module) flowing along solid boundaries. In each fluid, the compressible Navier-Stokes equations are solved to obtain the velocity $$\vec{u}$$ and pressure $p$ fields:

$$
\frac{\partial \left( \rho \vec{u} \right)}{\partial t} + \nabla \cdot \left( \rho \vec{u} \right) = 0
$$

$$
\frac{\partial \left( \rho \vec{u} \right)}{\partial t} + \vec{u} \cdot \nabla \left( \rho \vec{u} \right) = -\nabla p + \nabla \cdot \left( \mu \nabla \vec{u} \right) + \rho \vec{g} + F_{st}
$$

where $$\rho$$ is the fluid density, $$t$$ is the simulation time, $$\mu$$ is the dynamic viscosity, $$\vec{g}$$ is the gravity vector, and $$F_{st}$$ is the interface surface tension force.

In each phase a transport equation for the energy of the system is solved, represented either by the sensible enthalpy $$h$$ or internal energy $$e$$:

$$
\frac{\partial \left( \rho h + K \right)}{\partial t} + \vec{u} \cdot \nabla \left( \rho h + K \right) - \frac{\partial p}{\partial t} + \nabla \cdot q = S_R
$$

where $$K$$ is the kinetic energy of the flow, $$q$$ is the heat flux, and $$S_R$$ is the heat release associated with chemical reactions.

Finally, in each phase, a transport equation is solved for each specie $$Y_i$$ present in the fluid:

$$
\frac{\partial \left( \rho Y_i \right)}{\partial t} + \vec{u} \cdot \nabla \left( \rho Y_i \right) = D_i \nabla^2 Y_i + R_i + \dot{m}_i
$$

where $$D_i$$ is the diffusivity of species $$i$$ in the fluid, $$R_i$$ is the rate of consumption of the species by chemical reactions, and $$\dot{m}_i$$ is the rate of interphase mass transfer.

---

## Interphase Mass Transfer

The rate of mass transfer of a species $$i$$ between phases is related to the interfacial area density $$a_i$$, the molecular weight of the species $$MW_i$$, and the interfacial mass transfer rate coefficient $$N_i$$:

$$
\dot{m}_i = a_i MW_i N_i
$$

Calculation of the interfacial rate coefficient is the main task of the CaptureFOAM code base, and is done through a pair of `fvModels` which apply appropriate source terms to each of the governing equations above, based upon the calculated interfacial rate coefficient. The interfacial rate coefficient $$N_i$$ is related to the overall mass transfer coefficient $$K_i$$, actual partial pressure of the species in the gas phase $$P_i$$, and the equilibrium partial pressure of the species in the gas phase $$P^*$$:

$$
N_i = K_i \left( P_i - P^* \right)
$$

The equilibrium partial pressure is calculated using Henry's law and appropriate solubility/volatility constants for the species and liquid film solvent in question. The overall mass transfer coefficient $$K_i$$ considers effects from both the gas- and liquid-side mass transfer processes, employing a harmonic mean to balance the contributions:

$$
K_i = \left( \frac{RT}{k_g} + \frac{H_i}{E k_l} \right)^{-1}
$$

Here, $$R$$ is the universal gas constant, $$T$$ is the temperature of the gas at the interface, $$H_i$$ is the Henry constant of the specie, $$E$$ is the enhancement factor from chemical reactions, and $$k_g$$ and $$k_l$$ are the gas- and liquid-side mass transfer coefficients, respectively. A variety of options exist for the calculation of the enhancement factor, as well as for each of the mass transfer coefficients. These options will be discussed in the following sections.

---

## Mass Transfer Rate Coefficient Models

The following options are available for the calculation of the mass transfer rate coefficients:

### Higbie

The Higbie model is derived from analytical penetration theory, and is applicable to gas- or liquid-side mass transfer. It prescribes the rate $k$ as follows:

$$
k = 2 \left( \frac{D}{\pi \tau} \right)^{0.5}
$$

where $$\tau$$ is the contact time, defined as:

$$
\tau = \frac{L}{U}
$$

where $$L$$ is the distance of the fluid from the inlet, and $$U$$ is the magnitude of the velocity of the fluid.

The Higbie model is generally best reserved for simple scenarios, such as flow over a vertical or inclined plate.

### Blasius

The Blasius model is derived from boundary layer theory, and relies on two non-dimensional numbers, the Reynolds number $$Re$$ and the Schmidt number $$Sc$$:

$$
Re = \frac{\rho U L}{\mu}
$$

$$
Sc = \frac{\mu}{\rho D}
$$

Here, $$L$$ is the characteristic length of the system. From these quantities, the rate coefficient is calculated as:

$$
k = 0.332 \frac{D}{L} Re^{0.5} Sc^{0.333}
$$

For the special case of pipe flow, the `pipeBlasius` submodel prescribes the rate coefficient as:

$$
k = 3.66 \frac{D}{L}
$$

where $$L$$ is the diameter of the pipe.

### Porous Media

A specialized rate coefficient model is implemented for the gas-side mass transfer coefficient in porous media, following the work of [1]. The gas-side mass transfer coefficient is calculated as:

$$
k_g = C_p \frac{0.021 D_{i,g}}{d_h} \beta \gamma
$$

where $C_p$ is a model coefficient with a default value of 0.021, $$D_{i,g}$$ is the diffusivity of species $$i$$ in the gas phase, $$d_h$$ is the hydraulic diameter of the pores, and $$\beta$$ and $$\gamma$$ are quantities defined as:

$$
\beta = \left( \frac{\mu_g}{\rho_g D_{i,g}} \right)^{0.333}
$$

$$
\gamma = \left( \frac{\rho_g U_g d_h}{\mu_g \left( 1 - \varepsilon \right)} \right)^{0.8}
$$

Where $$\varepsilon$$ is the porosity of the media, on a scale of 0 to 1.

The porous media model is not applicable to calculate the liquid-side mass transfer rate coefficient.

---

## Enhancement Factor Models

The enhancement factor represents physics occuring within the liquid film which are spatially and temporally unresolved, but play an important part in accelerating mass transfer when solvent in the film is reacting with species which are being absorbed from the gaseous phase. These models are generally based upon the value of the Hatta number $$Ha$$, which is defined as:

$$
Ha = \frac{\left( D_{i,l} k_{app} \right)^{0.5}}{k_{i,l}}
$$

where $$D_{i,l}$$ is the diffusivity of species $$i$$ in the liquid film, $$k_{i,l}$$ is the liquid-side mass transfer rate coefficient of the species, and $$k_{app}$$ is the apparent reaction rate of the species in the liquid film. The apparent reaction rate is dependent upon the exact reaction mechanism implemented in the model, and incorporates effects from all reactions present. For a detailed presentation of the reaction mechanisms and apparent reaction rates implemented within CaptureFOAM for the MEA-CO2 and KSAR-CO2 reactions common in direct air capture applications, see reference [2].

Once calculated, the Hatta number is used to calculate the enhancement factor according to one of the following options.

### No Enhancement

For systems in which no reactions occur, the `noEnhancement` model is used to set the enhancement factor to unity:

$$
E = 1
$$

### Low Hatta Number
The `lowHa` model is generally applicable when $$Ha < 1$$, and simply sets the enhancement factor equal to the Hatta number:

$$
E = Ha
$$

### Film Psuedo-First Order
The `filmPseudoFirstOrder` model, taken from reference [3], is generally applicable when $$Ha > 3$$, and uses a hyperbolic tangent function to calculate the enhancement factor from the Hatta number:

$$
E = \frac{Ha}{\tanh \left( Ha \right)}
$$

---

## References

[1] Wang, Z., Gupta, M., Warudkar, S. S., Cox, K. R., Hirasaki, G. J., & Wong, M. S. (2016). Improved CO2 Absorption in a Gas-Liquid Countercurrent Column Using a Ceramic Foam Contactor. Industrial and Engineering Chemistry Research, 55(5), 1387–1400. https://doi.org/10.1021/acs.iecr.5b03600

[2] Kincaid, K., Brandao, F. L., Chuahy, F. D. F., and Nawaz, K., 2025, “Numerical Assessment of Triply Periodic Minimal Surfaces for Direct Air Capture of Carbon Dioxide,” International Journal of Greenhouse Gas Control, 146, p. 104457. https://doi.org/10.1016/j.ijggc.2025.104457.

[3] Putta, K. R., Tobiesen, F. A., Svendsen, H. F., & Knuutila, H. K. (2017). Applicability of enhancement factor models for CO2 absorption into aqueous MEA solutions. Applied Energy, 206, 765-783.
  

