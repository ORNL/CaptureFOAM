---
parent: CaptureFOAM Tutorials
nav_order: 1
usemathjax: true
---

# Falling Film Tutorial

## Case Description
The `fallingFilm` tutorial case is the most basic of the examples included with the CaptureFOAM repository. It describes a film flowing down a plate, inclined at 60*, with a gas flowing over top in a counter-flow configuration. The `constant` interphase mass transfer model is used to force the gas phase to transfer CO2 to the film phase at a rate of 0.1 g/m$$^2$$/s. The plate has a width of 0.5 m and a height of 1.0 m, leading to a surface area of 0.5 m$$^2$$ and hence a theoretical total mass transfer of 0.05 g/s of CO2 from the gas phase to the film phase.

## Mesh and Solver Settings
The `blockMesh` utility is used to generate a mesh with a width of 0.5 m in the x-direction, 0.1 m in the y-direction, and 1.0 m in the z-direction. The mesh is divided into 25, 20, and 50 elements in each respective direction with uniform grading. The film mesh is extuded from the wall boundary to a thickness of 1.0 mm. The gravity vector is rotated to simulate a plate with an inclination of 60* from horizontal. Film with zero CO2 content enters from the top boundary and flows freely out the bottom, while gas enters the bottom boundary and exits the top. The sides of the domain employ no-slip conditions for both phases. The simulation is allowed to run for 200 s, which was determined to be sufficient to reach a steady state. During the simulation, the mass flow rate of CO2 across each inflow/outflow boundary is monitored for each phase using a `coded` `functionObject`.

## Results
At steady state, the mass flow rate of CO2 in the film phase at the outlet is 0.05 g/s, and the difference in CO2 mass flow in the gas phase is 0.045 g/s between the inlet and outlet. Both figures are within 10% of the prescribed value of 0.05 g/s.
