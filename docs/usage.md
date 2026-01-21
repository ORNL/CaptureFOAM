---
parent: CaptureFOAM v1.0 User Guide
title: Usage
nav_order: 3
usemathjax: true
---

# Usage

---

## Running a Simulation
CaptureFOAM simulations can be run using any of the included tutorial case as a template. We will demonstrate here using the `meaTube` tutorial. To begin, source OpenFOAM and copy the case files from the tutorial to a new directory:
1. Prepare the case directory structure using a provided template:
   ```bash
   mkdir -p $FOAM_RUN/capturefoam
   cd $FOAM_RUN/capturefoam
   cp -r $WM_PROJECT_USER_DIR/CaptureFOAM/tutorials/meaTube testCase
   cd testCase
   ```

2. If desired, modify the input files to reflect the geometry, flow conditions, and chemical processes present in your problem. Important files are described in [Case Files](#case-files).

3. Generate a mesh, perform appropriate pre-processing steps, decompose your case to run in parallel if desired, and then run CaptureFOAM. These steps are illustrated in the `Allrun` scripts which are present in each tutorial case. For the `meaTube` case, the `Allrun` script contains the following:
   ```bash
    #!/bin/sh
    cd ${0%/*} || exit 1    # Run from this directory
    
    . $WM_PROJECT_DIR/bin/tools/RunFunctions
    
    # Create bulk mesh
    runApplication -s fluid blockMesh -region fluid
    
    # Create film mesh
    runApplication extrudeToRegionMesh -region fluid -overwrite
    
    # Change wedge patches to symmetry in film mesh
    sed -i 's/wedge/symmetry/g' constant/film/polyMesh/boundary
    
    # Run case
    runApplication $(getApplication)
    
    # Post-processing
    runApplication reconstructPar -allRegions
   ```

5. Visualize and post-process the results using [ParaView](https://www.paraview.org/)
   ```bash
   paraFoam -builtin
   ```
   
## CaptureFOAM File Structure

A typical CaptureFOAM case directory is organized as follows:
```
|-- case-directory/
|   |-- 0/
|   |   |-- film/
|   |   |   |-- U
|   |   |   |-- p
|   |   |   |-- T
|   |   |   |-- delta
|   |   |   |-- species_1
|   |   |   |-- species_2
|   |   |   |-- species_3
|   |   |-- fluid/
|   |   |   |-- U
|   |   |   |-- p
|   |   |   |-- p_rgh
|   |   |   |-- T
|   |   |   |-- species_1
|   |   |   |-- species_2
|   |   |   |-- species_4
|   |
|   |-- constant/
|   |   |-- film/
|   |   |   |-- chemistryProperties
|   |   |   |-- combustionProperties
|   |   |   |-- fvModels
|   |   |   |-- g
|   |   |   |-- momentumTransport
|   |   |   |-- physicalProperties
|   |   |   |-- reactions
|   |   |   |-- speciesThermo
|   |   |-- fluid/
|   |   |   |-- fvModels
|   |   |   |-- g
|   |   |   |-- momentumTransport
|   |   |   |-- physicalProperties
|   |   |   |-- speciesThermo
|   |
|   |
|   |-- system/
|   |   |-- film/
|   |   |   |-- fvSchemes
|   |   |   |-- fvSolution
|   |   |-- fluid/
|   |   |   |-- blockMeshDict
|   |   |   |-- extrudeToRegionMeshDict
|   |   |   |-- fvSchemes
|   |   |   |-- fvSolution
|   |   |-- controlDict
|   |   |-- fvSolution
|   |
|   +-- ..
```

### Case Files

- `0`: Directory containing the initial conditions and boundary conditions for each field in each phase. The fields are specified by their file name:
  - `phase_type/U`:            velocity
  - `phase_type/p_rgh`:        reduced pressure
  - `fluid/p`:                 thermodynamic pressure (gas phase only)
  - `phase_type/T`:            temperature
  - `film/delta`:              film thickness (film phase only)
  - `phase_type/species_i`:    mass fraction of species $$i$$ (each phase may have unique species)

- `constant`: Directory containing the definitions for thermophysical properties, turbulence modeling, and reactions for each phase, including:

  - `film/chemistryProperties`: ODE solver settings for the chemical reaction kinetics (film phase only)
  - `film/combustionProperties`: turbulence modeling type for reaction kinetics (film phase only)
  - `film/reactions`: description of reaction mechanism following OpenFOAM's standard format (film phase only)
 
    {: .custom }
    Reactions may be enabled in the gas phase by including the appropriate files, but are typically not used for the applications of interest for CaptureFOAM.
    
  - `film/fvModels`: Optional `fvModels` for film phase. Interphase mass transfer is enabled by setting the following in the `film/fvModels` file:
    ```cpp
    filmMassTransfer
    {
        type                filmMassTransfer;
        
        select              all;
        
        specie              CO2;
    }
    ```
  - `fluid/fvModels`: Optional `fvModels` for fluid phase. The bulk of the definitions for interphase mass transfer are contained in this file. An example is taken from the `meaTube` tutorial case:
    ```cpp
    bulkMassTransfer
    {
        type                bulkMassTransfer;
        
        select              all;
        
        transferPatch       film; // Name of patch on which mass transfer occurs
            
        specie              CO2; // Name of specie to transfer
        
        // Mass transfer model type (options: constant, physical)
        interphaseMassTransferModel physical;
        
        // Coefficients for physical mass transfer model
        physicalCoeffs
        {        
            H1                  0.0; // 1st coefficient for Henry volatility constant, Pa-m^3/kmol/K
    
            H2                  2.27e+6; // 2nd coefficient for Henry volatility constant, Pa-m^3/kmol
            
            Dl1                 0.0; // 1st coefficient for Diffusivity of specie in film, m^2/s/K
    
            Dl2                 1.66e-9; // 2nd coefficient for Diffusivity of specie in film, m^2/s 
            
            bulkMassTransferRateModel
            {
                massTransferRateCoefficientModel Higbie;
                
                HigbieCoeffs
                {            
                    D1              0.0;
                
                    D2           	1.5e-5;
    
                    inlet           (0.0 0.0 0.0); // Point anywhere along bulk inlet
    
                    direction       (0 0 1); // Bulk flow direction
                }
            }
            
            filmMassTransferRateModel
            {
                massTransferRateCoefficientModel Higbie;
                
                HigbieCoeffs
                {            
                   D1               $Dl1;
    
                   D2               $Dl2;
                
                   inlet            (0.0 0.0 1.0); // Point anywhere along film inlet
            
                   direction        (0 0 -1); // Film flow direction
                }
            }
            
            enhancementModel    lowHa; // Enhancement factor model type
            
            tStart              10.0; // When to start chemically enhanced xfer
        }
    }
    ```
  - `phase_type/g`: Definition of gravity vector for each phase
  - `phase_type/momentumTransport`: Sets turbulence model type for each phase
  - `phase_type/physicalProperties`: Mixture type and energy equation settings for each phase
  - `phase_type/speciesThermo`: Thermophysical properties for each specie contained in each phase

- `system`: Contains the discretization and simulation control settings.
  - `controlDict`: Simulation time settings and input/output options.
  - `fvSolution`: The base `fvSolution` file sets the number of outer correctors performed by the multi-region solver
  - `phase_type/fvSchemes`: Finite volume discretization scheme settings for each phase
  - `phase_type/fvSolution`: Linear solver and PIMPLE algorithm settings for each phase
