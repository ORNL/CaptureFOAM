/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                Copyright (C) 2023 Oak Ridge National Laboratory                
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::massTransferRateCoefficientModels::porousMedia

Description
    Mass transfer rate coefficient model for porous media using the formulation
    described by Wang et al. (2016).
    
References
    Wang, Z., Gupta, M., Warudkar, S. S., Cox, K. R., Hirasaki, G. J., & Wong,
    M. S. (2016). Improved CO2 Absorption in a Gas-Liquid Countercurrent Column
    Using a Ceramic Foam Contactor. Industrial and Engineering Chemistry
    Research, 55(5), 1387–1400. https://doi.org/10.1021/acs.iecr.5b03600

SourceFiles
    porousMedia.C

\*---------------------------------------------------------------------------*/

#ifndef porousMedia_H
#define porousMedia_H

#include "massTransferRateCoefficientModel.H"
#include "multicomponentFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferRateCoefficientModels
{

/*---------------------------------------------------------------------------*\
                                 Class porousMedia
\*---------------------------------------------------------------------------*/

class porousMedia
:
    public massTransferRateCoefficientModel
{
    // Private Data
    
        //- The multi-component fluid model
        const solvers::multicomponentFluid& multicomponentFluid_;
        
        //- Reference to the transferred specie thermo package
        const fluidMulticomponentThermo& thermo_;
        
        //- Model coefficient
        const scalar C_;

        //- Diffusivity of specie
        const scalar D_;
        
        //- Pore diameter
        const scalar dp_;
        
        //- Porosity
        const scalar eps_; 

public:

    //- Runtime type information
    TypeName("porousMedia");


    // Constructors

        //- Construct from components
        porousMedia
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const label& patchID
        );


    //- Destructor
    virtual ~porousMedia()
    {}


    // Member Functions
    
        //- Calculate and return rate coefficient field
        const volScalarField& k();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace massTransferRateCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
