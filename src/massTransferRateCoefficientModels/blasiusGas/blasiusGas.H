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
    Foam::massTransferRateCoefficientModels::blasiusGas

Description
    Mass transfer rate coefficient model using Blasius theory
    
References

SourceFiles
    blasiusGas.C

\*---------------------------------------------------------------------------*/

#ifndef blasiusGas_H
#define blasiusGas_H

#include "massTransferRateCoefficientModel.H"
#include "multicomponentFluid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferRateCoefficientModels
{

/*---------------------------------------------------------------------------*\
                                 Class blasiusGas
\*---------------------------------------------------------------------------*/

class blasiusGas
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

        //- Point on inlet patch
        const vector inlet_;

        //- Direction of flow
        const vector dir_;

public:

    //- Runtime type information
    TypeName("blasiusGas");


    // Constructors

        //- Construct from components
        blasiusGas
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const label& patchID
        );


    //- Destructor
    virtual ~blasiusGas()
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
