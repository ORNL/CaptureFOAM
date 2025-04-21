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
    Foam::interphaseMassTransferModels::physical

Description
    Enforces a physical mass transfer rate across a bulk-film interface.
    The overall mass transfer rate is calculated as the harmonic mean of
    the bulk-side mass transfer coefficient k_b and the film-side mass transfer
    coefficient k_l. Presently, the only option for k_b is a constant value,
    and k_l is calculated from Higbie's slip penetration theory.

SourceFiles
    physical.C

\*---------------------------------------------------------------------------*/

#ifndef physical_H
#define physical_H

#include "interphaseMassTransferModel.H"
#include "massTransferRateCoefficientModel.H"
#include "enhancementModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interphaseMassTransferModels
{

/*---------------------------------------------------------------------------*\
                               Class physical
\*---------------------------------------------------------------------------*/

class physical
:
    public interphaseMassTransferModel
{
    // Private Data

        //- Coeff1 for Henry constant for species and solvent (Pa-m^3/mol/K)
        const dimensionedScalar H1_;

	//- Coeff2 for Henry constant for species and solvent (Pa-m^3/mol)
        const dimensionedScalar H2_;

        //- Reference to film-side mass transfer rate coefficient model
        autoPtr<massTransferRateCoefficientModel> filmRateModel_;
        
        //- Reference to bulk-side mass transfer rate coefficient model
        autoPtr<massTransferRateCoefficientModel> bulkRateModel_;
        
        //- Reference to enhancement factor model
        autoPtr<enhancementModel> enhancementModel_;
        
        //- Interfacial mass transfer rate
        volScalarField N_;

public:

    //- Runtime type information
    TypeName("physical");


    // Constructors

        //- Construct from components
        physical
        (
            const dictionary& dict,
            const solvers::multicomponentFluid& fluid,
            const label& bulkPatchID,
            const mappedPatchBase& patchMap
        );


    //- Destructor
    virtual ~physical()
    {}


    // Member Functions
    
        //- Update mass transfer field
        virtual void update();

        //- Read the dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interphaseMassTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
