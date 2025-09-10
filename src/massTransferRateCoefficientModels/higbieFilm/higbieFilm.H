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
    Foam::massTransferRateCoefficientModels::higbieFilm

Description
    Mass transfer rate coefficient model based on Higbie penetration theory.

SourceFiles
    higbieFilm.C

\*---------------------------------------------------------------------------*/

#ifndef higbieFilm_H
#define higbieFilm_H

#include "OFstream.H"
#include "massTransferRateCoefficientModel.H"
#include "multicomponentFilm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferRateCoefficientModels
{

/*---------------------------------------------------------------------------*\
                                 Class higbieFilm
\*---------------------------------------------------------------------------*/

class higbieFilm
:
    public massTransferRateCoefficientModel
{
    // Private Data

        //- The multi-component fluid model
        const solvers::multicomponentFilm& multicomponentFilm_;

        //- Reference to the transferred specie thermo package
        const fluidMulticomponentThermo& thermoFilm_;

        //- 1st coefficient of Diffusivity of specie
        const dimensionedScalar D1_;
        
        //- 2nd coefficient of Diffusivity of specie
        const dimensionedScalar D2_;

        // Dictionary to get value of sigma
        IOdictionary physicalProperties_;

        // Surface tension coefficient
        dimensionedScalar sigma_;
	
public:

    //- Runtime type information
    TypeName("higbieFilm");


    // Constructors

        //- Construct from components
        higbieFilm
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const label& patchID
        );


    //- Destructor
    virtual ~higbieFilm()
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
