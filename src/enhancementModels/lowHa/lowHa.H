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
    Foam::interphaseMassTransferModels::lowHa

Description
    Enhancement factor for low Hatta number (Ha <= 1.0) where E = Ha

SourceFiles
    lowHa.C

\*---------------------------------------------------------------------------*/

#ifndef lowHa_H
#define lowHa_H

#include "enhancementModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace enhancementModels
{

/*---------------------------------------------------------------------------*\
                                  Class lowHa
\*---------------------------------------------------------------------------*/

class lowHa
:
    public enhancementModel
{
    // Private Data

        //- 1st coefficient of Diffusivity of specie in liquid
        const dimensionedScalar D1_;

	//- 2nd coefficient of Diffusivity of specie in liquid
        const dimensionedScalar D2_;

        //- Time to turn enhancement on
        const scalar tStart_;

public:

    //- Runtime type information
    TypeName("lowHa");


    // Constructors

        //- Construct from components
        lowHa
        (
            const dictionary& dict,
            const solvers::multicomponentFilm& film,
	    const solvers::multicomponentFluid& fluid,
	    const label& filmSpecieID
        );


    //- Destructor
    virtual ~lowHa()
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
