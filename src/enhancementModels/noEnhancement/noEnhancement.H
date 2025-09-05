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
    Foam::interphaseMassTransferModels::noEnhancement

Description
    Returns unity-valued enhancement factor field.

SourceFiles
    noEnhancement.C

\*---------------------------------------------------------------------------*/

#ifndef noEnhancement_H
#define noEnhancement_H

#include "enhancementModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace enhancementModels
{

/*---------------------------------------------------------------------------*\
                               Class noEnhancement
\*---------------------------------------------------------------------------*/

class noEnhancement
:
    public enhancementModel
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from components
        noEnhancement
        (
            const dictionary& dict,
            const solvers::multicomponentFilm& film,
	    const solvers::multicomponentFluid& fluid,
            const label& filmSpecieID
        );


    //- Destructor
    virtual ~noEnhancement()
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
