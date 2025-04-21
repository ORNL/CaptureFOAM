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
    Foam::interphaseMassTransferModels::constant

Description
    Enforces a constant mass transfer rate across a bulk-film interface.
    The user provides the rate in units of kg/m^2/s. Mass transfer is taken 
    to be negative when a species is absorbed into the film from the bulk,
    and positive when a species is released from the film into the bulk.
    
    Note that the mass transfer rate is limited when the film thickness is
    below some minimum value, and also by the maximum possible rate calculated
    for each cell if it were to lose all mass in a single time step. This
    limiting calculation occurs in the interphaseMassTransferModel base class.

SourceFiles
    constant.C

\*---------------------------------------------------------------------------*/

#ifndef constant_H
#define constant_H

#include "interphaseMassTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interphaseMassTransferModels
{

/*---------------------------------------------------------------------------*\
                               Class constant
\*---------------------------------------------------------------------------*/

class constant
:
    public interphaseMassTransferModel
{
    // Private Data
        
        //- Constant interfacial mass transfer rate [kg/m^2/s]
        dimensionedScalar transferRate_;

public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from components
        constant
        (
            const dictionary& dict,
            const solvers::multicomponentFluid& fluid,
            const label& bulkPatchID,
            const mappedPatchBase& patchMap
        );


    //- Destructor
    virtual ~constant()
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
