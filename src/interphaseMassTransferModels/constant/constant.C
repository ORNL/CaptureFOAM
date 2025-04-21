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

\*---------------------------------------------------------------------------*/

#include "constant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interphaseMassTransferModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(interphaseMassTransferModel, constant, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interphaseMassTransferModels::constant::constant
(
    const dictionary& dict,
    const solvers::multicomponentFluid& fluid,
    const label& bulkPatchID,
    const mappedPatchBase& patchMap
)
:
    interphaseMassTransferModel
    (
        typeName,
        dict,
        fluid,
        bulkPatchID,
        patchMap
    ),
    transferRate_
    (
        "transferRate",
        dimMass/dimArea/dimTime,
        massTransferModelCoeffs_.lookup<scalar>("transferRate")
    )
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interphaseMassTransferModels::constant::update()
{    
    //- Set interfacial mass transfer rate to be equal to the user-defined
    //  constant rate
    tmp<scalarField> R(new scalarField(R_.size(), transferRate_.value()));
    
    R_ = R;
    
    //- Call base class update to perform limiting based on maximum possible 
    //  mass transfer and minimum film thickness, and update the volume fields
    interphaseMassTransferModel::update();
}


bool Foam::interphaseMassTransferModels::constant::read()
{
    if (interphaseMassTransferModel::read())
    {
        return true;
    }
    
    else
    {
        return false;
    }
}

// ************************************************************************* //
