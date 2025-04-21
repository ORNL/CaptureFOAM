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

#include "lowHa.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace enhancementModels
{
    defineTypeNameAndDebug(lowHa, 0);
    addToRunTimeSelectionTable(enhancementModel, lowHa, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::enhancementModels::lowHa::lowHa
(
    const dictionary& dict,
    const solvers::multicomponentFilm& film,
    const label& filmSpecieID
)
:
    enhancementModel
    (
        typeName,
        dict,
        film,
        filmSpecieID
    ),
    
    D1_(dimArea/dimTime/dimTemperature, massTransferModelCoeffs_.lookup<scalar>("Dl1")),
    D2_(dimArea/dimTime, massTransferModelCoeffs_.lookup<scalar>("Dl2")),
    tStart_(massTransferModelCoeffs_.lookupOrDefault<scalar>("tStart", 0.0))
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::enhancementModels::lowHa::update()
{
    //- Look up film-side mass transfer rate coefficient field
    const volScalarField& k_l = filmMesh_.lookupObject<volScalarField>("k");
    const volScalarField& Tf = filmMesh_.lookupObject<volScalarField>("T");

    const volScalarField klLim
        = max(k_l, dimensionedScalar(dimVelocity, 1e-8));

    dimensionedScalar D1(dimArea/dimTime/dimTemperature, D1_.value());
    dimensionedScalar D2(dimArea/dimTime, D2_.value());
    const volScalarField D = (D1 * Tf) + D2;

    //- Set E = Ha
    if (filmMesh_.time().value() >= tStart_)
    {
        E_ = Foam::sqrt(D * enhancementModel::kApp()) / klLim;
    }
}


bool Foam::enhancementModels::lowHa::read()
{
    if (enhancementModel::read())
    {
        return true;
    }
    
    else
    {
        return false;
    }
}

// ************************************************************************* //
