/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "multicomponentFilm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(multicomponentFilm, 0);
    addToRunTimeSelectionTable(solver, multicomponentFilm, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::multicomponentFilm::multicomponentFilm(fvMesh& mesh)
:
    isothermalFilm(mesh),
    
    thermo_(refCast<fluidMulticomponentThermo>(isothermalFilm::thermo_)),
    
    Y_(thermo_.Y()),

    reaction_(combustionModel::New(thermo_, momentumTransport())),

    thermophysicalTransport
    (
        filmMulticomponentThermophysicalTransportModel::New
        (
            momentumTransport(),
            thermo_
        )
    ),
    
    thermo(thermo_),
    Y(Y_),
    reaction(reaction_)
{
    thermo.validate(type(), "h", "e");
    
    forAll(Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(thermo.he());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::multicomponentFilm::~multicomponentFilm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::multicomponentFilm::prePredictor()
{
    isothermalFilm::prePredictor();

    if (pimple.predictTransport())
    {
        thermophysicalTransport->predict();
    }
}

void Foam::solvers::multicomponentFilm::postCorrector()
{
    isothermalFilm::postCorrector();

    if (pimple.correctTransport())
    {
        thermophysicalTransport->correct();
    }
}

// ************************************************************************* //
