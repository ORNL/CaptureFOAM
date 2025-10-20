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

#include "pipeBlasius.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferRateCoefficientModels
{
    defineTypeNameAndDebug(pipeBlasius, 0);
    addToRunTimeSelectionTable(massTransferRateCoefficientModel, pipeBlasius, 
		    dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferRateCoefficientModels::pipeBlasius::pipeBlasius
(
    const dictionary& dict,
    const fvMesh& mesh,
    const label& patchID
)
:
    massTransferRateCoefficientModel
    (
        typeName,
        dict,
        mesh,
        patchID
    ),
    
    D1_(dimArea/dimTime/dimTemperature, rateModelCoeffs_.lookup<scalar>("D1")),
    D2_(dimArea/dimTime, rateModelCoeffs_.lookup<scalar>("D2")),
    diam_(rateModelCoeffs_.lookup<scalar>("Pipe Diameter"))
{

        if (mesh.foundObject<solvers::multicomponentFilm>("solver"))
        {
           const auto& solver = 
		   mesh.lookupObject<solvers::multicomponentFilm>("solver");
           thermo_ = &solver.thermo;
	   isFilm_ = true;
        }
        else if (mesh.foundObject<solvers::multicomponentFluid>("solver"))
        {
           const auto& solver = 
		   mesh.lookupObject<solvers::multicomponentFluid>("solver");
           thermo_ = &solver.thermo;
	   isFilm_ = false;
        }
        else
        {
           FatalErrorInFunction
                << "Cannot find either multicomponentFilm or"
                << " multicomponentFluid object in the mesh registry"
                << abort(FatalError);
        }

}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField&
Foam::massTransferRateCoefficientModels::pipeBlasius::k()
{
    //- Reset rate coefficient to zero
    k_ = Zero;

    //- Look up fields from mesh
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    const vectorField& Uboun = U.boundaryField()[patchID_];

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    const scalarField& Tboun = T.boundaryField()[patchID_];

    //- Loop over cells adjacent to transfer patch
    const labelList& transferCells = mesh_.boundary()[patchID_].faceCells();

    // Create an array to hold the area of the faces
    scalarField faceAreas(transferCells.size());

    scalar patchVelocity = 0.0;
    scalar patchMassTransferCoefficient = 0.0;
    scalar patchArea = 0.0;

    forAll(transferCells, facei)
    {
        const label celli = transferCells[facei];
        const scalar D = (D1_.value() * Tboun[facei]) + D2_.value();
        
        //- Set rate coefficient
        k_[celli] = 3.66 * (D / diam_);

        if (!isFilm_)
        {
           const label faceIndex = 
                               mesh_.boundary()[patchID_].faceCells()[facei];
           // Calculate face area using the mesh object
           const scalar area = mag(mesh_.faceAreas()[faceIndex]);
           faceAreas[facei] = area;

           patchVelocity += mag(Uboun[facei]) * area;
           patchMassTransferCoefficient += k_[celli] * area;
           patchArea += area;
        }
    }

    if (!isFilm_)
    {
        Info << "Surf averaged mag(Ub): " << patchVelocity / patchArea << endl;
        Info << "Surf averaged kg: " << patchMassTransferCoefficient 
                               / patchArea << endl;
    }

    return k_;
}

// ************************************************************************* //
