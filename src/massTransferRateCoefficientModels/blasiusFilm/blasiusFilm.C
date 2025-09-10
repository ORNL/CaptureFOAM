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

#include "blasiusFilm.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferRateCoefficientModels
{
    defineTypeNameAndDebug(blasiusFilm, 0);
    addToRunTimeSelectionTable(massTransferRateCoefficientModel, blasiusFilm, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferRateCoefficientModels::blasiusFilm::blasiusFilm
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
    
    multicomponentFilm_(mesh.lookupObject<solvers::multicomponentFilm>(solver::typeName)),
    thermoFilm_(multicomponentFilm_.thermo),
    D1_(dimArea/dimTime/dimTemperature, rateModelCoeffs_.lookup<scalar>("D1")),
    D2_(dimArea/dimTime, rateModelCoeffs_.lookup<scalar>("D2")),
    C_(rateModelCoeffs_.lookupOrDefault<scalar>("C", 0.332)),
    inlet_(rateModelCoeffs_.lookup<vector>("inlet")),
    dir_(rateModelCoeffs_.lookup<vector>("direction"))
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField&
Foam::massTransferRateCoefficientModels::blasiusFilm::k()
{
    //- Reset rate coefficient to zero
    k_ = Zero;

    //- Look up fields from mesh
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    const vectorField& Uboun = U.boundaryField()[patchID_];

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    const scalarField& Tboun = T.boundaryField()[patchID_];

    //- Get density and viscosity fields from solver
    const volScalarField& rho = thermoFilm_.rho();
    const volScalarField& mu = thermoFilm_.mu();

    //- Loop over cells adjacent to transfer patch
    const labelList& transferCells = mesh_.boundary()[patchID_].faceCells();
    
    forAll(transferCells, facei)
    {
        const label celli = transferCells[facei];
        
	//- Calculate correlation coefficients
	const scalar D = (D1_.value() * Tboun[facei]) + D2_.value();
	const scalar facePosition = mag((mesh_.C()[celli] - inlet_) & dir_);
        const scalar faceRe = rho[facei] * mag(Uboun[facei]) * facePosition / mu[facei];
	const scalar faceSc = mu[facei] / (rho[facei] * D);

        //- Set rate coefficient
        k_[celli] = C_ * (D / facePosition) * Foam::pow(faceRe, 0.5)
                * Foam::pow(faceSc, 0.333);
    }

    return k_;
}

// ************************************************************************* //
