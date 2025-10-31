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

#include "Dobbins.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferRateCoefficientModels
{
    defineTypeNameAndDebug(Dobbins, 0);
    addToRunTimeSelectionTable(massTransferRateCoefficientModel, Dobbins, 
		    dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferRateCoefficientModels::Dobbins::Dobbins
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
Foam::massTransferRateCoefficientModels::Dobbins::k()
{
    //- Reset rate coefficient to zero
    k_ = Zero;

    //- Look up velocity field from mesh
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    const vectorField& Uboun = U.boundaryField()[patchID_];

    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    const scalarField& Tboun = T.boundaryField()[patchID_];

    //- Get density field from solver
    const volScalarField& rho = thermo_->rho();
    
    //- Loop over cells adjacent to transfer patch
    const labelList& transferCells = mesh_.boundary()[patchID_].faceCells();
    
    forAll(transferCells, facei)
    {
        const label celli = transferCells[facei];
        const scalar D = (D1_.value() * Tboun[facei]) + D2_.value();

        //- Calculate surface renewal time and coefficient
	// kl = sqrt(D*s)*coth(sqrt(s * delta^2 / D))
        if (isFilm_)
        {
           const volScalarField& delta = 
		   mesh_.lookupObject<volScalarField>("delta");

           const scalar Ui = max(mag(Uboun[facei]), 1e-12);
           const scalar s = Ui / delta[celli];
           const scalar arg = Foam::sqrt(s * delta[celli] * delta[celli] / D);

           //- Set rate coefficient
           k_[celli] = Foam::sqrt(D * s) * Foam::cosh(arg) / Foam::sinh(arg);
        }
        else
        {
           FatalErrorInFunction
                << "Model not available for the bulk phase"
                << abort(FatalError);   
        }
    }

    return k_;
}

// ************************************************************************* //
