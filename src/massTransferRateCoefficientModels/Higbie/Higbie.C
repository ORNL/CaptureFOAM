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

#include "Higbie.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferRateCoefficientModels
{
    defineTypeNameAndDebug(Higbie, 0);
    addToRunTimeSelectionTable(massTransferRateCoefficientModel, Higbie, 
		    dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferRateCoefficientModels::Higbie::Higbie
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
    physicalProperties_
    (
        IOobject
        (
            "physicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    sigma_("sigma", dimForce/dimLength, 0.0)

{

	if (mesh.foundObject<solvers::multicomponentFilm>("solver"))
        {
           const auto& solver = 
		   mesh.lookupObject<solvers::multicomponentFilm>("solver");
           thermo_ = &solver.thermo;

	   if (physicalProperties_.found("sigma"))
           {
              sigma_ = dimensionedScalar("sigma", dimForce/dimLength, 
			      physicalProperties_.subDict("sigma"));
           }
           else
           {
	      FatalErrorInFunction
                << "sigma not found in Film physicalProperties"
                << abort(FatalError);
           }
	   
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
Foam::massTransferRateCoefficientModels::Higbie::k()
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
    
    // Create an array to hold the area of the faces
    scalarField faceAreas(transferCells.size());

    scalar patchVelocity = 0.0;
    scalar patchMassTransferCoefficient = 0.0;
    scalar patchArea = 0.0;

    forAll(transferCells, facei)
    {
        const label celli = transferCells[facei];
        const scalar D = (D1_.value() * Tboun[facei]) + D2_.value();

        //- Calculate contact time based on smallest eddies
        if (isFilm_)
        {
           const volScalarField& delta = 
		   mesh_.lookupObject<volScalarField>("delta");

           const scalar L = delta[celli];
           const scalar Ui = max(mag(Uboun[facei]), 1e-12);
           const scalar tau_conv = max(L / Ui, 1e-12);
           const scalar tau_cap = Foam::pow(rho[facei] * Foam::pow(L, 3.0) 
			/ sigma_.value(), 0.5);
           const scalar tau = min(tau_conv, tau_cap);

           //- Set rate coefficient
           k_[celli] = 2.0 * Foam::sqrt(D / pi / tau);
        }
        else
        {
           const label faceIndex = 
		   mesh_.boundary()[patchID_].faceCells()[facei];
           // Calculate face area using the mesh object
           const scalar area = mag(mesh_.faceAreas()[faceIndex]);
           faceAreas[facei] = area;

           const scalar L = Foam::pow(area, 0.5);
           const scalar Ui = mag(Uboun[facei]);
           const scalar tau = max(L / Ui, 1e-12);

           //- Set rate coefficient
           k_[celli] = 2.0 * Foam::sqrt(D / pi / tau);

           patchVelocity += Ui * area;
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
