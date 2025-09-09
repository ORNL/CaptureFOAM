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

#include "blasiusGas.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferRateCoefficientModels
{
    defineTypeNameAndDebug(blasiusGas, 0);
    addToRunTimeSelectionTable(massTransferRateCoefficientModel, blasiusGas, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferRateCoefficientModels::blasiusGas::blasiusGas
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
    
    multicomponentFluid_(mesh.lookupObject<solvers::multicomponentFluid>(solver::typeName)),
    thermo_(multicomponentFluid_.thermo),
    C_(rateModelCoeffs_.lookupOrDefault<scalar>("C", 0.332)),
    D_(rateModelCoeffs_.lookup<scalar>("D")),
    inlet_(rateModelCoeffs_.lookup<vector>("inlet")),
    dir_(rateModelCoeffs_.lookup<vector>("direction"))
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField&
Foam::massTransferRateCoefficientModels::blasiusGas::k()
{
    //- Reset rate coefficient to zero
    k_ = Zero;
    
    //- Look up velocity field from mesh
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");
    const vectorField& Ub = U.boundaryField()[patchID_];
    
    //- Get density and porosity fields from solver
    const volScalarField& rho = thermo_.rho();
    const volScalarField& mu = thermo_.mu();
    
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
	// Get the face index for the current cell
        const label faceIndex = mesh_.boundary()[patchID_].faceCells()[facei];

        // Calculate face area using the mesh object
        const scalar area = mag(mesh_.faceAreas()[faceIndex]);
        faceAreas[facei] = area;
        
        //- Calculate correlation quantities
	const scalar facePosition = mag((mesh_.C()[celli] - inlet_) & dir_);
	const scalar faceRe = rho[facei] * mag(Ub[facei]) * facePosition / mu[facei];
	const scalar faceSc = mu[facei] / (rho[facei] * D_);
        
        //- Set rate coefficient
	k_[celli] = C_ * (D_ / facePosition) * Foam::pow(faceRe, 0.5) 
		* Foam::pow(faceSc, 0.333);

	patchVelocity += mag(Ub[facei]) * area;
	patchMassTransferCoefficient += k_[celli] * area;
	patchArea += area;
    }
    
    Info << "Surf averaged mag(Ub): " << patchVelocity / patchArea << endl;
    Info << "Surf averaged kg: " << patchMassTransferCoefficient / patchArea << endl;

    return k_;
}

// ************************************************************************* //
