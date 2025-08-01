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

#include "porousMedia.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace massTransferRateCoefficientModels
{
    defineTypeNameAndDebug(porousMedia, 0);
    addToRunTimeSelectionTable(massTransferRateCoefficientModel, porousMedia, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massTransferRateCoefficientModels::porousMedia::porousMedia
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
    C_(rateModelCoeffs_.lookupOrDefault<scalar>("C", 0.054)),
    D_(rateModelCoeffs_.lookup<scalar>("D")),
    dp_(rateModelCoeffs_.lookup<scalar>("poreDiameter")),
    eps_(rateModelCoeffs_.lookup<scalar>("porosity"))
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField&
Foam::massTransferRateCoefficientModels::porousMedia::k()
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
    
    scalar dum1 = 0.0;
    scalar dum2 = 0.0;
    scalar dum3 = 0.0;
    scalar dum4 = 0.0;
    scalar dum5 = 0.0;
    forAll(transferCells, facei)
    {
        const label celli = transferCells[facei];
        
        //- Calculate correlation quantities
        const scalar beta = Foam::pow(mu[facei] / (rho[facei] * D_), 0.333);
        const scalar gamma
            = Foam::pow(rho[facei] * mag(Ub[facei]) * dp_
                        / (mu[facei] * (1.0 - eps_)), 0.8);
        
        //- Set rate coefficient
        k_[celli] = C_ * beta * gamma * D_ / dp_;

	// Get the face index for the current cell
        const label faceIndex = mesh_.boundary()[patchID_].faceCells()[facei];

        // Calculate face area using the mesh object
        const scalar area = mag(mesh_.faceAreas()[faceIndex]);
        faceAreas[facei] = area;

	dum1 += mag(Ub[facei]) * area;
	dum2 += k_[celli] * area;
	dum3 += gamma * area;
	dum4 += beta * area;
	dum5 += area;
    }
    
    Info << "Surf averaged mag(Ub): " << dum1/dum5 << endl;
    Info << "Surf averaged gamma: " << dum3/dum5 << endl;
    Info << "Surf averaged beta: " << dum4/dum5 << endl;
    Info << "Surf averaged kg: " << dum2/dum5 << endl;

    return k_;
}

// ************************************************************************* //
