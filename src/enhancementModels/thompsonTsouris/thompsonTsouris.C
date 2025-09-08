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

#include "thompsonTsouris.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include <iostream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace enhancementModels
{
    defineTypeNameAndDebug(thompsonTsouris, 0);
    addToRunTimeSelectionTable(enhancementModel, thompsonTsouris, dictionary);

    Foam::volScalarField& thompsonTsouris::Einf()
    {
	dimensionedScalar eps(dimMoles/dimVolume, 1.0e-16);
	dimensionedScalar eps2(dimVolume/dimMoles, 1.0e-16);

    	const volScalarField& rho = filmMesh_.lookupObject<volScalarField>("rho");
	const volScalarField& CO2 = filmMesh_.lookupObject<volScalarField>("CO2");
	const dimensionedScalar Wco2(dimMass/dimMoles, 44.0);
	const volScalarField Cco2 = (rho * CO2 / Wco2) + eps;

	const volScalarField& rhobulk = bulkMesh_.lookupObject<volScalarField>("rho");
	const volScalarField& CO2bulk = bulkMesh_.lookupObject<volScalarField>("CO2");
	const volScalarField Cco2b = (rhobulk * CO2bulk / Wco2) + eps;

	if (filmMesh_.foundObject<volScalarField>("MEA"))
	{
	  const volScalarField& MEA = filmMesh_.lookupObject<volScalarField>("MEA");
	  const volScalarField& MEAp = filmMesh_.lookupObject<volScalarField>("MEA+");
	  const volScalarField& MEACOO = filmMesh_.lookupObject<volScalarField>("MEACOO-");

	  const dimensionedScalar Wmea(dimMass/dimMoles, 61.0);
	  const dimensionedScalar Wmeap(dimMass/dimMoles, 62.0);
	  const dimensionedScalar Wmeacoo(dimMass/dimMoles, 104.0);

	  const volScalarField C1 = rho * MEA / Wmea;
	  const volScalarField C2 = rho * MEAp / Wmeap;
	  const volScalarField C3 = rho * MEACOO / Wmeacoo;

	  const volScalarField Keq = (C2 * C3 / (Cco2 * Foam::pow(C1, 2.0))) + eps2;

	  dimensionedScalar Dmea(dimArea/dimTime, 1.61e-9);
	  dimensionedScalar Dmeacoo(dimArea/dimTime, 1.5e-9);	  

	  volScalarField num = Foam::sqrt(Keq) * C1 * (Dmeacoo / D_);
	  volScalarField denom = (1.0 + (2.0 * Dmeacoo / Dmea)) * Foam::sqrt(Keq * Cco2) * (Foam::sqrt(Cco2) + Foam::sqrt(Cco2));

	  Einf_ = 1.0 + (num / denom);
	}

	else if (filmMesh_.foundObject<volScalarField>("KSAR"))
	{
	  const volScalarField& KSAR = filmMesh_.lookupObject<volScalarField>("KSAR");
	  const volScalarField& KSARp = filmMesh_.lookupObject<volScalarField>("KSAR+");
	  const volScalarField& KSARCOO = filmMesh_.lookupObject<volScalarField>("KSARCOO-");

	  const dimensionedScalar Wksar(dimMass/dimMoles, 127.0);
	  const dimensionedScalar Wksarp(dimMass/dimMoles, 128.0);
	  const dimensionedScalar Wksarcoo(dimMass/dimMoles, 171.0);

	  const volScalarField C1 = rho * KSAR / Wksar;
	  const volScalarField C2 = rho * KSARp / Wksarp;
	  const volScalarField C3 = rho * KSARCOO / Wksarcoo;

	  const volScalarField Keq = C2 * C3/ (Cco2 * Foam::pow(C1, 2.0)) + eps2;

	  dimensionedScalar Dksar(dimArea/dimTime, 1.0e-9);
          dimensionedScalar Dksarcoo(dimArea/dimTime, 1.0e-9);

	  volScalarField num = Foam::sqrt(Keq) * C1 * (Dksarcoo / D_);
	  volScalarField denom = (1.0 + (2.0 * Dksarcoo / Dksar)) * Foam::sqrt(Keq * Cco2) * (Foam::sqrt(Cco2) + Foam::sqrt(Cco2));

	  Einf_ = 1.0 + (num / denom);
	}

	else
    	{
          //- Throw fatal error for unknown chemistry
          FatalErrorInFunction
          << "Neither MEA nor KSAR found in mixture." << nl
          << "Please add known solvent for enhancement model." << endl
          << exit(FatalError);
    	}

	return Einf_;
    }
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::enhancementModels::thompsonTsouris::thompsonTsouris
(
    const dictionary& dict,
    const solvers::multicomponentFilm& film,
    const solvers::multicomponentFluid& fluid,
    const label& filmSpecieID
)
:
    enhancementModel
    (
        typeName,
        dict,
        film,
        fluid,
        filmSpecieID
    ),

    D1_(dimArea/dimTime/dimTemperature, massTransferModelCoeffs_.lookup<scalar>("Dl1")),
    D2_(dimArea/dimTime, massTransferModelCoeffs_.lookup<scalar>("Dl2")),

    D_
    (
	IOobject
        (
            "D",
            filmMesh_.time().name(),
            filmMesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        filmMesh_,
        dimensionedScalar(dimArea/dimTime, 0.0)
    ),

    tStart_(massTransferModelCoeffs_.lookupOrDefault<scalar>("tStart", 0.0)),

    Einf_
    (
        IOobject
	(
	    "Einf",
	    filmMesh_.time().name(),
	    filmMesh_,
	    IOobject::READ_IF_PRESENT,
	    IOobject::AUTO_WRITE
	),
	filmMesh_,
	dimensionedScalar(dimless, 0.01)
    )

{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::enhancementModels::thompsonTsouris::update()
{

    //- Look up film-side mass transfer rate coefficient field
    const volScalarField& k_l = filmMesh_.lookupObject<volScalarField>("k");
    const volScalarField& Tf = filmMesh_.lookupObject<volScalarField>("T");

    const volScalarField klLim
        = max(k_l, dimensionedScalar(dimVelocity, 1e-8));

    dimensionedScalar D1(dimArea/dimTime/dimTemperature, D1_.value());
    dimensionedScalar D2(dimArea/dimTime, D2_.value());
   // const volScalarField D = (D1 * Tf) + D2;
    D_ = (D1 * Tf) + D2;

    //- Set E
    if (filmMesh_.time().value() >= tStart_)
    {
	Einf();

	const volScalarField Ha = Foam::sqrt(D_ * enhancementModel::kApp()) / klLim;

	volScalarField denom = (1.0 / Foam::pow(Ha - 1.0, 1.35))
		+ (1.0 / Foam::pow(Einf_ - 1.0, 1.35));

        E_ = 1.0 + (1.0 / Foam::pow(denom, 0.74));
  
//	E_ = Einf_;
//	Info << "E_: " << E_ << " at time: " << filmMesh_.time().value() << endl;
    }
}

bool Foam::enhancementModels::thompsonTsouris::read()
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
