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

#include "enhancementModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(enhancementModel, 0);
    defineRunTimeSelectionTable(enhancementModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::enhancementModel::enhancementModel
(
    const word& type,
    const dictionary& dict,
    const solvers::multicomponentFilm& film,
    const solvers::multicomponentFluid& fluid,
    const label& filmSpecieID
)
:
    massTransferModelCoeffs_(dict.optionalSubDict(type + "Coeffs")),
    film_(film),
    filmMesh_(film_.delta.mesh()),
    filmSpecieID_(filmSpecieID),
    E_
    (
        IOobject
        (
            "E",
            filmMesh_.time().name(),
            filmMesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        filmMesh_,
        dimensionedScalar(dimless, 1.0)
    ),
    kApp_
    (
        IOobject
        (
            "kApp",
            filmMesh_.time().name(),
            filmMesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        filmMesh_,
        dimensionedScalar(dimless / dimTime, 1.0)
    ),
    fluid_(fluid),
    bulkMesh_(fluid_.p.mesh())

{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::enhancementModel::kApp()
{
    if (filmMesh_.foundObject<volScalarField>("MEA")
		      && filmMesh_.foundObject<volScalarField>("KSAR"))
    {
	//- MEA and KSAR cannot be present as of now
	//  (don't have mechanism for combination)    
	FatalErrorInFunction
	<< "Both MEA and KSAR found in mixture." << nl
	<< "Mixed reaction model is not implemented." << nl
	<< "Please use single solvent chemistry." << endl
	<< exit(FatalError);
    }
		
    else if (filmMesh_.foundObject<volScalarField>("MEA"))
    {
        //- Calculate herereaction rate coeffs
        const volScalarField& T = filmMesh_.lookupObject<volScalarField>("T");
        
        const dimensionedScalar dimK
        (
            Foam::pow(dimLength, 6) / Foam::pow(dimMoles, 2) / dimTime,
            1.0
        );
        
        const volScalarField k1
            = dimK * 4.61e+9
              * Foam::exp(dimensionedScalar(dimTemperature, -4412)/T);
              
        const volScalarField k2
            = dimK * 4.55e+6
              * Foam::exp(dimensionedScalar(dimTemperature, -3287)/T);
        
        //- Calculate concentration of species
        const volScalarField& rho
            = filmMesh_.lookupObject<volScalarField>("rho");
        const volScalarField& MEA
            = filmMesh_.lookupObject<volScalarField>("MEA");
        const volScalarField& H2O
            = filmMesh_.lookupObject<volScalarField>("H2O");
        
        const dimensionedScalar Wmea(dimMass/dimMoles, 61.0);
        const dimensionedScalar Wh2o(dimMass/dimMoles, 18.0);
        
        const volScalarField Cmea = rho * MEA / Wmea;
        const volScalarField Ch2o = rho * H2O / Wh2o;

        //- Update apparent reaction rate
        kApp_ = k1 * Foam::pow(Cmea, 2.0) + k2 * Cmea * Ch2o;

	//- Calculate instantaneous enhancement factor
//	const volScalarField& MEAp = filmMesh_.lookupObject<volScalarField>("MEA+");
//	const volScalarField& MEACOOm = filmMesh_.lookupObject<volScalarField>("MEACOO-");
//	const volScalarField& CO2 = filmMesh_.lookupObject<volScalarField>("CO2");

//	const dimensionedScalar Wmeap(dimMass/dimMoles, 62.0);
//	const dimensionedScalar Wmeacoo(dimMass/dimMoles, 104.0);
//	const dimensionedScalar Wco2(dimMass/dimMoles, 44.0);
	
//	const volScalarField Cmeap = rho * MEAp / Wmeap;
//	const volScalarField Cmeacoo = rho * MEACOOm / Wmeacoo;
//	const volScalarField Cco2 = rho * CO2 / Wco2;

//	Keq_ = Cmeap * Cmeacoo/(Cco2 * Foam::pow(Cmea,2.));
    }

    else if (filmMesh_.foundObject<volScalarField>("KSAR"))
    {	    
                      
        //- Calculate concentration of species
        const volScalarField& rho
            = filmMesh_.lookupObject<volScalarField>("rho");
        const volScalarField& KSAR
            = filmMesh_.lookupObject<volScalarField>("KSAR");
        const volScalarField& H2O
            = filmMesh_.lookupObject<volScalarField>("H2O");
        
        const dimensionedScalar Wksar(dimMass/dimMoles, 127.0);
        const dimensionedScalar Wh2o(dimMass/dimMoles, 18.0);
        
        const volScalarField Cksar = rho * KSAR / Wksar;
        const volScalarField Ch2o = rho * H2O / Wh2o;

        //- Calculate reaction rate coeffs
        const volScalarField& T = filmMesh_.lookupObject<volScalarField>("T");

	const dimensionedScalar dimK
	(
	    Foam::pow(dimLength, 6) / Foam::pow(dimMoles, 2) / dimTime,
	    1.0
	);

	const volScalarField k1
           = dimK * 6.35e+6
	     * Foam::exp(dimensionedScalar(dimTemperature, -1590)/T);

        const volScalarField k2
 	   = dimK * 3.98e+8 
	     * Foam::exp(dimensionedScalar(dimTemperature, -3924)/T);
   
        const volScalarField Cksar_ion = Foam::exp(
		dimensionedScalar(dimVolume/dimMoles, 0.38) * Cksar);

        //- Update apparent reaction rate
        kApp_ = (k1 * Foam::pow(Cksar, 2.0) + k2 * Cksar * Ch2o) * Cksar_ion;
        
        //- Add hydroxide contribution if present
        if (filmMesh_.foundObject<volScalarField>("OH"))
        {
            //- Calculate hydroxide concentration
            const volScalarField& OH
                = filmMesh_.lookupObject<volScalarField>("OH");
            const dimensionedScalar Woh(dimMass/dimMoles, 17.0);
            const volScalarField Coh = rho * OH / Woh;

            //- Calculate rate coefficient for hydroxide reaction
            const volScalarField kOHinf
                = 3.28e+13
                  * Foam::exp(dimensionedScalar(dimTemperature, -6612) / T);
            
            //- Temperature dependence valid from 20*C through 50*C
            const volScalarField betaCO2
                = dimensionedScalar(Foam::pow(dimless/dimTemperature, 2.0), 3.40e-4)
                  * Foam::pow(T, 2.0)
                  - dimensionedScalar(dimless/dimTemperature, 2.12e-1) * T
                  + 33.51;
                  
            const volScalarField kOH = kOHinf * Foam::exp(betaCO2 * Coh);
            
            kApp_ += (kOH * Coh);
        }
    }

    else
    {
        //- Throw fatal error for unknown chemistry
        FatalErrorInFunction
        << "Neither MEA nor KSAR found in mixture." << nl
        << "Please add known solvent for enhancement model." << endl
        << exit(FatalError);
    }
    
    return kApp_;
}

bool Foam::enhancementModel::read()
{
    return true;
}

// ************************************************************************* //
