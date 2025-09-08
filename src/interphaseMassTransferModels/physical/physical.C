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

#include "physical.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interphaseMassTransferModels
{
    defineTypeNameAndDebug(physical, 0);
    addToRunTimeSelectionTable(interphaseMassTransferModel, physical, dictionary);
}
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interphaseMassTransferModels::physical::physical
(
    const dictionary& dict,
    const solvers::multicomponentFluid& fluid,
    const label& bulkPatchID,
    const mappedPatchBase& patchMap
)
:
    interphaseMassTransferModel
    (
        typeName,
        dict,
        fluid,
        bulkPatchID,
        patchMap
    ),
    H1_(dimPressure*dimVolume/dimTemperature/dimMoles, massTransferModelCoeffs_.lookup<scalar>("H1")),
    H2_(dimPressure*dimVolume/dimMoles, massTransferModelCoeffs_.lookup<scalar>("H2")),
    filmRateModel_(nullptr),
    bulkRateModel_(nullptr),
    enhancementModel_(nullptr),
    N_
    (
        IOobject
        (
            "N",
            bulkMesh_.time().name(),
            bulkMesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        bulkMesh_,
        dimensionedScalar(dimMoles/dimArea/dimTime, 0.0)
    )
{
    //- Create mass transfer rate coefficient model in film
    filmRateModel_ =
        massTransferRateCoefficientModel::New
        (
            massTransferModelCoeffs_.subDict("filmMassTransferRateModel"),
            filmMesh_,
            filmPatchID_
        );
        
    //- Create mass transfer rate coefficient model in bulk
    bulkRateModel_ =
        massTransferRateCoefficientModel::New
        (
            massTransferModelCoeffs_.subDict("bulkMassTransferRateModel"),
            bulkMesh_,
            bulkPatchID_
        );
        
    //- Create enhancement model in film
    enhancementModel_ =
        enhancementModel::New
        (
            massTransferModelCoeffs_,
            film_,
            fluid_,
            filmSpecieID_
        );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interphaseMassTransferModels::physical::update()
{    
    //- Calculate the bulk-side mass transfer rate
    tmp<scalarField> Rbt(new scalarField(R_.size(), 0.0));
    scalarField& Rb = Rbt.ref();
        
    const volScalarField& Tb = bulkMesh_.lookupObject<volScalarField>("T");

    dimensionedScalar Ru(dimEnergy/dimTemperature/dimMoles, 8314.0);
    
    const volScalarField& kg = bulkRateModel_->k();
    
    const labelList& bulkCells
        = bulkMesh_.boundary()[bulkPatchID_].faceCells();
        
    forAll(bulkCells, facei)
    {
        const label celli = bulkCells[facei];
        
        Rb[facei] = kg[celli] / Ru.value() / max(Tb[celli], 1e-12);
    }

    //- Calculate the film-side mass transfer rate
    tmp<scalarField> Rft(new scalarField(R_.size(), 0.0));
    scalarField& Rf = Rft.ref();

    const volScalarField& kl = filmRateModel_->k();
    
    //- Temperature dependent Henry constant
    tmp<scalarField> Ht(new scalarField(R_.size(), 0.0));
    scalarField& H = Ht.ref();

    const volScalarField& Tf = filmMesh_.lookupObject<volScalarField>("T");

    dimensionedScalar H1(dimPressure*dimVolume/dimTemperature/dimMoles, H1_.value());
    dimensionedScalar H2(dimPressure*dimVolume/dimMoles, H2_.value());

    H = (H1 * Tf) + H2;

    //- Update enhancement factor
    enhancementModel_->update();
    const volScalarField& E = enhancementModel_->E();
    
    const labelList& filmCells = film_.surfacePatch().faceCells();
    
    //- Calculate film-side rate contribution
    forAll(filmCells, facei)
    {
        const label celli = filmCells[facei];

        Rf[facei] = max(E[celli] * kl[celli] / H[celli], 1e-20);
    }
    
    //- Overall mass transfer rate coeff is the harmonic mean of bulk-
    //  and film-side contributions
    const scalarField K = 1.0 / (1.0 / Rb + 1.0 / patchMap_.fromNeighbour(Rf));
    
    //- Calculate equilibrium partial pressure in gas from film concentration
    tmp<scalarField> Pstart(new scalarField(R_.size(), 0.0));
    scalarField& Pstar = Pstart.ref();
    
    const scalarField& Yif = film_.thermo.Y(filmSpecieID_);
    const volScalarField& rhof = film_.thermo.rho();
    const scalarField& rhofb = rhof.boundaryField()[filmPatchID_];
    const dimensionedScalar& W = film_.thermo.Wi(filmSpecieID_);
    const scalarField Cif = Yif * rhofb / W.value();
    
    Pstar = H * max(Cif, 1e-12);
    
    //- Get pointers to fields for bulk phase partial pressure
    const scalarField& Pb
        = bulkMesh_.lookupObject<volScalarField>("p").boundaryField()[bulkPatchID_];
    const scalarField& Yib
        = fluid_.thermo.Y(bulkSpecieID_).boundaryField()[bulkPatchID_];
    
    //- Calculate overall mass transfer rate
    R_ = -1.0 * W.value() * K * (Pb * Yib - patchMap_.fromNeighbour(Pstar));
    
    //- Call base class update to perform limiting based on maximum possible 
    //  mass transfer and minimum film thickness, and update the volume fields
    interphaseMassTransferModel::update();
    
    //- Update mass transfer debug field
    N_ = Zero;
    
    forAll(bulkCells, facei)
    {
        N_[bulkCells[facei]] = R_[facei];
    }
}


bool Foam::interphaseMassTransferModels::physical::read()
{
    if (interphaseMassTransferModel::read())
    {
        return true;
    }
    
    else
    {
        return false;
    }
}

// ************************************************************************* //
