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

#include "interphaseMassTransferModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interphaseMassTransferModel, 0);
    defineRunTimeSelectionTable(interphaseMassTransferModel, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interphaseMassTransferModel::interphaseMassTransferModel
(
    const word& type,
    const dictionary& dict,
    const solvers::multicomponentFluid& fluid,
    const label& bulkPatchID,
    const mappedPatchBase& patchMap
)
:
    massTransferModelCoeffs_(dict.optionalSubDict(type + "Coeffs")),
    fluid_(fluid),
    patchMap_(patchMap),
    film_(patchMap_.nbrMesh().lookupObject<solvers::multicomponentFilm>(solver::typeName)),
    bulkMesh_(fluid_.p.mesh()),
    filmMesh_(film_.delta.mesh()),
    bulkPatchID_(bulkPatchID),
    filmPatchID_(film_.surfacePatch().index()),
    specieName_(dict.lookup("specie")),
    bulkSpecieID_(fluid_.thermo.species()[specieName_]),
    filmSpecieID_(film_.thermo.species()[specieName_]),
    mDotBulk_
    (
        IOobject
        (
            "mDotBulk",
            bulkMesh_.time().name(),
            bulkMesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        bulkMesh_,
        dimensionedScalar(dimMass/dimVolume/dimTime, 0.0)
    ),
    mDotFilm_
    (
        IOobject
        (
            "mDotFilm",
            filmMesh_.time().name(),
            filmMesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        filmMesh_,
        dimensionedScalar(dimMass/dimVolume/dimTime, 0.0)
    ),
    R_(film_.delta.size(), 0.0)
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::interphaseMassTransferModel::update()
{
    //- Calculate maximum possible mass transfer from bulk to film
    const volScalarField& rhob = fluid_.thermo.rho();
    const volScalarField& Yib = fluid_.thermo.Y(bulkSpecieID_);
    const volScalarField mDotMax = rhob * Yib / bulkMesh_.time().deltaT();
    const scalarField& Sb
        = bulkMesh_.magSf().boundaryField()[bulkPatchID_];
    
    //- Limit mass transfer in bulk phase
    const labelList& bulkCells
        = bulkMesh_.boundary()[bulkPatchID_].faceCells();

    forAll(bulkCells, facei)
    {
        const label celli = bulkCells[facei];
        
        const scalar Rmax
            = -1.0 * mDotMax[celli] * bulkMesh_.V()[celli] / Sb[facei];
        
        R_[facei] = max(R_[facei], Rmax);
    }
    
    //- Map interfacial mass transfer rate to film phase for limiting
    scalarField Rf(patchMap_.toNeighbour(R_));
    
    //- Calculate maximum possible mass transfer from film to bulk
    const volScalarField& rhof = film_.thermo.rho();
    const volScalarField& Yif = film_.thermo.Y(filmSpecieID_);
    const volScalarField mDotMin = rhof * Yif / filmMesh_.time().deltaT();
    const scalarField& Sf
        = filmMesh_.magSf().boundaryField()[filmPatchID_];
    
    //- Limit mass transfer in film phase
    const labelList& filmCells = film_.surfacePatch().faceCells();
        
    forAll(filmCells, facei)
    {
        const label celli = filmCells[facei];
        
        if (film_.delta[celli] > 1.0e-6)
        {
            const scalar Rmax
                = mDotMin[celli] * film_.delta.mesh().V()[celli] / Sf[facei];
                
            Rf[facei] = min(Rf[facei], Rmax);
        }
        else
        {
            Rf[facei] = 0.0;
        }
    }
    
    //- Map mass transfer back to bulk phase
    R_ = patchMap_.fromNeighbour(Rf);
    
    //- Update volumetric mass transfer fields using limited mass transfer rate
    mDotBulk_ = Zero;
    mDotFilm_ = Zero;
    
    forAll(bulkCells, facei)
    {
        const label celli = bulkCells[facei];
        
        mDotBulk_[celli] = R_[facei] * Sb[facei] / bulkMesh_.V()[celli];
    }
    
    forAll(filmCells, facei)
    {
        const label celli = filmCells[facei];
        
        mDotFilm_[celli]
            = -1.0 * Rf[facei] * Sf[facei] / filmMesh_.V()[celli];
    }
    
    dimensionedScalar mDotb = fvc::domainIntegrate(mDotBulk_);
    dimensionedScalar mDotf = fvc::domainIntegrate(mDotFilm_);
    
    Info << "Integrated bulk-side mass transfer: " << mDotb.value() << endl;
    Info << "Integrated film-side mass transfer: " << mDotf.value() << endl;
}

bool Foam::interphaseMassTransferModel::read()
{
    return true;
}

// ************************************************************************* //
