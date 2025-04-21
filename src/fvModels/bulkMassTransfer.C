/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "bulkMassTransfer.H"
#include "filmMassTransfer.H"
#include "mappedPatchBase.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(bulkMassTransfer, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            bulkMassTransfer,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::bulkMassTransfer::bulkMassTransfer
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    
    multicomponentFluid_(mesh.lookupObject<solvers::multicomponentFluid>(solver::typeName)),
    thermo_(multicomponentFluid_.thermo),
    
    specieName_(dict.lookup("specie")),
    
    filmPatchName_(dict.lookup("transferPatch")),
    filmPatchi_(mesh.boundaryMesh().findIndex(filmPatchName_)),
    filmPatch_(mesh.boundaryMesh()[filmPatchi_]),
    filmPatchMap_(refCast<const mappedPatchBase>(filmPatch_)),
    
    interphaseMassTransfer_(nullptr),
    
    transferRate_
    (
        volScalarField::Internal::New
        (
            "transferRate",
            mesh,
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    ),
    
    curTimeIndex_(-1)
{
    Info << "creating mass transfer autoptr" << endl;
    Info << "filmPatchi = " << filmPatchi_ << endl;
    Info << "filmPatch = " << filmPatch_.name() << endl;
    interphaseMassTransfer_ =
        interphaseMassTransferModel::New
        (
            dict,
            multicomponentFluid_,
            filmPatchi_,
            filmPatchMap_
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::bulkMassTransfer::addSupFields() const
{
    //- These are the fields/eqns that a source term is added to
    return wordList
    {
        specieName_,
        thermo_.rho()().name(),
        thermo_.he().name(),
        multicomponentFluid_.U.name()
    };
}


void Foam::fv::bulkMassTransfer::correct()
{
    //- This logic makes sure the correct() function is only called
    //  once per time step
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    curTimeIndex_ = mesh().time().timeIndex();
    
    //- Update interphaseMassTransferModel
    interphaseMassTransfer_->update();
    
    //- Grab mass transfer rate from model
    volScalarField mDotBulk = interphaseMassTransfer_->mDotBulk();
    
    //- Only apply negative mass transfer (i.e., transfer from bulk to film)
    //  on the bulk-side fvModel. Positive mass transfer is accounted for in
    //  the film-side fvModel.
    transferRate_ = neg0(mDotBulk) * mDotBulk;
}


template<class Type, class TransferRateFunc>
Foam::tmp<Foam::VolInternalField<Type>>
inline Foam::fv::bulkMassTransfer::filmToBulkMassTransferRate
(
    TransferRateFunc transferRateFunc,
    const dimensionSet& dimProp
) const
{
    //- Get reference to film side fvModels
    const Foam::fvModels& fvModels
    (
        fvModels::New
        (
            refCast<const fvMesh>(filmPatchMap_.nbrMesh())
        )
    );

    //- Get pointer to the film side mass transfer model
    const filmMassTransfer* filmBulkPtr = nullptr;

    forAll(fvModels, i)
    {
        if (isType<filmMassTransfer>(fvModels[i]))
        {
            filmBulkPtr = &refCast<const filmMassTransfer>(fvModels[i]);
        }
    }

    if (!filmBulkPtr)
    {
        FatalErrorInFunction
            << "Cannot find filmMassTransfer fvModel for the film region "
            << filmPatchMap_.nbrMesh().name()
            << exit(FatalError);
    }

    //- Create temporary field to hold transfer rate field from film side
    tmp<VolInternalField<Type>> tSu
    (
        VolInternalField<Type>::New
        (
            "Su",
            mesh(),
            dimensioned<Type>(dimProp/dimTime, Zero)
        )
    );

    //- Map data from film side to bulk side
    UIndirectList<Type>(tSu.ref(), filmPatch_.faceCells()) =
        filmPatchMap_.fromNeighbour
        (
            (filmBulkPtr->*transferRateFunc)()
        );

    return tSu/mesh().V();
}


void Foam::fv::bulkMassTransfer::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }
    
    //- Add source to density equation
    if (&rho == &thermo_.rho()())
    {
        eqn +=
            //- Contribution from mass transfer from film to bulk
            filmToBulkMassTransferRate<scalar>
            (
                &filmMassTransfer::rhoTransferRate,
                dimMass
            )
          //- Contribution from mass transfer from bulk to film
          + fvm::SuSp(transferRate_, eqn.psi());
    }
    //- Add source to species concentration equations
    else if (rho.name() == specieName_)
    {
        Info << "Adding source for specie " << specieName_ << endl;
        
        const volScalarField::Internal Yi
        (
            min(max(eqn.psi(), scalar(1e-6)), scalar(1))
        );
        
        eqn +=
            //- Contribution from mass transfer from film to bulk
            filmToBulkMassTransferRate<scalar>
            (
                &filmMassTransfer::rhoTransferRate,
                dimMass
            )
          //- Contribution from mass transfer from bulk to film
          + fvm::SuSp(transferRate_ / Yi, eqn.psi());
        
    }
    //- Add source to species concentration equations
    else if (&rho == &thermo_.he())
    {
        eqn +=
            //- Contribution from mass transfer from film to bulk
            filmToBulkMassTransferRate<scalar>
            (
                &filmMassTransfer::heTransferRate,
                dimEnergy
            )
          //- Contribution from mass transfer from bulk to film
          + fvm::SuSp(transferRate_, eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << rho.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::bulkMassTransfer::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    FatalErrorInFunction
        << "Support for field " << he.name() << " is not implemented"
        << exit(FatalError);
}


void Foam::fv::bulkMassTransfer::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&U == &multicomponentFluid_.U)
    {
        //- Add source to momentum equation
        eqn +=
            //- Contribution from mass transfer from film to bulk
            filmToBulkMassTransferRate<vector>
            (
                &filmMassTransfer::UTransferRate,
                dimMass*dimVelocity
            )      
          //- Contribution from mass transfer from bulk to film
          + fvm::SuSp(transferRate_, eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << U.name() << " is not implemented"
            << exit(FatalError);
    }
}


template<class Type, class FieldType>
inline Foam::tmp<Foam::Field<Type>> Foam::fv::bulkMassTransfer::TransferRate
(
    const FieldType& f
) const
{
    const labelList& faceCells = filmPatch_.faceCells();

    return tmp<Field<Type>>
    (
        new Field<Type>
        (
            UIndirectList<Type>
            (
                -1.0 * transferRate_ * mesh().V() * f,
                faceCells
            )
        )
    );
}


Foam::tmp<Foam::scalarField>
Foam::fv::bulkMassTransfer::rhoTransferRate() const
{
    return TransferRate<scalar>(oneField());
}


Foam::tmp<Foam::scalarField>
Foam::fv::bulkMassTransfer::heTransferRate() const
{
    return TransferRate<scalar>(thermo_.he()());
}


Foam::tmp<Foam::vectorField>
Foam::fv::bulkMassTransfer::UTransferRate() const
{
    return TransferRate<vector>(multicomponentFluid_.U());
}


void Foam::fv::bulkMassTransfer::topoChange(const polyTopoChangeMap&)
{
    transferRate_.setSize(mesh().nCells());
}


void Foam::fv::bulkMassTransfer::mapMesh(const polyMeshMap& map)
{
    transferRate_.setSize(mesh().nCells());
}


void Foam::fv::bulkMassTransfer::distribute(const polyDistributionMap&)
{
    transferRate_.setSize(mesh().nCells());
}


bool Foam::fv::bulkMassTransfer::movePoints()
{
    return true;
}


// ************************************************************************* //
