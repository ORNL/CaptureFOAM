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

#include "filmMassTransfer.H"
#include "bulkMassTransfer.H"
#include "mappedFvPatchBaseBase.H"
#include "multicomponentFluid.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(filmMassTransfer, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            filmMassTransfer,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::filmMassTransfer::filmMassTransfer
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),
    film_(mesh.lookupObject<solvers::multicomponentFilm>(solver::typeName)),
    specieName_(dict.lookup("specie")),
    specieIndex_(film_.thermo.species()[specieName_]),
    specie_(film_.thermo.Y(specieIndex_)),
    curTimeIndex_(-1),
    transferRate_
    (
        volScalarField::Internal::New
        (
            "transferRate",
            mesh,
            dimensionedScalar(dimMass/dimVolume/dimTime, 0)
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::filmMassTransfer::addSupFields() const
{
    return wordList
    {
        film_.alpha.name(),
        specie_.name(),
        film_.thermo.he().name(),
        film_.U.name()
    };
}


void Foam::fv::filmMassTransfer::correct()
{
    //- This logic makes sure the correct() function is only called
    //  once per time step
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    curTimeIndex_ = mesh().time().timeIndex();
    
    const solvers::multicomponentFluid& fluid_
    (
        film_.surfacePatchMap().nbrMesh().lookupObject<solvers::multicomponentFluid>
        (
            solver::typeName
        )
    );
    
    const Foam::fv::bulkMassTransfer& bulkFvModel(this->bulkFilm(fluid_.fvModels()));
    
    //- Grab mass transfer rate from model
    const volScalarField mDotFilm = bulkFvModel.interphaseMassTransfer().mDotFilm();
    
    //- Only apply positive mass transfer (i.e., transfer from film to bulk)
    //  on the film-side fvModel. Negative mass transfer is accounted for in
    //  the bulk-side fvModel.
    transferRate_ = neg0(mDotFilm) * mDotFilm;
}


const Foam::fv::bulkMassTransfer& Foam::fv::filmMassTransfer::bulkFilm
(
    const Foam::fvModels& fvModels
) const
{
    const bulkMassTransfer* bulkFilmPtr = nullptr;

    forAll(fvModels, i)
    {
        if (isType<bulkMassTransfer>(fvModels[i]))
        {
            const bulkMassTransfer& bulkFilm
            (
                refCast<const bulkMassTransfer>(fvModels[i])
            );

            if
            (
                bulkFilm.filmPatchIndex()
             == film_.surfacePatchMap().nbrFvPatch().index()
            )
            {
                bulkFilmPtr = &bulkFilm;
            }
        }
    }

    if (!bulkFilmPtr)
    {
        FatalErrorInFunction
            << "Cannot find bulkMassTransfer fvModel for this film "
               "in bulk region " << film_.surfacePatchMap().nbrMesh().name()
            << exit(FatalError);
    }

    return *bulkFilmPtr;
}


template<class Type, class TransferRateFunc>
Foam::tmp<Foam::VolInternalField<Type>>
inline Foam::fv::filmMassTransfer::bulkToFilmTransferRate
(
    TransferRateFunc transferRateFunc,
    const dimensionSet& dimProp
) const
{    
    const Foam::fvModels& fvModels =
        fvModels::New(film_.surfacePatchMap().nbrMesh());

    //- Create temporary field to hold transfer rate field from bulk side
    tmp<VolInternalField<Type>> tSu
    (
        VolInternalField<Type>::New
        (
            "Su",
            mesh(),
            dimensioned<Type>(dimProp/dimTime, Zero)
        )
    );

    //- Map transfer field from bulk to film
    UIndirectList<Type>(tSu.ref(), film_.surfacePatch().faceCells()) =
        film_.surfacePatchMap().fromNeighbour
        (
            (bulkFilm(fvModels).*transferRateFunc)()
        );

    return tSu/mesh().V();
}


void Foam::fv::filmMassTransfer::addSup
(
    const volScalarField& rho,
    const volScalarField& alpha,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    //- Add source term to film continuity equation
    if (&alpha == &film_.alpha)
    {
        Info << "Adding film mass transfer source to continuity equation" << endl;
        
        volScalarField::Internal bulkToFilmMassXfer =
            bulkToFilmTransferRate<scalar>
            (
                &bulkMassTransfer::rhoTransferRate,
                dimMass
            );
        
        eqn +=
            //- Get bulk-to-film contribution from bulk fvModel
            bulkToFilmTransferRate<scalar>
            (
                &bulkMassTransfer::rhoTransferRate,
                dimMass
            )
          //- Subtract film-to-bulk contribution (zero since we
          //  set transferRate_ = 0.0)
          + fvm::Sp(transferRate_, eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << alpha.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::filmMassTransfer::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }


    if (&he == &specie_)
    {
        Info << "Adding film mass transfer source to species equation for specie "
             << specie_.name() << endl;
        
        eqn +=
            //- Get bulk-to-film contribution from bulk fvModel
            bulkToFilmTransferRate<scalar>
            (
                &bulkMassTransfer::rhoTransferRate,
                dimMass
            )
          //- Subtract film-to-bulk contribution (zero since we
          //  set transferRate_ = 0.0)
          //  NOTE: if we ever want non-zero film-to-bulk mass transfer
          //  we will need to double check this formulation.
          + fvm::Sp(alpha() * transferRate_, eqn.psi());
    }
    //- Add source term to film energy equation
    else if (&he == &film_.thermo.he())
    {
        Info << "Adding film mass transfer source to energy equation" << endl;
        
        eqn +=
            //- Get bulk-to-film contribution from bulk fvModel
            bulkToFilmTransferRate<scalar>
            (
                &bulkMassTransfer::heTransferRate,
                dimEnergy
            )
          //- Subtract film-to-bulk contribution (zero since we
          //  set transferRate_ = 0.0)
          //  NOTE: if we ever want non-zero film-to-bulk mass transfer
          //  we will need to double check this formulation.
          + fvm::Sp(alpha() * transferRate_, eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << he.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::filmMassTransfer::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }
    
    if (&U == &film_.U)
    {
    Info << "Adding film mass transfer source to momentum equation" << endl;

    //- Add source term to film momentum equation
    eqn +=
        //- Get bulk-to-film contribution from bulk fvModel
        bulkToFilmTransferRate<vector>
        (
            &bulkMassTransfer::UTransferRate,
            dimMomentum
        )
      //- Subtract film-to-bulk contribution (zero since we
      //  set transferRate_ = 0.0)
      //  NOTE: if we ever want non-zero film-to-bulk mass transfer
      //  we will need to double check this formulation.
      + fvm::Sp(alpha() * transferRate_, eqn.psi());
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << U.name() << " is not implemented"
            << exit(FatalError);
    }
}


template<class Type, class FieldType>
inline Foam::tmp<Foam::Field<Type>> Foam::fv::filmMassTransfer::TransferRate
(
    const FieldType& f
) const
{
    const labelList& faceCells = film_.surfacePatch().faceCells();

    return tmp<Field<Type>>
    (
        new Field<Type>
        (
            UIndirectList<Type>
            (
                //  NOTE: if we ever want non-zero film-to-bulk mass transfer
                //  we will need to double check this formulation.
                transferRate_ * film_.alpha() * mesh().V() * f,
                faceCells
            )
        )
    );
}


Foam::tmp<Foam::scalarField>
Foam::fv::filmMassTransfer::transferRate() const
{
    return TransferRate<scalar>(oneField());
}


Foam::tmp<Foam::scalarField>
Foam::fv::filmMassTransfer::rhoTransferRate() const
{
    return TransferRate<scalar>(oneField());
}


Foam::tmp<Foam::scalarField>
Foam::fv::filmMassTransfer::heTransferRate() const
{
    return TransferRate<scalar>(film_.thermo.he()());
}


Foam::tmp<Foam::vectorField>
Foam::fv::filmMassTransfer::UTransferRate() const
{
    return TransferRate<vector>(film_.U());
}


void Foam::fv::filmMassTransfer::topoChange(const polyTopoChangeMap&)
{
    transferRate_.setSize(mesh().nCells());
}


void Foam::fv::filmMassTransfer::mapMesh(const polyMeshMap& map)
{
    transferRate_.setSize(mesh().nCells());
}


void Foam::fv::filmMassTransfer::distribute(const polyDistributionMap&)
{
    transferRate_.setSize(mesh().nCells());
}


bool Foam::fv::filmMassTransfer::movePoints()
{
    return true;
}


// ************************************************************************* //
