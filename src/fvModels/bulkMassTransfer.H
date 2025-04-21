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

Class
    Foam::fv::bulkMassTransfer

Description
    Model for mass transfer between a bulk phase and a liquid film. This class
    computes quantities and returns source terms on the bulk phase side, and
    requires a corresponding twin model in the liquid film phase.

Usage
    Example usage:
    \verbatim
    bulkMassTransfer
    {
        type                bulkMassTransfer;
        
        select              all;
        
        filmPatch           film;
        
        specie              CO2;
        
        transferRateCoeff   1e-6;
    }
    \endverbatim

SourceFiles
    bulkMassTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef bulkMassTransfer_H
#define bulkMassTransfer_H

#include "fvModel.H"
#include "fvc.H"
#include "multicomponentFluid.H"
#include "interphaseMassTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class bulkMassTransfer Declaration
\*---------------------------------------------------------------------------*/

class bulkMassTransfer
:
    public fvModel
{
    // Private Data
    
        //- The multi-component fluid model
        const solvers::multicomponentFluid& multicomponentFluid_;
        
        //- Reference to the transferred specie thermo package
        const fluidMulticomponentThermo& thermo_;
        
        //- The name of the transferred specie
        const word specieName_;

        //- Film patch name
        const word filmPatchName_;

        //- Index of the film patch in the mesh
        const label filmPatchi_;
        
        //- Reference to mass transfer patch
        const polyPatch& filmPatch_;
        
        //- Patch mapper between two phases
        const mappedPatchBase& filmPatchMap_;
        
        //- Reference to interphase mass transfer model
        autoPtr<interphaseMassTransferModel> interphaseMassTransfer_;
    
        //- Transfer rate field
        volScalarField::Internal transferRate_;
        
        //- Current time index (used for updating)
        mutable label curTimeIndex_;


    // Private Member Functions

        //- Return the transfer rate from the film transferRateFunc
        template<class Type, class TransferRateFunc>
        inline tmp<VolInternalField<Type>> filmToBulkMassTransferRate
        (
            TransferRateFunc transferRateFunc,
            const dimensionSet& dimProp
        ) const;

        //- Return the bulk-to-film transfer rate of field f
        template<class Type, class FieldType>
        inline tmp<Field<Type>> TransferRate(const FieldType& f) const;


public:

    //- Runtime type information
    TypeName("bulkMassTransfer");


    // Constructors

        //- Construct from explicit source name and mesh
        bulkMassTransfer
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        bulkMassTransfer
        (
            const bulkMassTransfer&
        ) = delete;


    // Member Functions

            label filmPatchIndex() const
            {
                return filmPatchi_;
            }
            
            const interphaseMassTransferModel& interphaseMassTransfer() const
            {
                return interphaseMassTransfer_();
            }


        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Correct

            //- Solve the film and update the sources
            virtual void correct();


        // Sources

            //- Add implicit contribution to phase continuity equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add implicit contribution to phase energy equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add implicit contribution to mixture momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;


        // Transfer to film

            //- Return the mass transfer rate
            tmp<scalarField> rhoTransferRate() const;

            //- Return the energy transfer rate
            tmp<scalarField> heTransferRate() const;

            //- Return the momentum transfer rate
            tmp<vectorField> UTransferRate() const;


        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const bulkMassTransfer&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
