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

Class
    Foam::interphaseMassTransferModels::interphaseMassTransferModel

Description
    Abstract base class for interphaseMassTransfer models.
    
    This class of models is designed to transfer mass from a designated species
    across an interface between a multicomponentFluid and multicomponentFilm.
    Mass transfer fields are provided in both meshes.
    
    The interphaseMassTransferModel::update() function provides limiting
    functionality to prevent mass transfer to/from very thin films, and to
    limit the mass transfer rate based on the total mass in each phase and
    the current time step. This function should be called in any derived
    classes after the rates are defined in the derived class.

SourceFiles
    interphaseMassTransferModel.C

\*---------------------------------------------------------------------------*/

#ifndef interphaseMassTransferModel_H
#define interphaseMassTransferModel_H

#include "dictionary.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "multicomponentFluid.H"
#include "multicomponentFilm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class interphaseMassTransferModel Declaration
\*---------------------------------------------------------------------------*/

class interphaseMassTransferModel
{
    
protected:

    // Protected data
        
        //- Mass transfer subdictionary in fvModels dict
        const dictionary& massTransferModelCoeffs_;
        
        //- Pointer to bulk phase instance
        const solvers::multicomponentFluid& fluid_;
        
        //- Mapper between bulk and film
        const mappedPatchBase& patchMap_;
        
        //- Pointer to film phase instance
        const solvers::multicomponentFilm& film_;
        
        //- Pointer to bulk mesh instance
        const fvMesh& bulkMesh_;
        
        //- Pointer to film mesh instance
        const fvMesh& filmMesh_;
        
        //- Index of mass transfer patch in bulk mesh
        const label bulkPatchID_;
        
        //- Index of mass transfer patch in film mesh
        const label filmPatchID_;
        
        //- The name of the transferred specie
        const word specieName_;
        
        //- Index of transferred specie in bulk composition
        const label bulkSpecieID_;
        
        //- Index of transferred specie in film composition
        const label filmSpecieID_;
        
        //- Mass transfer field defined on the bulk mesh [kg/m^3/s]
        volScalarField mDotBulk_;
        
        //- Mass transfer field defined on the film mesh [kg/m^3/s]
        volScalarField mDotFilm_;
        
        //- Interfacial mass transfer rate [kg/m^2/s], defined on the
        //  bulk-side interface patch
        scalarField R_;

public:

    //- Runtime type information
    TypeName("interphaseMassTransferModel");
    
    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        interphaseMassTransferModel,
        dictionary,
        (
            const dictionary& dict,
            const solvers::multicomponentFluid& fluid,
            const label& bulkPatchID,
            const mappedPatchBase& patchMap
        ),
        (dict, fluid, bulkPatchID, patchMap)
    );


    // Constructors

        //- Construct from components
        interphaseMassTransferModel
        (
            const word& type,
            const dictionary& dict,
            const solvers::multicomponentFluid& fluid,
            const label& bulkPatchID,
            const mappedPatchBase& patchMap
        );
        
        //- Disallow default bitwise copy construction
        interphaseMassTransferModel(const interphaseMassTransferModel&);
        
    // Selectors

        //- Return a reference to the selected interphaseMassTransfer model
        static autoPtr<interphaseMassTransferModel> New
        (
            const dictionary& dict,
            const solvers::multicomponentFluid& fluid,
            const label& bulkPatchID,
            const mappedPatchBase& patchMap
        );

    //- Destructor
    virtual ~interphaseMassTransferModel()
    {}


    // Member Functions
        
        //- Return const access to the bulk-side mass transfer field
        const volScalarField& mDotBulk() const
        {
            return mDotBulk_;
        }
        
        //- Return const access to the film-side mass transfer field
        const volScalarField& mDotFilm() const
        {
            return mDotFilm_;
        }
        
        //- Return const access to the interfacial mass transfer field,
        //  defined on the gas-side interface patch
        const scalarField& R() const
        {
            return R_;
        }        

        //- Update the mass transfer field
        virtual void update() = 0;

        //- Read the mass transfer dictionary
        virtual bool read() = 0;
        
        
    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const interphaseMassTransferModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
