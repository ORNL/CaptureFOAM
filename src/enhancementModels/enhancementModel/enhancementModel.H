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
    Foam::enhancementModels::enhancementModel

Description
    Abstract base class for enhancement models.
    
    Provides an enhancement factor field for mass transfer where the effective
    rate is increased owing to the effects of chemical reactions.

SourceFiles
    enhancementModel.C

\*---------------------------------------------------------------------------*/

#ifndef enhancementModel_H
#define enhancementModel_H

#include "dictionary.H"
#include "volFields.H"
#include "multicomponentFilm.H"
#include "multicomponentFluid.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class enhancementModel Declaration
\*---------------------------------------------------------------------------*/

class enhancementModel
{
    
protected:

    // Protected data
        
        //- Mass transfer subdictionary in fvModels dict
        const dictionary& massTransferModelCoeffs_;
        
        //- Pointer to film phase instance
        const solvers::multicomponentFilm& film_;
        
        //- Pointer to film mesh instance
        const fvMesh& filmMesh_;
        
        //- Index of transferred specie in film composition
        const label filmSpecieID_;
        
        //- Enhancement factor
        volScalarField E_;
        
        //- Apparent reaction rate
        volScalarField kApp_;

	//- Pointer to bulk phase instance
        const solvers::multicomponentFluid& fluid_;

	//- Pointer to bulk mesh instance
        const fvMesh& bulkMesh_;

	//- Index of mass transfer patch in bulk mesh
        const label bulkPatchID_;

public:

    //- Runtime type information
    TypeName("enhancementModel");
    
    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        enhancementModel,
        dictionary,
        (
            const dictionary& dict,
            const solvers::multicomponentFilm& film,
            const solvers::multicomponentFluid& fluid,
            const label& filmSpecieID,
            const label& bulkPatchID
        ),
        (dict, film, fluid, filmSpecieID, bulkPatchID)
    );


    // Constructors

        //- Construct from components
        enhancementModel
        (
            const word& type,
            const dictionary& dict,
            const solvers::multicomponentFilm& film,
            const solvers::multicomponentFluid& fluid,
            const label& filmSpecieID,
            const label& bulkPatchID
        );
        
        //- Disallow default bitwise copy construction
        enhancementModel(const enhancementModel&);
        
    // Selectors

        //- Return a reference to the selected enhancement model
        static autoPtr<enhancementModel> New
        (
            const dictionary& dict,
            const solvers::multicomponentFilm& film,
            const solvers::multicomponentFluid& fluid,
            const label& filmSpecieID,
            const label& bulkPatchID
        );

    //- Destructor
    virtual ~enhancementModel()
    {}


    // Member Functions
        
        //- Return const access to the enhancement factor
        const volScalarField& E() const
        {
            return E_;
        }
        
        //- Calculate and return the apparent reaction rate
        const volScalarField& kApp();

        //- Update the mass transfer field
        virtual void update() = 0;

        //- Read the enhancement factor dictionary
        virtual bool read() = 0;
        
        
    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const enhancementModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
