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
    Foam::massTransferRateCoefficientModels::massTransferRateCoefficientModel

Description
    Abstract base class for mass transfer rate coefficient models.

SourceFiles
    massTransferRateCoefficientModel.C

\*---------------------------------------------------------------------------*/

#ifndef massTransferRateCoefficientModel_H
#define massTransferRateCoefficientModel_H

#include "dictionary.H"
#include "volFields.H"
#include "multicomponentFilm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class massTransferRateCoefficientModel Declaration
\*---------------------------------------------------------------------------*/

class massTransferRateCoefficientModel
{
    
protected:

    // Protected data
        
        //- Rate coefficient subdictionary in fvModels dict
        const dictionary& rateModelCoeffs_;
        
        //- Pointer to film mesh instance
        const fvMesh& mesh_;
        
        //- Index of transfer patch in mesh
        const label patchID_;
        
        //- Rate coefficient
        volScalarField k_;

public:

    //- Runtime type information
    TypeName("massTransferRateCoefficientModel");
    
    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        massTransferRateCoefficientModel,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const label& patchID
        ),
        (dict, mesh, patchID)
    );


    // Constructors

        //- Construct from components
        massTransferRateCoefficientModel
        (
            const word& type,
            const dictionary& dict,
            const fvMesh& mesh,
            const label& patchID
        );
        
        //- Disallow default bitwise copy construction
        massTransferRateCoefficientModel(const massTransferRateCoefficientModel&);
        
    // Selectors

        //- Return a reference to the selected rate coefficient model
        static autoPtr<massTransferRateCoefficientModel> New
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const label& patchID
        );

    //- Destructor
    virtual ~massTransferRateCoefficientModel()
    {}


    // Member Functions
        
        //- Calculate and return the rate coefficient
        virtual const volScalarField& k() = 0;        
        
    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const massTransferRateCoefficientModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
