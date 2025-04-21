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
    Foam::solvers::multicomponentFilm

Description
    Solver module for flow of compressible liquid films with species transport

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, Lagrangian particles,
    radiation, surface film etc. and constraining or limiting the solution.

SourceFiles
    multicomponentFilm.C

See also
    Foam::solvers::isothermalFilm

\*---------------------------------------------------------------------------*/

#ifndef multicomponentFilm_H
#define multicomponentFilm_H

#include "isothermalFilm.H"
#include "fluidMulticomponentThermo.H"
#include "combustionModel.H"
#include "filmMulticomponentThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                     Class multicomponentFilm Declaration
\*---------------------------------------------------------------------------*/

class multicomponentFilm
:
    public isothermalFilm
{

protected:

    // Thermophysical properties

        fluidMulticomponentThermo& thermo_;


    // Composition

        PtrList<volScalarField>& Y_;


    // Reactions

        autoPtr<combustionModel> reaction_;
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    // Thermophysical transport

        autoPtr<filmMulticomponentThermophysicalTransportModel> thermophysicalTransport;


public:

    // Public Data

        //- Reference to the fluid thermophysical properties
        const fluidMulticomponentThermo& thermo;

        //- Reference to the composition
        const PtrList<volScalarField>& Y;
        
        //- Reference to the combustion model
        const autoPtr<combustionModel>& reaction;
        

    //- Runtime type information
    TypeName("multicomponentFilm");


    // Constructors

        //- Construct from region mesh
        multicomponentFilm(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        multicomponentFilm(const multicomponentFilm&) = delete;


    //- Destructor
    virtual ~multicomponentFilm();


    // Member Functions

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();
        
        //- Predict the momentum and thermophysical transport modelling
        virtual void prePredictor();

        //- Correct the momentum and thermophysical transport modelling
        virtual void postCorrector();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const multicomponentFilm&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
