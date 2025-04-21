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

#include "multicomponentFilm.H"
#include "fvmDdt.H"
#include "fvmDiv.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::multicomponentFilm::thermophysicalPredictor()
{
    //- Adopted from multicomponentFluid::thermophysicalPredictor()
    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            alphaRhoPhi,
            mesh.schemes().div("div(alphaRhoPhi,Yi_h)")
        )
    );
    
    reaction_->correct();

    forAll(Y, i)
    {
        if (thermo_.solveSpecie(i))
        {
            volScalarField& Yi = Y_[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(alpha, rho, Yi)
              + mvConvection->fvmDiv(alphaRhoPhi, Yi)
              + thermophysicalTransport->divj(Yi)
             ==
                reaction_->R(Yi)
              + fvModels().source(alpha, rho, Yi)
            );

            YiEqn.relax();

            fvConstraints().constrain(YiEqn);

            YiEqn.solve("Yi");

            fvConstraints().constrain(Yi);
        }
    }

    thermo_.normaliseY();

    //- Rest is same as film::thermophysicalPredictor()
    volScalarField& he = thermo_.he();

    fvScalarMatrix heEqn
    (
        fvm::ddt(alpha, rho, he) + fvm::div(alphaRhoPhi, he)
      - fvm::Sp(contErr(), he)
      + thermophysicalTransport->divq(he)
     ==
        reaction_->Qdot()
      + fvModels().source(alpha, rho, he)
    );

    heEqn.relax();

    fvConstraints().constrain(heEqn);

    heEqn.solve();

    fvConstraints().constrain(he);

    thermo_.correct();
}


// ************************************************************************* //
