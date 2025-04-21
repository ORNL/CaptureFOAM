/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "ionizedArrheniusReactionRate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::ionizedArrheniusReactionRate::
ionizedArrheniusReactionRate
(
    const speciesTable& species,
    const objectRegistry& ob,
    const dimensionSet& dims,
    const dictionary& dict
)
:
    ArrheniusReactionRate(species, dims, dict),
    ob_(ob),
    ionName_(species["OH-"]),
    A1_(dict.lookup<scalar>("A1", dimVolume/dimMoles/dimTemperature/dimTemperature)),
    A2_(dict.lookup<scalar>("A2", dimVolume/dimMoles/dimTemperature)),
    A3_(dict.lookup<scalar>("A3", dimVolume/dimMoles))
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::ionizedArrheniusReactionRate::preEvaluate() const
{
    ArrheniusReactionRate::preEvaluate();
}


inline void Foam::ionizedArrheniusReactionRate::postEvaluate() const
{
    ArrheniusReactionRate::postEvaluate();
}


inline Foam::scalar Foam::ionizedArrheniusReactionRate::operator()
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    const scalar Cion = c[ionName_];
    const scalar beta_co2 = A1_*pow(T,2) + A2_*T + A3_;

    return ArrheniusReactionRate::operator()(p, T, c, li)*exp(beta_co2*Cion);
}


inline Foam::scalar Foam::ionizedArrheniusReactionRate::ddT
(
    const scalar p,
    const scalar T,
    const scalarField& c,
    const label li
) const
{
    const scalar Cion = c[ionName_];
    const scalar beta_co2 = A1_*pow(T,2) + A2_*T + A3_;
    const scalar dbetadT = 2.0*A1_*T + A2_;

    const scalar dk1dT = ArrheniusReactionRate::ddT(p, T, c, li)*exp(beta_co2*Cion);

     const scalar dk2dT = ArrheniusReactionRate::operator()(p, T, c, li)
	    * dbetadT * Cion * exp(beta_co2*Cion);

    return dk1dT + dk2dT;
}


inline void Foam::ionizedArrheniusReactionRate::write(Ostream& os) const
{
    ArrheniusReactionRate::write(os);
    writeEntry(os, "ion specie", ionName_);
}


inline Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ionizedArrheniusReactionRate& arr
)
{
    arr.write(os);
    return os;
}


// ************************************************************************* //
