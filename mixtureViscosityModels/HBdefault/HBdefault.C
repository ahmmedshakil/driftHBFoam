/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2017 OpenFOAM Foundation
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

#include "HBdefault.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(HBdefault, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        HBdefault,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::HBdefault::HBdefault
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word modelName
)
:
    mixtureViscosityModel(name, viscosityProperties, U, phi),
    HBdefaultCoeffs_(viscosityProperties.optionalSubDict(modelName +"Coeffs")),

    k_ 
    (
    	"k",
         dimensionSet(1, -1, -1, 0, 0),
         HBdefaultCoeffs_.lookup("k")	 
    ),
    n_
    (
     	"n",
	dimless,
	HBdefaultCoeffs_.lookup("n")
    ),
    tau0_
    (
        "tau0",
	dimensionSet(1, -1, -2, 0, 0),
	HBdefaultCoeffs_.lookup("tau0")
    ),
    mu0_
    (
    	"mu0",
        dimensionSet(1, -1, -1, 0, 0),
        HBdefaultCoeffs_.lookup("mu0")	
    ),
    alpha_
    (
        U.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName
            (
                viscosityProperties.lookupOrDefault<word>("alpha", "alpha"),
                viscosityProperties.dictName()
            )
        )
    ),
    U_(U)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::HBdefault::mu(const volScalarField& muc) const
{

	//-calculating the strain-rate
	volScalarField strainRate
	(
		sqrt(2.0)*mag(symm(fvc::grad(U_)))
	);


	 //-rtone is for dimensioned settings
	 dimensionedScalar tone("tone", dimTime, 1.0);
         dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

         return min
    	   	(

			mu0_,
            		(tau0_ + k_*rtone*pow(tone*strainRate, n_))
           		/(max(strainRate, dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)))

    		);
}



bool Foam::mixtureViscosityModels::HBdefault::read
(
    const dictionary& viscosityProperties
)
{
    mixtureViscosityModel::read(viscosityProperties);

    HBdefaultCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    HBdefaultCoeffs_.lookup("k") >> k_;
    HBdefaultCoeffs_.lookup("n") >> n_;
    HBdefaultCoeffs_.lookup("tau0") >> tau0_;
    HBdefaultCoeffs_.lookup("mu0") >> mu0_;
    return true;
}


// ************************************************************************* //
