/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
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

#include "HerschelBuckley.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(HerschelBuckley, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        HerschelBuckley,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::HerschelBuckley::HerschelBuckley
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    plastic(name, viscosityProperties, U, phi, typeName),

    k_
    (
        "k",
        dimensionSet(1, -1, -1, 0, 0),
        plasticCoeffs_
    ),


    n_
    (
        "n",
        dimless,
        plasticCoeffs_
    ),

    tau0_
    (
        "tau0",
        dimensionSet(1, -1, -2, 0, 0),
        plasticCoeffs_
    ),

    mu0_
    (
        "mu0",
        dimensionSet(1, -1, -1, 0, 0),
        plasticCoeffs_
    ),

    
    U_(U)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//-mixture mu----------------------

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::HerschelBuckley::mu
(
    const volScalarField& muc
) const
{

	//-calculating the strain-rate
	volScalarField strainRate
	(
		sqrt(2.0)*mag(symm(fvc::grad(U_)))
	);


	 //-rtone is for dimensioned settings
	 dimensionedScalar tone("tone", dimTime, 1.0);
         dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);


	//-temporary dispersedPhase viscosity
/*	volScalarField muHB

	(
		min
        	(
            		mu0_,
            		(tau0_ + k_*rtone*pow(tone*strainRate, n_))
           		/(max(strainRate, dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)))
        	)
	);*/


        //-plastic viscosity based on the continuous viscosity
        //-   volScalarField mup(plastic::mu(muc));

           

    return min
    	   	(

			mu0_,
            		(tau0_ + k_*rtone*pow(tone*strainRate, n_))
           		/(max(strainRate, dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)))

        	//	muHB
      		//	+ mup,
		//
        	//	muMax_
    		);
}


bool Foam::mixtureViscosityModels::HerschelBuckley::read
(
    const dictionary& viscosityProperties
)
{
    plastic::read(viscosityProperties);

    plasticCoeffs_.lookup("k") >> k_;
    plasticCoeffs_.lookup("tau0") >> tau0_;
    plasticCoeffs_.lookup("mu0") >> mu0_;
    plasticCoeffs_.lookup("n") >> n_;

    return true;
}


// ************************************************************************* //
