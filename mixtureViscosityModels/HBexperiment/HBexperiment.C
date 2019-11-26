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

#include "HBexperiment.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{
    defineTypeNameAndDebug(HBexperiment, 0);

    addToRunTimeSelectionTable
    (
        mixtureViscosityModel,
        HBexperiment,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixtureViscosityModels::HBexperiment::HBexperiment
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const word modelName
)
:
    mixtureViscosityModel(name, viscosityProperties, U, phi),
    HBexperimentCoeffs_(viscosityProperties.optionalSubDict(modelName +"Coeffs")),

    k_ 
    (
    	"k",
         dimensionSet(1, -1, -1, 0, 0),
         HBexperimentCoeffs_.lookup("k")	 
    ),
    n_
    (
     	"n",
	dimless,
	HBexperimentCoeffs_.lookup("n")
    ),
    tau0_
    (
        "tau0",
	dimensionSet(1, -1, -2, 0, 0),
	HBexperimentCoeffs_.lookup("tau0")
    ),
    tauCoeffs_
    (
        "tauCoeffs",
	dimensionSet(1, -1, -2, 0, 0),
	HBexperimentCoeffs_.lookup("tauCoeffs")
    ),
    tauExpo_
    (
     	"tauExpo",
	dimless,
	HBexperimentCoeffs_.lookup("tauExpo")
    ),
    kCoeffs_
    (
        "kCoeffs",
	dimensionSet(1, -1, -1, 0, 0),
	HBexperimentCoeffs_.lookup("kCoeffs")
    ),
    kExpo_
    (
     	"kExpo",
	dimless,
	HBexperimentCoeffs_.lookup("kExpo")
    ),
    mu0_
    (
    	"mu0",
        dimensionSet(1, -1, -1, 0, 0),
        HBexperimentCoeffs_.lookup("mu0")	
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
    checkAlpha_(checkParameters()),		//-checking the parameter
    U_(U)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mixtureViscosityModels::HBexperiment::mu(const volScalarField& muc) const
{

	//-calculating the strain-rate
	volScalarField strainRate
	(
		sqrt(2.0)*mag(symm(fvc::grad(U_)))
	);

	//-Yield stress calculation based on the experiment data
       /* volScalarField tauy
        (
		tauCoeffs_*exp(tauExpo_*alpha_)
	);*/
	volScalarField tauy		//-if alpha->0, tauy=0, tauCoeffs_= 4.78e7, tauExpo =5.5542
        (
		tauCoeffs_*pow(max(alpha_, scalar(0.0)),tauExpo_)
	);


        //-Consistency K calculation based on the experiment data
        volScalarField ky
        (
		kCoeffs_*exp(kExpo_*alpha_)
	);

	 //-rtone is for dimensioned settings
	 dimensionedScalar tone("tone", dimTime, 1.0);
         dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);
         

	 //------------------for switching to water----------------//
	// if (mag(alpha_) == 0.0)
	//    {
	//	n_ = ONE;
	//	ky = water;	//-viscosity of water
	//    }

         return min
    	   	(

			mu0_,
            		(tauy + ky*rtone*pow(tone*strainRate, n_))
           		/(max(strainRate, dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)))

    		);
}

bool Foam::mixtureViscosityModels::HBexperiment::checkParameters()
{
     	 dimensionedScalar zero("zero", dimless, 0);
         dimensionedScalar ONE("ONE", dimless, 1);
         dimensionedScalar water("water", dimensionSet(1, -1, -1, 0, 0), 0.001);
	 
	 forAll(alpha_, cellI)
		{
			if(alpha_[cellI] == 0.0)
			{
			//	ky = water;
				n_ = ONE;
			}
		}

	 return true;
}

bool Foam::mixtureViscosityModels::HBexperiment::read
(
    const dictionary& viscosityProperties
)
{
    mixtureViscosityModel::read(viscosityProperties);

    HBexperimentCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    HBexperimentCoeffs_.lookup("k") >> k_;
    HBexperimentCoeffs_.lookup("n") >> n_;
    HBexperimentCoeffs_.lookup("tau0") >> tau0_;
    HBexperimentCoeffs_.lookup("tauCoeffs") >> tauCoeffs_;
    HBexperimentCoeffs_.lookup("tauExpo") >> tauExpo_;
    HBexperimentCoeffs_.lookup("kCoeffs") >> kCoeffs_;
    HBexperimentCoeffs_.lookup("kExpo") >> kExpo_;
    HBexperimentCoeffs_.lookup("mu0") >> mu0_;
    return true;
}


// ************************************************************************* //
