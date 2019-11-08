/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::mixtureViscosityModels::HBexperiment

Description
     Viscosity correction model for a generic power-law HBexperiment.

SourceFiles
    HBexperiment.C
Default Herschel-Bulkley model where ty and K is dependent on alpha
ty = f(alpha) and K = f(alpha), mu = f(alpha, strRate)
\*---------------------------------------------------------------------------*/

#ifndef HBexperiment_H
#define HBexperiment_H

#include "mixtureViscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class incompressibleTwoPhaseInteractingMixture;

namespace mixtureViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class HBexperiment Declaration
\*---------------------------------------------------------------------------*/

class HBexperiment
:
    public mixtureViscosityModel
{
protected:

    // Protected data

        //- Dictionary
        dictionary HBexperimentCoeffs_;

	//-consistency index
	dimensionedScalar k_;

	//-flow/power-law index
	dimensionedScalar n_;

	//-Yield-stress
	dimensionedScalar tau0_;
	
        //-Yield stress Coeffs
	dimensionedScalar tauCoeffs_;

	//-Yield stress exponent
	dimensionedScalar tauExpo_;

	//-Viscosity at the zero shear rate
	dimensionedScalar mu0_;
	
        //- phase fraction
        const volScalarField& alpha_;

	//-velocity field 
	const volVectorField& U_;


public:

    //- Runtime type information
    TypeName("HBexperiment");


    // Constructors

        //- Construct from components
        HBexperiment
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const word modelName=typeName
        );


    //- Destructor
    ~HBexperiment()
    {}


    // Member Functions

        //- Return the mixture viscosity
        //  given the viscosity of the continuous phase
        tmp<volScalarField> mu(const volScalarField& muc) const;

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace mixtureViscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
