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

Class
    Foam::mixtureViscosityModels::HerschelBuckley

Description
     Viscosity correction model for Herschel Buckley plastics.

     strain-rate = sqrt(2.0)*mag(symm(grad(U)))

     Shakil: I am adding HB model considering the constant yield stress 
     Ref.: Radhakrishnan et al., Rheological characterisation of concentrated domestic slurry
           Water Research (141), 2018
           Dapelo et al., Advances in Engineering Software (117), 2018

	     Parameters I need as input:
	     mu0  = viscosity at zero shear-rate = [1 -1 -1 0 0 0 0]
	     tau0 = yield stress                 = [1 -1 -2 0 0 0 0]
             k   = consistency index		 = [1 -1 -1 0 0 0 0]
             n   = dimensionless
                 
	     

SourceFiles
    HerschelBuckley.C

\*---------------------------------------------------------------------------*/

#ifndef HerschelBulkley_H
#define HerschelBulkley_H

#include "plastic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class BinghamPlastic Declaration
\*---------------------------------------------------------------------------*/

class HerschelBulkley
:
    public plastic
{
protected:

    // Protected data

	//-consistency index
	dimensionedScalar k_;
	
	//-power/flow index
        dimensionedScalar n_;
	
	//-yield stress        
	dimensionedScalar tau0_;
        
	//-viscosity at zero shear-rate
	dimensionedScalar mu0_;
      
        //- Velocity
        const volVectorField& U_;

	//-for printing mu
	//const volScalarField& muMix;
  


public:

    //- Runtime type information
    TypeName("HerschelBulkley");


    // Constructors

        //- Construct from components
        HerschelBulkley
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~HerschelBulkley()
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
