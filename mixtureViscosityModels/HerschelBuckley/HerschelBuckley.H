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

     Shakil: I am adding HB model considering the constant yield stress at this version
	     similar to the OpenFOAM default version

	     Parameters I need as input:
	     mu0  = viscosity at zero shear-rate = [1 -1 -1 0 0 0 0]
	     tau0 = yield stress                 = [1 -1 -2 0 0 0 0]
             k   = consistency index		 = [1 -1 -1 0 0 0 0]
             n   = dimensionless
                 
	     

SourceFiles
    HerschelBuckley.C

\*---------------------------------------------------------------------------*/

#ifndef HerschelBuckley_H
#define HerschelBuckley_H

#include "plastic.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace mixtureViscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class BinghamPlastic Declaration
\*---------------------------------------------------------------------------*/

class HerschelBuckley
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

     


public:

    //- Runtime type information
    TypeName("HerschelBuckley");


    // Constructors

        //- Construct from components
        HerschelBuckley
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~HerschelBuckley()
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
