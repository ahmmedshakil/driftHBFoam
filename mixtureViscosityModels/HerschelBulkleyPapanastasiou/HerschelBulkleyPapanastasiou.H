/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::viscosityModels::HerschelBulkleyPapanastasiou

Description
     Herschel-Bulkley non-Newtonian viscosity model.

SourceFiles
    HerschelBulkleyPapanastasiou.C
Example

transportModel  HerschelBulkleyPapanastasiou;
HerschelBulkleyPapanastasiouCoeffs
{
    k               k		[ 0 2 -1 0 0 0 0 ]	1;
	n				n		[ 0 0  0 0 0 0 0 ]	0.5;
	tau0			tau0	[ 0 2 -2 0 0 0 0 ]	1.0; 
	nu0				nu0		[ 0 2 -1 0 0 0 0 ]	1e4; 
	m				m		[ 0 0 0 0 0 0 0 ]	1e2; 
	sr0				sr0		[0 0 -1 0 0 0 0]	 1.0;  // Relative shear rate. You can set it up equal to one.
}

\*---------------------------------------------------------------------------*/

#ifndef HerschelBulkleyPapanastasiou_H
#define HerschelBulkleyPapanastasiou_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class HerschelBulkleyPapanastasiou Declaration
\*---------------------------------------------------------------------------*/

class HerschelBulkleyPapanastasiou
:
    public viscosityModel
{
    // Private data

        dictionary HerschelBulkleyPapanastasiouCoeffs_;

        dimensionedScalar k_;
        dimensionedScalar n_;
        dimensionedScalar tau0_;
        dimensionedScalar nu0_;
        dimensionedScalar m_;
        dimensionedScalar sr0_; // Relative shear rate. You can set it up equal to one.

        volScalarField nu_;


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;


public:

    //- Runtime type information
    TypeName("HerschelBulkleyPapanastasiou");


    // Constructors

        //- Construct from components
        HerschelBulkleyPapanastasiou
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~HerschelBulkleyPapanastasiou()
    {}


    // Member Functions

        //- Return the laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        void correct()
        {
            nu_ = calcNu();
        }

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
