/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "functionObjectTemplate.H"
#include "fvCFD.H"
#include "unitConversion.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(averageFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    averageFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = c8ec27a90dcc8f0af0923ba8ed27a3ec5afe1243
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void average_c8ec27a90dcc8f0af0923ba8ed27a3ec5afe1243(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const fvMesh& averageFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

averageFunctionObject::averageFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

averageFunctionObject::~averageFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool averageFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        Info<<"read average sha1: c8ec27a90dcc8f0af0923ba8ed27a3ec5afe1243\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool averageFunctionObject::execute()
{
    if (false)
    {
        Info<<"execute average sha1: c8ec27a90dcc8f0af0923ba8ed27a3ec5afe1243\n";
    }

//{{{ begin code
    #line 155 "/home/shakil/OpenFOAM/shakil-5.x/run/case/ARC/rheology/Herschel/HBexperiment/dahl_knvary/system/controlDict.functions.extraInfo"
const volVectorField& Umix = mesh().lookupObject<volVectorField>("U");
             //const volVectorField& Uwater = mesh().lookupObject<volVectorField>("U.water");
	     Info << "max U = " << max(mag(Umix)).value() << ", min  U = " << min(mag(Umix)).value() << endl;
             const volScalarField& p = mesh().lookupObject<volScalarField>("p");
	     Info << "p min/max = " << min(p).value() << ", " << max(p).value() << endl;
	     //const volScalarField& kmix = mesh().lookupObject<volScalarField>("km");
	     //Info << "p min/max = " << min(kmix).value() << ", " << max(kmix).value() << endl;
//}}} end code

    return true;
}


bool averageFunctionObject::write()
{
    if (false)
    {
        Info<<"write average sha1: c8ec27a90dcc8f0af0923ba8ed27a3ec5afe1243\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool averageFunctionObject::end()
{
    if (false)
    {
        Info<<"end average sha1: c8ec27a90dcc8f0af0923ba8ed27a3ec5afe1243\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

