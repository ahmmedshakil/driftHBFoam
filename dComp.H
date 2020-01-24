/*-solving the transport equation for the tracer test: added by shakil
Calculation of the compression term based on
R.Burger et al., Water Science Technology 68(1) 2013s
*/
Info<< "Calculating the compression term\n" << endl;

   dComp = pos(alpha1-alphaComp)*(mixture.rhod()/(gComp*(mixture.rhod()-mixture.rhoc())))*mag(UdmModel.Udm())*(lambdaComp/(betaComp+(alpha1-alphaComp)*mixture.rhod()))*dummyComp;
Info<< "Maximum value of compression term:" << max(dComp).value()<< endl;

