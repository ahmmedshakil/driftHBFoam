/*-solving the transport equation for the tracer test: added by shakil
-Ref. Simick et al., Smart RTD for multiphase flow systems, Chemical Eng. Rese
-and Design 90 (2012), 1739-1749. 
Sc = Schmdit number = 0.7
Dm = molecular diffusion of the tracer 
DT = turbulent diffusion + molecular diffusion
*/
Info<< "Solving transport equation C\n" << endl;

 {
   //-Calcuate the diffusivity of the tracer
    DT = turbulence->nut() / Sc + Dm;   	//-

    fvScalarMatrix CEqn
    (
        fvm::ddt(C)                         //-C[0,0,0,0,0,0,0]
        + fvm::div(phi, C)
        - fvm::laplacian(DT, C)             //-DT [kg m s K mol A cd] = [0,0,-1,0,0,0,0]
        ==
        fvOptions(C)
    );

    CEqn.relax();
    fvOptions.constrain(CEqn);
    CEqn.solve();
    fvOptions.correct(C);
}
