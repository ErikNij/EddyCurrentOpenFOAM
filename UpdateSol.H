//Let alpINew and alpRNew be set in, will return sosI and R

TI = TI + (alpINew * HI);
TR = TR + (alpRNew * HR); 


/* projectionI = HI & RestrictVec; //no need to divide by mag(RestrictVec) because it has a mag of 0
projectionR = HR & RestrictVec;

forAll(mesh.cells(), cellI) {
    if (EdgeCell[cellI] == 1){
        TI[cellI] = TI[cellI] - (alpINew * HI[cellI]) +  (alpI * projectionI[cellI] * RestrictVec[cellI]);
        TR[cellI] = TR[cellI] - (alpRNew * HR[cellI]) +  (alpR * projectionR[cellI] * RestrictVec[cellI]);

        TINew[cellI] = TI[cellI] + (alpINew * projectionI[cellI] * RestrictVec[cellI]);
        TRNew[cellI] = TR[cellI] + (alpRNew * projectionR[cellI] * RestrictVec[cellI]);
    }
}
*/
//TI and TINew are now correct.
ImgRes = fvc::curl(1/sigma * fvc::curl(TI)) - HI;
RealRes = fvc::curl(1/sigma * fvc::curl(TR)) - HR;


Info<< alpINew << " | " << sosINew << nl;
