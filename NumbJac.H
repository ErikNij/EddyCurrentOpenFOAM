//Evalute the numerical jacobian? 
        
//X dimetion
//(re)setting the potential 

TI = TI + (alpI * HI);
TR = TR + (alpR * HR);

alpINew = alpI * 0.1;
alpRNew = alpR * 0.1;

TINew = TI + (alpINew * HI);
TRNew = TR + (alpRNew * HR); 

projectionI = HI & RestrictVec; //no need to divide by mag(RestrictVec) because it has a mag of 0
projectionR = HR & RestrictVec;

forAll(mesh.cells(), cellI) {
    if (EdgeCell[cellI] == 1){
        TI[cellI] = TI[cellI] - (alpINew * HI[cellI]) +  (alpI * projectionI[cellI] * RestrictVec[cellI]);
        TR[cellI] = TR[cellI] - (alpRNew * HR[cellI]) +  (alpR * projectionR[cellI] * RestrictVec[cellI]);

        TINew[cellI] = TI[cellI] + (alpINew * projectionI[cellI] * RestrictVec[cellI]);
        TRNew[cellI] = TR[cellI] + (alpRNew * projectionR[cellI] * RestrictVec[cellI]);
    }
}
//TI and TINew are now correct.
ImgResNew = fvc::curl(1/sigma * fvc::curl(TINew)) - HI;
RealResNew = fvc::curl(1/sigma * fvc::curl(TRNew)) - HR;

ImgRes = fvc::curl(1/sigma * fvc::curl(TI)) - HI;
RealRes = fvc::curl(1/sigma * fvc::curl(TR)) - HR;
float sosI = 0;
float sosR = 0;

float sosINew = 0;
float sosRNew = 0;
forAll(mesh.cells(), cellI) {
    if (sigma[cellI] > 1){
        sosI =+ mag(ImgRes[cellI]);
        sosR =+ mag(RealRes[cellI]);

        sosINew =+ mag(ImgResNew[cellI]);
        sosRNew =+ mag(RealResNew[cellI]);
    }
}


float dFitnessI = sosINew - sosI; 
float dFitnessR = sosRNew - sosR;

float NumJacReal = dFitnessR / (alpR);
float NumJacImg = dFitnessI / (alpI);

Info<< alpI << " | " << sosI << " | " << NumJacImg << nl;




