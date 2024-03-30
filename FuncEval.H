//Let alpINew and alpRNew be set in, will return sosI and R

TINew = TI + (alpINew * HR);
TRNew = TR + (alpRNew * HI); 


projectionI = HI & RestrictVec; //no need to divide by mag(RestrictVec) because it has a mag of 1
projectionR = HR & RestrictVec;

forAll(mesh.cells(), cellI) {
    if (EdgeCell[cellI] == 1){
        TINew[cellI] = TI[cellI] - (alpINew * HR[cellI]) + (alpINew * projectionR[cellI] * RestrictVec[cellI]);
        TRNew[cellI] = TR[cellI] - (alpRNew * HI[cellI]) + (alpRNew * projectionI[cellI] * RestrictVec[cellI]);
    }
}

//TI and TINew are now correct.
dBIdt = HR * muMag * w;
dBRdt = - HI * muMag * w;

dBIndIdt = fvc::curl(1/sigma * fvc::curl(TINew));
dBIndRdt = fvc::curl(1/sigma * fvc::curl(TRNew));

ImgResNew = dBIndIdt - dBIdt;
RealResNew = dBIndRdt - dBRdt;

sosINew = 0;
sosRNew = 0;

forAll(mesh.cells(), cellI) {
    if (sigma[cellI] > 1){

        sosINew =+ mag(ImgResNew[cellI]);
        sosRNew =+ mag(RealResNew[cellI]);
    }
}

Info<< typeIter << " | " << alpINew << " | " << sosINew << " I"<<nl;
Info<< typeIter << " | " << alpRNew << " | " << sosRNew << " R"<<nl;
