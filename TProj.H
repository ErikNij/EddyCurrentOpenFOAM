projectionI = TI & RestrictVec;
projectionR = TR & RestrictVec;
forAll(mesh.cells(), cellI) {
        if (EdgeCell[cellI] == 1){
            TI[cellI] = RestrictVec[cellI] * projectionI[cellI];
            TR[cellI] = RestrictVec[cellI] * projectionR[cellI];
        }
}