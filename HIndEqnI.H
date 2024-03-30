forAll(mesh.cells(), cellI){
        Vector<double> H_PI = Vector<double>::zero; // Initialize for accumulation
        if ((quick == 0) or (sigma[cellI] > 1)){
                forAll(mesh.cells(), cellJ){
                        if (cellI != cellJ){
                                if (sigma[cellJ] > 1){
                                        r = mesh.C()[cellJ] - mesh.C()[cellI]; 
                                        H_PI += (JI[cellJ] ^ r) / pow(mag(r), 3) * mesh.V()[cellJ];
                                }
                        }
                        
                }
        }
        HIndI[cellI] = (H_PI * (-1/(4*3.1415)));
}





