    // Solve for the H feild created by the coil 
Vector<double> r;
Vector<double> J_Q;
Vector<double> H_new;
forAll(mesh.cells(), cellI){
        Vector<double> H_P = Vector<double>::zero;
        if ((quick == 0) or (sigma[cellI] > 1)){
                forAll(mesh.cells(), cellJ){
                        if (cellI != cellJ){
                                if (mag(JCoil[cellJ]) > 0){
                                        r = mesh.C()[cellJ] - mesh.C()[cellI]; 
                                        J_Q = JCoil[cellJ];
                                        H_new = ((J_Q ^ r) * pow(mag(r), -3) * mesh.V()[cellJ]);
                                        H_P = H_P + H_new;
                                }
                        }
                        
                }
        }
        HCoil[cellI] = H_P * (-1/(4*3.1415));
}
    





