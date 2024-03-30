/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    MyEmFoam3

Description
    Solver for (quasi-)stationary electromagnetic fields produced by 
    alternating currents with a given frequency. The potential (T) is 
    solved in the frequency domain using  matrix solver.

    Note: "JCoil" refers to the source current producing the electromagnetic 
    field I kept this naming convention to be consistent with the documentation 
    in Christian Busse's  master thesis, where you can read more about the theory behind the 
    solver: 
    Busse, Christian: Numerical Modeling of an Inductively Coupled 
    Plasma (ICP). Ilmenau 2019.  https://doi.org/10.22032/dbt.40314

    In the case setup the following inputs need to be specified
    in "../case/0":
        - Source current density JCoil in A/m^2
        - Electrical conductivity sigma in A/Vm
        - Intial guess for potential (T)
    in "../case/constant/physicalProperties":
        - Magnetic permeability muMag (default muMag=mu0) in Vs/Am
        - Current frequency w in 1/s

    The output quantities are:
        - Magnetic vector potential T in 
        - Magnetic flux density H in 
        - Magnetic field strength B in A/m
        - Induced current density JI and JR in A/m^2
        Eventully
        - (Time-averaged) Joule heat density in W/m^3  
        - (Time-averaged) Lorentz-force density fL in N/m^3

Author
    Erik Nijkamp

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initConvergenceCheck.H"
#   include<list>

//Vector<double> r;
//Vector<double> J_Q;
Vector<double> H_newI;
Vector<double> H_newR;
int counter;
float TRelaxationFactor;
//std::cout << "TRelaxationFactor: ";
            //std::cin >> TRelaxationFactor;
float largeNumber;
int quick;
float DHI;
float DTI;
float DHINorm;

std::ifstream file("system/SolverSettings.csv");
std::string line;
if (file.is_open()) {
    while (std::getline(file, line)) {
        std::stringstream lineStream(line);
        std::string field, value;
        std::getline(lineStream, field, ',');
        std::getline(lineStream, value, ',');
        if (field == "quick") {
            quick = std::stoi(value);
        } else if (field == "largeNumber") {
            largeNumber = std::stod(value);
        } else if (field == "outerLoops") {
            counter = std::stoi(value);
        } else if(field == "TRelaxationFactor"){
            TRelaxationFactor = std::stod(value);
        }
    }
    file.close();
} else {
    std::cerr << "Error opening file" << std::endl;
    return 1; // Indicate failure
}

//volVectorField HIOld;
//volVectorField HROld;
labelUList next;

//Determing edge and non edge cells
        forAll(mesh.cells(), cellI) {
            if (sigma[cellI] > 1){
                const cell& c = mesh.cells()[cellI];
                const labelList& neighbors = mesh.cellCells()[cellI];
                forAll(neighbors, i) {
                    label neighborCellIndex = neighbors[i];
                    if (sigma[neighbors[i]] < 1){
                            EdgeCell[cellI] = 1;
                            RestrictVec[cellI] = RestrictVec[cellI] + mesh.C()[neighbors[i]] - mesh.C()[cellI];
                    }
                } 
            }
        }
        forAll(mesh.cells(), cellI) {
            if (mag(RestrictVec[cellI]) > 0){
                RestrictVec[cellI] = largeNumber * RestrictVec[cellI] / mag(RestrictVec[cellI]);
            }
            //else{RestrictVec[cellI] = Foam::Vector<double>(1.0, 1.0, 1.0);}
            
        }

        # include "HCoilEqn.H"
        int i = 0;
        Info<< "Calculated HCoil" << nl;

    Info<< "\nStarting time loop\n" << endl;
std::cout << typeid(TI).name() << '\n';
    // Getting edge cells
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        int i = 0;
        
        while (i < counter){


            HR = HIndR + HCoil;
            HI = HIndI;
            //TNew is really TOld but dont worry about it :)
            //TI = ((TI & (TI + RestrictVec)) /pow(mag(TI + RestrictVec),2)) * (TI + RestrictVec)
            //TR = ((TR & (TR + RestrictVec)) /pow(mag(TR + RestrictVec),2)) * (TR + RestrictVec)


            //TINew = TI;
            //TRNew = TR;

            forAll(mesh.cells(), cellI) {
                TI[cellI] = Foam::Vector<double>(0.0, 0.0, 0.0);
                TR[cellI] = Foam::Vector<double>(0.0, 0.0, 0.0);
            }
            solve(relax(fvm::laplacian(TI)) == 
            fvc::grad(fvc::div(((TI & (TI + RestrictVec)) /pow(mag(TI + RestrictVec)+1e-10,2)) * (TI + RestrictVec))) 
            + HR*muMag*w*sigma 
            + (sigma * (fvc::grad(1/sigma) ^ fvc::curl(((TI & (TI + RestrictVec)) /pow(mag(TI + RestrictVec)+1e-10,2)) * (TI + RestrictVec)))));
            //+ gamma*(TI - (mag(TI) * RestrictVec * (RestrictVec & TI))) );
            
            TI = ((TI & (TI + RestrictVec)) /pow(mag(TI + RestrictVec)+1e-10,2)) * (TI + RestrictVec);
            DTI = sum(mag(TI-TINew)).value();
            TI = (TRelaxationFactor * TI) + ((1 - TRelaxationFactor)*TINew);
            JI = fvc::curl(TI);

            # include "HIndEqnI.H"
            DHI = sum(mag((HI - HIndI))).value();
            DHINorm = sum(mag((HI - HIndI))).value() / sum(mag(HIndI)).value();
            HI = HIndI;
            solve(relax(fvm::laplacian(TR)) ==
            fvc::grad(fvc::div(((TR & (TR + RestrictVec)) /pow(mag(TR + RestrictVec)+1e-10,2)) * (TR + RestrictVec))) 
            - HI*muMag*w*sigma 
            + (sigma * (fvc::grad(1/sigma) ^ fvc::curl(((TR & (TR + RestrictVec)) /pow(mag(TR + RestrictVec)+1e-10,2)) * (TR + RestrictVec)))));
            //+ gamma* (TR - (mag(TR) * RestrictVec * (RestrictVec & TR))));

            

            //#include "TProj.H"
            //TI = ((TI & (TI + RestrictVec)) /pow(mag(TI + RestrictVec)+1e-10,2)) * (TI + RestrictVec);
            TR = ((TR & (TR + RestrictVec)) /pow(mag(TR + RestrictVec)+1e-10,2)) * (TR + RestrictVec);

            //DTI = sum(mag(TI-TINew)).value();

            TR = (TRelaxationFactor * TR) + ((1 - TRelaxationFactor)*TRNew); 

            //JI = fvc::curl(TI);
            JR = fvc::curl(TR);
            divJI = fvc::div(JI);
            # include "HIndEqnR.H"
            //# include "HIndEqn.H"
            Info << "norm Delta HI:" << DHINorm << nl;
            Info << "Delta HI:" << DHI << nl;
            Info << "Delta TI:" << DTI << nl;
            HR = HIndR + HCoil;
            HI = HIndI;
            i = i + 1;
        }
        
        Hmag = mag(HR);
        Bmag = Hmag * muMag;

        // Compute the time-averaged Joule heat (power density) and Lorentz-force 
        //qJ = sigma/2.0 * sqr(w) * sqr(Amag); // Dissipated power density in [W/m^3]   
        //fL = 0.5 * ((JR ^ BR) + (JI ^ BI)); // Lorentz-force in [N/m^3]

        // Compute the total dissipated power in [W]
        //dimensionedScalar QJ = fvc::domainIntegrate(qJ);       

        //Info<< "----- Total dissipated power: " << nl << endl;
        //Info<< "QJ  = " << QJ.value() << nl << endl;

        //DHCoil = - muMag * w * HCoil;

        Info<< "End\n" << endl;

        Info<< "EM Solver: ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
           << nl << endl;

        #include "convergenceCheck.H"

        runTime.write();
    }


    return 0;
}
