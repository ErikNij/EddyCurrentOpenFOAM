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
    emFoam

Description
    Solver for (quasi-)stationary electromagnetic fields produced by steady or 
    alternating currents with a given frequency. The vector potential is 
    solved in the frequency domain using a block coupled matrix solver so that 
    even high-frequency electromagnetic phenomena can be treated well. 

    Note: "Jcoil" refers to the source current producing the electromagnetic 
    field I kept this naming convention to be consistent with the documentation 
    in my master thesis, where you can read more about the theory behind the 
    solver: 
    Busse, Christian: Numerical Modeling of an Inductively Coupled 
    Plasma (ICP). Ilmenau 2019.  https://doi.org/10.22032/dbt.40314

    In the case setup the following inputs need to be specified
    in "../case/0":
        - Source current density Jcoil in A/m^2
        - Electrical conductivity sigma in A/Vm
    in "../case/constant/physicalProperties":
        - Magnetic permeability muMag (default muMag=mu0) in Vs/Am
        - Current frequency w in 1/s

    The output quantities are:
        - Magnetic vector potential A in Vs/m
        - Magnetic flux density B in Vs/m^2
        - Magnetic field strength H in A/m
        - Induced current density Jind in A/m^2
        - (Time-averaged) Joule heat density in W/m^3  
        - (Time-averaged) Lorentz-force density fL in N/m^3

Author
    Christian Busse

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

Vector<double> r;
Vector<double> J_Q;
Vector<double> H_newI;
Vector<double> H_newR;

float learningRate;
float alpI;
float alpR;
float alpINew;
float alpRNew;

//std::cout << "Learning Rate";
            //std::cin >> learningRate;
//std::cout << "Init scalar";
            //std::cin >> alpI;
//alpI = alpI * muMag.value();
alpR = -alpI;
float sosRNew,sosINew;
float maxImg;
float maxReal;
float maxRealRaw;
float maxImgRaw;
float maxSlope;
float change;
int counter;
float jacStep = 0.1;
std::cout << "Number of outer loops: ";
            std::cin >> counter;
float TRelaxationFactor = 0.2;
std::cout << "TRelaxationFactor: ";
            std::cin >> TRelaxationFactor;
float HRelaxationFactor = 0.1;
float NumJacCutOff = 0.01; // If NumJac of a certian cell is less than NumJacCutOff * maxNumJac, the numJac of that cell will be treated as 0
int innerCounterSave;
int nCells = 0;
int nCellCond = 0;
//std::cout << "Number of inner loops:";
            //std::cin >> innerCounterSave;
int innerCounter;
float maxmax = 2;
float sumMax = 0;
//float numbJacCut = 0.25;
float localMax;
int quick = 1;
labelUList next;
    Info<< "\nStarting time loop\n" << endl;
std::cout << typeid(TI).name() << '\n';
    // Getting edge cells
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        //Determing edge and non edge cells
        forAll(mesh.cells(), cellI) {
            nCells = nCells +1;
            if (sigma[cellI] > 1){
                nCellCond = nCellCond + 1;
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
                RestrictVec[cellI] = RestrictVec[cellI] / mag(RestrictVec[cellI]);
            }
        }

        # include "HCoilEqn.H"
        int i = 0;

        while (i < counter){

            Info<< "Calculated HCoil" << nl;

            HR = HIndR + HCoil;
            HI = HIndI;
            //TNew is really TOld but dont worry about it :)

            TINew = TI;
            TRNew = TR;

            solve(fvm::laplacian(TI) == fvc::grad(fvc::div(TI)) + HR*muMag*w*sigma);
            if (i!= 0){
                solve(fvm::laplacian(TR) == fvc::grad(fvc::div(TR)) - HI*muMag*w*sigma);
            }
            #include "TProj.H"
            
            TI = (TRelaxationFactor * TI) + ((1 - TRelaxationFactor)*TINew);
            TR = (TRelaxationFactor * TR) + ((1 - TRelaxationFactor)*TRNew); 

            JI = fvc::curl(TI);
            JR = fvc::curl(TR);
            # include "HIndEqn.H"
            HR = HIndR + HCoil;
            HI = HIndI;
            i = i + 1;
        }

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
