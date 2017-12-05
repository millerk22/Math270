//
//  VarRelaxOp2D.h
//  Assignment9
//
//  Created by Kevin Miller on 12/2/17.
//  Copyright © 2017 Kevin Miller. All rights reserved.
//
#include <iostream>   // For C++ input/output
#include <vector>     // For STL vectors
#include <cstdio>     // For C input/output
#include <cmath>      // For math functions
#include <cstdlib>
using namespace std; // Using the "standard" (std) standard template library components


#include "GridFun2D.h"
#include "TriSolver.h"
#include "GNUplotUtility.h"

void print_diags(vector<double> loDiag, vector<double> diag, vector<double> upDiag){
    long n = diag.size();
    cout << loDiag.size() << " " << n << " " << upDiag.size() << endl;
    for(long i =0; i < n-1; i++){
        if(i == 0){
            cout << "\t\t" << setw(5) << diag[i] << "\t" << setw(5) << upDiag[i] << endl;
        }
        else{
            cout << setw(5) << loDiag[i-1] << "\t" << setw(5) << diag[i] << "\t" << setw(5) << upDiag[i] << endl;
        }
    }
    cout << setw(5) << loDiag[n-2] << "\t" << setw(5) << diag[n-1] << endl;
    
}

#ifndef _VarRelaxOp2D_
#define _VarRelaxOp2D_


class VarRelaxOp2D {
    
public:
    vector<TriSolver> triSolverXX;
    vector<TriSolver> triSolverYY;
    double dt;
    GridFun2D alpha;
    GridFun2D F;
    
    // initialize function for RelaxOp2D
    void initialize(double dt, const GridFun2D& alpha, const GridFun2D& f){
        this->dt = dt;
        this->alpha = alpha;
        this->F = f;
        F *= -dt; // Already apply the dt multiplication to the constant F, so only calculate once
        
        vector<TriSolver> triSolverXX; // vectors to store all the different TriSolvers
        vector<TriSolver> triSolverYY;
        triSolverXX.resize(alpha.xPanel-1);
        triSolverYY.resize(alpha.yPanel-1);
        
        // instantiate vectors, to be used in the different TriSolvers
        vector<double> loDiagX;
        vector<double> diagX;
        vector<double> upDiagX;
        
        vector<double> loDiagY;
        vector<double> diagY;
        vector<double> upDiagY;
        
        // Resize the vectors
        long Mx = F.xPanel + 1;
        long My = F.yPanel + 1;
        loDiagX.resize(Mx-1);
        diagX.resize(Mx);
        upDiagX.resize(Mx-1);
        loDiagY.resize(My-1);
        diagY.resize(My);
        upDiagY.resize(My-1);
        
        
        // Instantiate the X and Y TriSolvers with different alpha values, as specified in Douglas ADI scheme
        
        double Xval = 0.5*dt/(alpha.hx*alpha.hx);
        diagX[0] = 1.0; // set the boundary values to remain constant
        upDiagX[0] = 0.0;
        for(long j = 1; j < F.yPanel; j++){
            for(long i = 1; i < F.yPanel; i++){
                diagX[i] = Xval*(0.5*alpha.values(i+1,j) + alpha.values(i,j) + 0.5*alpha.values(i-1,j));
                diagX[i] += 1.0;
                loDiagX[i-1] = -Xval*(0.5*alpha.values(i,j) + 0.5*alpha.values(i-1,j));
                upDiagX[i] = -Xval*(0.5*alpha.values(i+1,j) + 0.5*alpha.values(i,j));
            }
            loDiagX[Mx-2] = 0.0; // set the boundary values to remain constant
            diagX[Mx-1] = 1.0;
            triSolverXX[j-1].initialize(Mx, loDiagX, diagX, upDiagX);
            
            
            
            
        }
        
        double Yval = 0.5*dt/(alpha.hy*alpha.hy);
        diagY[0] = 1.0;
        upDiagY[0] = 0.0;
        for(long i = 1; i < F.yPanel; i++){
            for(long j = 1; j < F.yPanel; j++){
                diagY[j] = Yval*(0.5*alpha.values(i,j+1) + alpha.values(i,j) + 0.5*alpha.values(i,j-1));
                diagY[j] += 1.0;
                loDiagY[j-1] = -Yval*(0.5*alpha.values(i,j) + 0.5*alpha.values(i,j-1));
                upDiagY[j] = -Yval*(0.5*alpha.values(i,j+1) + 0.5*alpha.values(i,j));
            }
            loDiagY[My-2] = 0.0;
            diagY[My-1] = 1.0;
            triSolverYY[i-1].initialize(My, loDiagY, diagY, upDiagY);
            
            
        }
        
        
        this->triSolverXX = triSolverXX;
        this->triSolverYY = triSolverYY;
        
    };
    
    
    // apply function for RelaxOp2D
    void apply(const GridFun2D& uIn, GridFun2D& uOut){
        long Mx = uIn.xPanel + 1;
        long My = uIn.yPanel + 1;
        uOut = uIn; // uOut gets the right boundary values set right away
        
        
        GridFun2D ustar; // make copy of uIn, use as intermediate step
        ustar = uIn; // copy initializer to get same sizes and place to start
        
        
        // Step 1. This, with ustar = uIn to start, gives (I + 0.5*dt*D_D+(X) + dt*D_D+(Y)
        double valx = 0.5*dt/(uIn.hx*uIn.hx);
        double valy = dt/(uIn.hy*uIn.hy);
        double alphax_1 = 0.0; // instantiate the alpha variables to make the coding in the double for loop easier
        double alphax = 0.0;
        double alphaxp1 = 0.0;
        double alphay_1 = 0.0;
        double alphay = 0.0;
        double alphayp1 = 0.0;
        
        
        for(long i=1; i < Mx - 1; i++){
            for(long j = 1; j < My - 1; j++){
                alphax_1 = 0.5*valx*(alpha.values(i-1,j)+ alpha.values(i,j));
                alphax = -valx*(0.5*alpha.values(i+1,j) + alpha.values(i,j) + 0.5*alpha.values(i-1,j));
                alphaxp1 = 0.5*valx*(alpha.values(i+1,j) + alpha.values(i,j));
                ustar.values(i,j) += alphax_1*uIn.values(i-1,j) + alphax*uIn.values(i,j) + alphaxp1*uIn.values(i+1,j); // D_D+ X operator on uIn
                
                alphay_1 = 0.5*valy*(alpha.values(i,j-1) + alpha.values(i,j));
                alphay = -valy*(0.5*alpha.values(i,j+1) + alpha.values(i,j) + 0.5*alpha.values(i,j-1));
                alphayp1 = 0.5*valy*(alpha.values(i,j+1) + alpha.values(i,j));
                ustar.values(i,j) += alphay_1*uIn.values(i,j-1) + alphay*uIn.values(i,j) + alphayp1*uIn.values(i,j+1); // D_D+ Y operator on uIn
            }
        }
        
        
        ustar += F; // note that F is really -dt*F. This finishes the forward step by subtracting -dt*F from ustar
        
        
        GridFun2D ustar2;
        ustar2 = uIn;
        
        // solve (I - 0.5*dt*alphaX*D_D+ X)ustar2 = ustar for each row (for each j)
        
        vector<double> ustar_j;
        vector<double> ustar2_j;
        
        ustar_j.resize(Mx);
        ustar2_j.resize(Mx);
        
        for (long j = 1; j < My-1; j++){
            for (long i=0; i < Mx; i++){
                ustar_j[i] =  ustar.values(i,j);
            }
            
            triSolverXX[j-1].apply(Mx, ustar_j, ustar2_j);
            
            for (long i=0; i < Mx; i++){
                ustar2.values(i,j) = ustar2_j[i];
            }
        }
        
        
        
        
        // Step 2 - Subtract 0.5*dt*alpha(i,j)*D_D+(Y) on uIn from previous calculation.
        
        for(long i = 1; i < Mx-1; i++){
            for(long j = 1; j < My-1; j++){
                alphay_1 = 0.25*valy*(alpha.values(i,j-1) + alpha.values(i,j));
                alphay = -0.5*valy*(0.5*alpha.values(i,j+1) + alpha.values(i,j) + 0.5*alpha.values(i,j-1));
                alphayp1 = 0.25*valy*(alpha.values(i,j+1) + alpha.values(i,j));
                ustar2.values(i,j) -= alphay_1*uIn.values(i,j-1) + alphay*uIn.values(i,j) + alphayp1*uIn.values(i,j+1); // subtract 0.5*dt*D_D+ Y operator on uIn from ustar2
            }
        }
        
        
        
        // Solve (I - 0.5*dt*alphaY*D_D+ Y uOut = ustar2
        vector<double> ustar2_i;
        vector<double> uOut_i;
        ustar2_i.resize(My);
        uOut_i.resize(My);
        
        for (long i = 1; i < Mx-1; i++){
            for (long j=0; j < My; j++){
                ustar2_i[j] =  ustar2.values(i,j);
            }
            
            triSolverYY[i-1].apply(My, ustar2_i, uOut_i);
            for (long j=0; j < My; j++){
                uOut.values(i,j) = uOut_i[j];
            }
        }
        
        
        
        return;
    };
    
    
};



#endif /* _VarRelaxOp2D_ */
