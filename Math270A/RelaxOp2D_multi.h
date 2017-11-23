//
//  RelaxOp2D.h
//  Assignment4
//
//  Created by Kevin Miller on 10/30/17.
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

#ifdef _OPENMP
#include<omp.h>
#endif

#ifndef _RelaxOp2Dm_
#define _RelaxOp2Dm_



class RelaxOp2Dm {
    
public:
    TriSolver triSolverX;
    TriSolver triSolverY;
    double dt;
    double alphaX;
    double alphaY;
    GridFun2D F;
    
    // initialize function for RelaxOp2D
    void initialize(double dt, double alphaX, double alphaY, const GridFun2D& f){
        this->dt = dt;
        this->alphaX = alphaX;
        this->alphaY = alphaY;
        this->F = f;
        F *= -dt; // Already apply the dt multiplication to the constant F, so only calculate once
        
        
        TriSolver triSolverX; // instantiate the X direction TriSolver
        TriSolver triSolverY; // instantiate the Y direction TriSolver
        
        // instantiate vectors
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
        
        double valx = -0.5*dt*alphaX/(F.hx*F.hx);
        double valy = -0.5*dt*alphaY/(F.hy*F.hy);
        
        diagX[0] = 1.0; // set the boundary values to remain constant
        diagY[0] = 1.0;
        upDiagX[0] = 0.0;
        upDiagY[0] = 0.0;
        
        // Fill in the diagonal vectors values for the X and Y Trisolvers, as specified in Douglas ADI scheme
        for(long i = 1; i < F.xPanel; i++){
            for(long j = 1; j < F.yPanel; j++){
                diagX[i] = -2.0*valx;
                diagX[i] += 1.0;
                loDiagX[i-1] = valx;
                upDiagX[i] = valx;
                diagY[i] = -2.0*valy;
                diagY[i] += 1.0;
                loDiagY[i-1] = valy;
                upDiagY[i] = valy;
            }
        }
        
        loDiagX[Mx-2] = 0.0; // set the boundary vlaues to remain constant
        loDiagY[My-2] = 0.0;
        diagX[Mx-1] = 1.0;
        diagY[My-1] = 1.0;
        
        
        
        // Initializing the forward operator
        triSolverX.initialize(Mx,loDiagX,diagX,upDiagX);
        triSolverY.initialize(My,loDiagY,diagY,upDiagY);
        
        this->triSolverX = triSolverX;
        this->triSolverY = triSolverY;
    };
    
#ifdef _OPENMP
    // apply function for RelaxOp2D, with OPENMP enabled
    void apply(const GridFun2D& uIn, GridFun2D& uOut){
        long Mx = uIn.xPanel + 1;
        long My = uIn.yPanel + 1;
        uOut = uIn; // uOut gets the right boundary values set right away
        
        
        GridFun2D ustar; // make copy of uIn, use as intermediate step
        ustar = uIn; // copy initializer to get same sizes and place to start
       
	// Forward step
        double Xval = 0.5*alphaX*dt/(uIn.hx*uIn.hx); // Note the constant in front of D_D+ X is
        double Yval = alphaY*dt/(uIn.hy*uIn.hy); // different than D_D+ Y....
        for(long i=1; i < Mx - 1; i++){
            long j;
#pragma omp parallel for default(shared) private(j) schedule(static)
            for(j = 1; j < My - 1; j++){
                ustar.values(i,j) += Xval*(uIn.values(i-1,j) - 2.0*uIn.values(i,j) + uIn.values(i+1,j)); // D_D+ X operator
                ustar.values(i,j) += Yval*(uIn.values(i,j-1) - 2.0*uIn.values(i,j) + uIn.values(i,j+1)); // D_D+ Y operator
            }
        }
        
        ustar += F;
        
        
        
        GridFun2D ustar2;
        ustar2 = uIn;
        
        // solve (I - 0.5*dt*alphaX*D_D+ X)ustar2 = ustar for each row (for each j)
        vector<double> ustar_j;
        vector<double> ustar2_j;
        
        ustar_j.resize(Mx);
        ustar2_j.resize(Mx);
        
        
        TriSolver TSX = triSolverX;
        long j;
#pragma omp parallel for default(shared) private(j) firstprivate(TSX, ustar_j, ustar2_j) schedule(static)
        for (j = 1; j < My-1; j++){
            for (long i=0; i < Mx; i++){
                ustar_j[i] =  ustar.values(i,j);
            }
            
            TSX.apply(Mx, ustar_j, ustar2_j);
            
            for (long i=0; i < Mx; i++){
                ustar2.values(i,j) = ustar2_j[i];
            }
        }
        
        
        for(long i = 1; i < Mx-1; i++){
            long j;
#pragma omp parallel for default(shared) private(j) schedule(static)
            for(j = 1; j < My-1; j++){
                ustar2.values(i,j) -= 0.5*Yval*(uIn.values(i,j-1) - 2.0*uIn.values(i,j) + uIn.values(i, j+1));
            }
        }
        
        // Solve (I - 0.5*dt*alphaY*D_D+ Y uOut = ustar2
        vector<double> ustar2_i;
        vector<double> uOut_i;
        ustar2_i.resize(My);
        uOut_i.resize(My);
        
        long i;
        
#pragma omp parallel for default(shared) private(i) firstprivate(triSolverY, ustar2_i, uOut_i) schedule(static)
        for (i = 1; i < Mx-1; i++){
            for (long j=0; j < My; j++){
                ustar2_i[j] =  ustar2.values(i,j);
            }
            
            triSolverY.apply(My, ustar2_i, uOut_i);
            
            for (long j=0; j < My; j++){
                uOut.values(i,j) = uOut_i[j];
            }
        }
        
        
        return;
    };
    
    
#else
    // apply function for RelaxOp2D WITHOUT OPENMP
    void apply(const GridFun2D& uIn, GridFun2D& uOut){
        long Mx = uIn.xPanel + 1;
        long My = uIn.yPanel + 1;
        uOut = uIn; // uOut gets the right boundary values set right away
        
        
        GridFun2D ustar; // make copy of uIn, use as intermediate step
        ustar = uIn; // copy initializer to get same sizes and place to start
        
        
        // Forward step
        double Xval = 0.5*alphaX*dt/(uIn.hx*uIn.hx); // Note the constant in front of D_D+ X is
        double Yval = alphaY*dt/(uIn.hy*uIn.hy); // different than D_D+ Y....
        for(long i=1; i < Mx - 1; i++){
            for(long j = 1; j < My - 1; j++){
                ustar.values(i,j) += Xval*(uIn.values(i-1,j) - 2.0*uIn.values(i,j) + uIn.values(i+1,j)); // D_D+ X operator
                ustar.values(i,j) += Yval*(uIn.values(i,j-1) - 2.0*uIn.values(i,j) + uIn.values(i,j+1)); // D_D+ Y operator
            }
        }
        
        
        ustar += F;
        
        
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
            
            triSolverX.apply(Mx, ustar_j, ustar2_j);
            
            for (long i=0; i < Mx; i++){
                ustar2.values(i,j) = ustar2_j[i];
            }
        }
        
        
        
        for(long i = 1; i < Mx-1; i++){
            for(long j = 1; j < My-1; j++){
                ustar2.values(i,j) -= 0.5*Yval*(uIn.values(i,j-1) - 2.0*uIn.values(i,j) + uIn.values(i, j+1));
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
            
            triSolverY.apply(My, ustar2_i, uOut_i);
            for (long j=0; j < My; j++){
                uOut.values(i,j) = uOut_i[j];
            }
        }
        
        
        return;
    };
#endif /* _OPENMP */
    
};



#endif /* RelaxOp2Dm_ */
