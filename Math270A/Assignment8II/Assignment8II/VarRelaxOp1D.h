//
//  VarRelaxOp1D.h
//  Assignment8II
//
//  Created by Kevin Miller on 11/29/17.
//  Copyright © 2017 Kevin Miller. All rights reserved.
//
#include <cstdio>
#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

#include "GridFun1D.h"
#include "TriSolver.h"

// Set Up System for the variable coefficient function option
void setUpSystemVar(double h, vector<double> Alpha, long M,
                    vector<double>& loDiag, vector<double>& diag, vector<double>& upDiag)
{
    loDiag.resize(M);
    upDiag.resize(M);
    diag.resize(M+1);
    
    long i;
    
    // First equation is the identity
    
    i = 0;
    diag[0]   =  1.0;
    upDiag[0] =  0.0;
    
    // Interior equation coefficients associated with a standard
    // second order finite difference approximation
    
    for(long i = 1; i < M; i++)
    {
        loDiag[i-1] =  (Alpha[i] + Alpha[i-1])/(2.0*h*h);
        diag[i] =  -(0.5*Alpha[i+1] + Alpha[i] + 0.5*Alpha[i-1])/(h*h);
        upDiag[i] = (Alpha[i+1] + Alpha[i])/(2.0*h*h);
    }
    
    
    // Last equation is the identity
    i          =    M;
    diag[M]    =  1.0;
    loDiag[M-1] = 0.0;
}


#ifndef _VarRelaxOp1D_
#define _VarRelaxOp1D_

class VarRelaxOp1D {
    
public:
    TriSolver tSolver;
    double dt;
    GridFun1D alpha;
    GridFun1D f;
    
    // initialize function for VarRelaxOp2D
    void initialize(double dt, const GridFun1D& alpha, const GridFun1D& f){
        this->dt = dt;
        this->alpha = alpha;
        this->f = f;
        
        // instantiate the TriSolver member variable
        TriSolver tSolver;
        
        long M = f.xPanel;
        
        // instantiate vectors, picking out the values of alpha for setting up the system
        vector<double> loDiag;
        vector<double> diag;
        vector<double> upDiag;
        vector<double> Alpha;
        Alpha = alpha.values;
        
        // Set up the tridiag system's matrix A, in Au = f
        setUpSystemVar(f.hx, Alpha, M, loDiag, diag, upDiag);
        
        // now adjust the diagonals to satisfy (I - 0.5*dt*A), ie the Crank-Nicholson step
        // note that the identity remains unchanged on the boundary
        double val = -0.5*dt;
        for(long i=1; i < M; i ++){
            diag[i] *= val;
            diag[i] += 1;
            loDiag[i-1] *= val;
            upDiag[i] *= val;
        }
        
        // initialize the forward operator
        tSolver.initialize(M+1, loDiag, diag, upDiag);
        this->tSolver = tSolver;
        
        
        return;
    };
    
    void apply(const GridFun1D& uIn, GridFun1D& uOut){
        // foward step calculation (I + 0.5*dt*A)uIn - dt*f
        GridFun1D ustar; // make copy of uIn
        ustar = uIn;
        double val = 0.5*dt/(uIn.hx*uIn.hx);
        
        ustar.values[0] = 0.0; // D_D+ operator, with the variable coeff (0.5*dt*A)*uIn
        for(long i=1; i < ustar.xPanel; i++){
            ustar.values[i] = 0.5*(alpha.values[i-1] + alpha.values[i])*uIn.values[i-1] + -(0.5*alpha.values[i+1] + alpha.values[i] + 0.5*alpha.values[i-1])*uIn.values[i] + 0.5*(alpha.values[i+1] + alpha.values[i])*uIn.values[i+1];
            ustar.values[i] *= val;
        }
        ustar.values[ustar.xPanel] = 0.0; // boundary is 0, the previous operator only applies to interior pts
        
        // add I*uIn
        ustar += uIn;
        
        // subtract dt*f
        GridFun1D fstar;
        fstar = f;
        fstar *= dt;
        ustar -= fstar;
        
        // solve (I - 0.5*dt*A)uOut = ustar
        vector<double> ustar_vector = ustar.values;
        vector<double> uOut_vector = uOut.values;
        tSolver.apply(f.xPanel+1, ustar_vector, uOut_vector);
        
        for(long i=0; i < ustar.xPanel+1; i++){
            uOut.values[i] = uOut_vector[i];
        }
        
        return;
        
    }
};


#endif /* _VarRelaxOp1D_ */

