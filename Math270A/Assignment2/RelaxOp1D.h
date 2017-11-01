//
//  RelaxOp1D.h
//  Assignment2
//
//  Created by Kevin Miller on 10/17/17.
//  Copyright © 2017 Kevin Miller. All rights reserved.
//
#include <iostream>   // For C++ input/output
#include <vector>     // For STL vectors
#include <cstdio>     // For C input/output
#include <cmath>      // For math functions
#include <cstdlib>
using namespace std; // Using the "standard" (std) standard template library components


#include "GridFun1D.h"
#include "TriSolver.h"


void setUpSystem(double h, double alpha, double beta, long M,
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
        loDiag[i-1] =  alpha/(h*h);
        upDiag[i]   =  alpha/(h*h);
        diag[i]     = -2.0*alpha/(h*h) + beta;
    }
    
    
    // Last equation is the identity
    
    i          =    M;
    diag[M]    =  1.0;
    loDiag[M-1] = 0.0;
}


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


#ifndef _RelaxOp1D_
#define _RelaxOp1D_


class RelaxOp1D {
    
public:
    TriSolver tSolver;
    double dt;
    double alpha;
    GridFun1D f;
    
    // initialize function for RelaxOp1D
    void initialize(double dt, double alpha, const GridFun1D& f){
        this->dt = dt;
        this->alpha = alpha;
        this->f = f;
        
        //cout << f;
        
        TriSolver tSolver; // instantiate the TriSolver
        
        // instantiate vectors
        vector<double> loDiag;
        vector<double> diag;
        vector<double> upDiag;
        
        // Setting up system. Gets us the matrix A's diagonals
        setUpSystem(f.hx, alpha, 0.0, f.xPanel, loDiag, diag, upDiag);
        
        //print_diags(loDiag, diag, upDiag);
        
        // adjust diagonals to satisfy (I - 0.5*dt*A), noting that 1 is already in the first and last entries of diag.
        double val = -0.5*dt;
        for(long i = 1; i < f.xPanel; i++){ //Note only iterate thru one less because loDiag and upDiag are one less in length
            diag[i] *= val;
            diag[i] += 1.0;
            loDiag[i-1] *= val;
            upDiag[i] *= val;
            
        }
        
        //print_diags(loDiag, diag, upDiag);
        
        
        // Initializing the forward operator
        tSolver.initialize(f.xPanel+1,loDiag,diag,upDiag);
        
        this->tSolver = tSolver;
    };
    
    
    // apply function for RelaxOp1D
    void apply(const GridFun1D& uIn, GridFun1D& uOut){
        // forward step calculation (I + 0.5*dt*A)uIn - dt*f
        
        GridFun1D ustar; // make copy of uIn
        ustar = uIn;
        double val = alpha/(uIn.hx * uIn.hx);
        double halfdt = 0.5*dt;
        
        ustar.values[0] = 0.0; // D_D+ operator
        for(long i=1; i < ustar.xPanel; i++){
            ustar.values[i] = uIn.values[i-1] - 2.0*uIn.values[i] + uIn.values[i+1];
            ustar.values[i] *= val*halfdt;
        }
        ustar.values[ustar.xPanel] = 0.0;
        
        ustar += uIn;
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
    };
    
    
};






#endif /* _RelaxOp1D_ */


