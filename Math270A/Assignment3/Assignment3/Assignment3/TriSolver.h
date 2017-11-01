//
//  TriSolver.h
//  Assignment2
//
//  Created by Kevin Miller on 10/17/17.
//  Copyright © 2017 Kevin Miller. All rights reserved.
//
#include <iostream>   // For C++ input/output
#include <vector>     // For STL vectors
#include <cstdio>     // For C input/output
#include <cmath>      // For math functions
using namespace std; // Using the "standard" (std) standard template library components


#ifndef _TriSolver_
#define _TriSolver_


class TriSolver
{
public:
    
    // Null Constructor
    TriSolver()
    { initialize();};
    
    TriSolver(const TriSolver& T) // called when you declare an instance with an existing instance
    {initialize(T);};
    
    TriSolver(long systemSize, vector<double>& loDiag, vector<double>& diag, vector<double>& upDiag){
        initialize(systemSize, loDiag, diag, upDiag);
    };
    
    // Null initializer. Set the
    
    void initialize(){
        systemSize = 0;
        L_loDiag.clear();
        L_diag.clear();
        U_upDiag.clear();
    };
    
    ///  Copy initializer. Duplicates the entries of T
    
    void initialize(const TriSolver& T)
    {
        systemSize = T.systemSize;
        L_loDiag = T.L_loDiag;
        L_diag = T.L_diag;
        U_upDiag = T.U_upDiag;
    };
    
    
    // Constructor -- Here we solve for the diagonals of the LU factorization. Thus, when we apply this solver for different boundary conditions (different f vectors), we don't have to refactor into LU factorization.
    
    void initialize(long systemSize, vector<double>& loDiag, vector<double>& diag, vector<double>& upDiag)
    {
        this->systemSize = systemSize;
        
        // Compute LU factorization, as discussed in class
        
        this->L_loDiag = loDiag; // the subdiagonal of L is just the subdiagonal of the original matrix
        
        L_diag.resize(systemSize);
        U_upDiag.resize(systemSize-1);
        
        L_diag[0] = diag[0];
        
        U_upDiag[0] = upDiag[0]/L_diag[0];
        
        for (long i=1; i <= systemSize-2; i ++){
            L_diag[i] = diag[i] - loDiag[i-1]*U_upDiag[i-1]; // compute main diagonal of L
            U_upDiag[i] = upDiag[i]/L_diag[i]; // compute the super diagonal of U
        }
        
        L_diag[systemSize-1] = diag[systemSize-1] - loDiag[systemSize-2]*U_upDiag[systemSize-2]; // compute last entry in the main diagonal of L
        
        this->L_diag = L_diag; // assign the calculated diagonals to be parts of our Solver.
        this->U_upDiag = U_upDiag;
        
        return;
    }
    
    // Default destructor
    virtual ~TriSolver(){};
    
    
    // given the vector f, we calculate a solution u.
    
    void apply(long systemSize, vector<double>& f, vector<double>& u)
    {
        // Have the LU factorization A = LU computed in constructor
        
        // Step 1 -- Solve the system Lv = f
        vector<double> v; // vector we use as intermediate step of solving this system Lv = f
        v.resize(systemSize);
        
        v[0] = f[0] / L_diag[0];
        for(long i=1; i <= systemSize-1; i ++){
            v[i] = (f[i] - L_loDiag[i-1]*v[i-1])/L_diag[i];
        }
        
        // Step 2 -- Solve the system Uu = v
        u[systemSize-1] = v[systemSize-1];
        for(long i=systemSize-2; i >= 0; i--){
            u[i] = v[i] - U_upDiag[i]*u[i+1];
        }
        
        // delete v, clean up space
        v.clear();
        
        return;
    }
    
    // internal data
    long systemSize ;
    vector<double> L_loDiag;
    vector<double> L_diag;
    vector<double> U_upDiag;
};

#endif /* TriSolver_h */

