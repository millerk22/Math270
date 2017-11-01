#include <cstdio>
#include <cmath>
#include <cstdlib>
using namespace std;

//
// Math 270A : Assignment 2
//
// Test program for fixed step relaxation to a solution of
//
// alpha * d^2 U/dx^2  = f  with Dirichlet boundary conditions
//
// for x in [xMin, xMax]
//
//
// Version : Tues. Oct. 10, 2017 06:48:20 PM PST
//
#include "GridFun1D.h"  // 1D grid function class
#include "RelaxOp1D.h"  // 1D Crank-Nicholson relaxation operator class

class TestProblem1D
{
public:
    
    TestProblem1D(double alpha)
    {
        this->alpha = alpha;
        this->pi    = 3.1415926535897932;
    }
    
    double operator()(double x)
    {
        return -cos(2.0*pi*x)/(4.0*pi*pi*alpha);
    }
    
    double pi;
    double alpha;
    
};


// include a class for applying the matrix A to a GridFun1D?

class ApplyOriginal
{
public:
    ApplyOriginal(double alpha, double h, long systemSize){
        this-> alpha = alpha;
        this-> h = h;
        this-> systemSize = systemSize;
    }
    
    void apply(GridFun1D& u_calc, GridFun1D& Au){
        double val = alpha/(h*h);
        
        // apply alpha * D_D+ to the interior points
        for (long i=1; i < systemSize-1; i++){
            Au.values[i] = val*(u_calc.values[i-1] - 2.0*u_calc.values[i] + u_calc.values[i+1]);
        }
        // note that we apply identity to endpoints
        Au.values[0] = u_calc.values[0];
        Au.values[systemSize-1] = u_calc.values[systemSize-1];
        return;
    }
    
    double alpha;
    double h;
    long systemSize;
};


int main()
{
    vector<long> test_Ms = {10, 20, 40};
    for(int k = 0; k < 3; k++){
        //
        // Set up test problem
        //
        double alpha  =     2.0;  // Laplace operator prefactor
        double dt     =     0.1;  // Relaxation timestep
        long   M      =     test_Ms[k];   // Panel count
        double tol    =  1.0e-6;  // Stopping tolerance
        
        
        double xMin   =     0.0;
        double xMax   =     1.0;
        double h      = (xMax-xMin)/(double)M;
        
        
        
        
        //
        // Instantiate the test problem solution
        //
        
        TestProblem1D      testSoln(alpha);
        
        
        // Instantiate grid functions for a discretization
        // of [xMin,xMax] with M panels.
        //
        
        GridFun1D f(M,xMin,xMax);
        GridFun1D uk(M,xMin,xMax);
        GridFun1D ukp1(M,xMin,xMax);
        GridFun1D uTmp(M,xMin,xMax);
        
        // Instantiate the Apply Original system operator for residual calculation after sweeps
        ApplyOriginal A(alpha, h, M+1);
        
        
        
        //
        // Initialize right hand side
        //
        
        double x;
        double pi =  3.1415926535897932;
        for(long i = 0; i <= M; i++)
        {
            x          =   xMin + i*h;
            f.values[i] = cos(2.0*pi*x);
        }
        
        // Set edge values of f to zero
        
        f.values[0] = 0.0;
        f.values[M] = 0.0;
        
        // make f_orig for the residual calculation for iterates
        GridFun1D f_orig = f;
        f_orig.values[0] = testSoln(xMin);
        f_orig.values[M] = testSoln(xMax);
        
        
        
        
        
        // Calculate all the dtj's once so we can just iterate thru a vector to get the timesteps for a sweep. Also make a vector of RelaxOp1D's, 1 for each different dtj.
        int levels = int(log2(double(M))) + 1;
        vector<float> dTJ ;
        vector<RelaxOp1D> RO1D_levels;
        dTJ.resize(levels);
        RO1D_levels.resize(levels);
        for (int i=0; i < levels; i++){
            dTJ[i] = 4.0*h*h*pow(2.0, i)*pow(2.0,i)/(pi*pi*alpha);
            RO1D_levels[i].initialize(dTJ[i], alpha, f);
        }
        
        
        
        /* // Old relax operator, for a fixed timestep
        //
        // Instantiate and initialize relaxation operator
        //
        
        RelaxOp1D relaxOp;
        relaxOp.initialize(dt,alpha,f);
        */
        
        
        
        //
        // Set up initial iterate: zero in the interior, boundary
        // values at the edges.
        //
        
        uk.setToValue(0.0);
        uk.values[0] = testSoln(xMin);
        uk.values[M] = testSoln(xMax);
        
        //
        // Initialize relaxation loop
        //
        
        
        double diffNorm = 2*tol;
        long   iterMax  = 4000;
        long   iter     = 0;
        long   iters_levels = 0;
        
        while((diffNorm > tol)&&(iter < iterMax))
        {
            // forward sweep
            uTmp = uk;
            
            for (int j=0; j < levels; j++){
                RO1D_levels[j].apply(uTmp, ukp1);
                uTmp = ukp1;
            }
            
            // backward sweep
            for (int j=levels-1; j >= 0; j --){
                RO1D_levels[j].apply(uTmp, ukp1);
                uTmp = ukp1;
            }
            
            iters_levels += 2*levels;
            
            //relaxOp.apply(uk,ukp1); // DONT NEED, from old
            
            /* // Old check for convergence, can't use because what's dt? and the other things said in class
            // Check relative difference between iterates
            
            uk       -= ukp1;
            diffNorm  = uk.normInf();   // Norm of time derivative
            diffNorm /= dt;
            */
            
            // Residual of original Au = f after a complete sweep.
            A.apply(ukp1, uk); // take the new iterate ukp1 and apply A to it, use prev uk to hold the values for A*ukp1
            uk -= f_orig;
            diffNorm = uk.normInf(); // norm of the residual
            
            
            // Update iterate
            
            uk  = ukp1;
            iter++;
        }
        
        
        //
        // Evaluate the error
        //
        
        double errorMax = 0.0;
        double errVal   = 0.0;
        for(long i = 0; i <= M; i++)
        {
            x = xMin + i*h;
            errVal = abs(uk.values[i] -testSoln(x));
            if(errorMax < errVal) {errorMax = errVal;}
        }
        
        //
        // Print out the results (using standard C I/O because I find it
        // easier to create nice output).
        //
        
        printf("XXXX   Sweep Timestep Relaxation Output XXXX \n\n");
        printf("Panel Count : %ld  \n",M);
        printf("Num levels    : %i \n",levels);
        printf("Residual Norm : %10.5e \n",diffNorm);
        printf("Number of sweeps       : %ld \n",iter);
        printf("Tot Number of iterations       : %ld \n",iters_levels);
        printf("Iterative solution error    : %10.5e \n",errorMax);
        
        cout << endl << endl;
        
        
        
        
        
        // Last week's assignment, output to compare
        
        
        
        
        //
        // Instantiate and initialize relaxation operator
        //
        
        RelaxOp1D relaxOp;
        relaxOp.initialize(dt,alpha,f);
        
        //
        // Set up initial iterate: zero in the interior, boundary
        // values at the edges.
        //
        
        uk.setToValue(0.0);
        ukp1.setToValue(0.0);
        uTmp.setToValue(0.0);
        uk.values[0] = testSoln(xMin);
        uk.values[M] = testSoln(xMax);
        
        //
        // Initialize relaxation loop
        //
        
        
        diffNorm = 2*tol;
        iterMax  = 4000;
        iter     = 0;
        
        while((diffNorm > tol)&&(iter < iterMax))
        {
            
            relaxOp.apply(uk,ukp1);
            
            
            // Check relative difference between iterates
            
            uk       -= ukp1;
            diffNorm  = uk.normInf();   // Norm of time derivative
            diffNorm /= dt;
            
            // Update iterate
            
            uk  = ukp1;
            iter++;
        }
        
        //
        // Evaluate the error
        //
        
        errorMax = 0.0;
        errVal   = 0.0;
        for(long i = 0; i <= M; i++)
        {
            x = xMin + i*h;
            errVal = abs(uk.values[i] -testSoln(x));
            if(errorMax < errVal) {errorMax = errVal;}
        }
        
        //
        // Print out the results (using standard C I/O because I find it
        // easier to create nice output).
        //
        
        printf("XXXX   Fixed Timstep Relaxation Output XXXX \n\n");
        printf("Panel Count : %ld  \n",M);
        printf("Timestep    : %10.5e \n",dt);
        printf("Difference between iterates : %10.5e \n",diffNorm);
        printf("Number of iterations        : %ld \n",iter);
        printf("Iterative solution error    : %10.5e \n",errorMax);
        
        cout << endl << endl << endl;
        
    
    
    }
    
    
    
    
    
    return 0;
}

