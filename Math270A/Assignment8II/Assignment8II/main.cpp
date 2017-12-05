#include <cstdio>
#include <cmath>
#include <cstdlib>
using namespace std;

//
// Math 270A : Assignment 8 Option II
//
// Test program Douglas ADI with a fixed step to obtain a solution to the variable coefficient problem
//
// d/dx( a(x) du/dx )   = f(x)  with Dirichlet boundary conditions
//
// for x in [xMin, xMax].
//
//
// Date: Mon 20 Nov 2017 08:32:41 AM PST
//

#include "GridFun1D.h"       // 1D grid function class
#include "VarRelaxOp1D.h"    // 1D Crank-Nicholson relaxation operator


//
// This test problem for
//
// d/dx( a(x) du/dx )   = f(x)  with Dirichlet boundary conditions
//
// is created by using the variable coefficient
//
// a(x) = (1 + ((x-xMin)/(xMax-xMin))
//
// and obtaining the right hand side that corresponds to the solution
//
// u(x) = cos(2*pi*waveNumberX*(x-xMin)/(xMax-xMin))
//
// by inserting this solution into the equation and analytically differentiating.
//
//
class VarCoefficientTestProblem1D
{
public:
    
    VarCoefficientTestProblem1D(double xMin, double xMax, double waveNumberX)
    {
        this->xMin        = xMin;
        this->xMax        = xMax;
        this->waveNumberX = waveNumberX;
        this->pi          = 3.1415926535897932;
    }
    
    /// Returns the solution
    
    double operator()(double x)
    {
        return  cos(2.0*pi*waveNumberX*((x-xMin)/(xMax-xMin)));
    }
    //
    // Returns d/dx ( (1 + ((x-xMin)/(xMax-xMin)) du/dx )
    //
    // where   L = (xMax-xMin)
    // and  u(x) = cos(2*pi*waveNumberX*(x-xMin)/L)
    //
    //
    double evaluateRHS(double x)
    {
        double L = (xMax-xMin);
        return    (1.0 + (x - xMin)/L)*(-4.0*pi*pi*waveNumberX*waveNumberX*(1.0/(L*L))*cos(2.0*pi*waveNumberX*((x-xMin)/L)))
        - ((1.0/L)*2.0*pi*(waveNumberX/L)*sin(2.0*pi*waveNumberX*((x-xMin)/L)));
    }
    
    /// Returns the values of the variable coefficient
    
    double evaluateCoeff(double x)
    {
        return (1.0 + (x - xMin)/(xMax-xMin));
    }
    
    
    double normInf()
    {
        return 1.0;
    }
    
    
    double pi;
    double xMin;
    double xMax;
    double waveNumberX;
};

/// Returns the size of the residual at interior points

double evaluateResidualNorm(const GridFun1D& a, const GridFun1D& u, const GridFun1D& f)
{
    double residualNorm = 0.0;
    double residualValue;
    
    double hx = u.hx;
    
    for(long i = 1; i <  u.xPanel; i++)
    {
        residualValue  = (  0.5*(a.values[i+1]  + a.values[i])*u.values[i+1]
                          -   (0.5*a.values[i+1]   + a.values[i] + 0.5*a.values[i-1])*u.values[i]
                          +    0.5*(a.values[i]    + a.values[i-1])*u.values[i-1]) /(hx*hx) - f.values[i];
        
        if(abs(residualValue) > residualNorm){residualNorm = abs(residualValue);}
    }
    return residualNorm;
}

int main()
{
    // Set up test problem
    
    double waveNumberX =  1.0;  // test problem x-coordinate wave number
    long   xPanel      =   40;  // X panel count
    double dt          =   .1;  // Relaxation timestep
    
    double tol        =  1.0e-6;  // Stopping tolerance
    
    double xMin        =  0.0;
    double xMax        =  1.0;
    double hx          = (xMax-xMin)/(double)xPanel;
    
    
    // Echo input parameters
    
    cout << "XXXX   Crank-Nicholson Program Start      XXXX " << endl;
    cout << "X Panel Count : " << xPanel << endl;
    cout << "Timestep      : " << dt << endl;
    
    
    // Instantiate the test problem solution
    
    VarCoefficientTestProblem1D testSoln(xMin, xMax, waveNumberX);
    
    // Extract coefficients, right hand side values, and exact
    // solution from the test problem instance
    
    
    GridFun1D aCoeff(xPanel,xMin,xMax);
    GridFun1D      f(xPanel,xMin,xMax);
    GridFun1D uExact(xPanel,xMin,xMax);
    GridFun1D residual(xPanel,xMin,xMax);
    
    double x;
    for(long i = 0; i <= xPanel; i++)
    {
        x                = xMin + i*hx;
        aCoeff.values[i] = testSoln.evaluateCoeff(x);
        uExact.values[i] = testSoln(x);
        f.values[i]      = testSoln.evaluateRHS(x);
    }
    
    // Set boundary values of f to zero.
    
    f.values[0]      = 0.0;
    f.values[xPanel] = 0.0;
    
    GridFun1D   uk(xPanel,xMin,xMax);
    GridFun1D ukp1(xPanel,xMin,xMax);
    
    // Set up initial iterate: zero in the interior, boundary
    // values at the edges.
    
    uk.setToValue(0.0);
    uk.values[0] = testSoln(xMin);
    uk.values[xPanel] = testSoln(xMax);
    
    // Initialize the variable coefficient relaxation operator
    
    VarRelaxOp1D relaxOp;
    relaxOp.initialize(dt,aCoeff,f);
    
    //
    // Initialize relaxation loop
    //
    
    double diffNorm = 2*tol;
    long   iterMax  = 4000;
    long   iter     = 0;
    
    while((diffNorm > tol)&&(iter < iterMax))
    {
        relaxOp.apply(uk,ukp1);
        
        // Check difference between iterates
        
        uk      -= ukp1;
        diffNorm = uk.normInf()/dt;
        
        // Update iterate
        
        uk  = ukp1;
        iter++;
        if(iter%10 == 0) cout << "Step : " << iter << " Relative Difference Norm : " << diffNorm << endl;
    }
    
    // Evaluate the error
    
    double errorMax = 0.0;
    double errVal   = 0.0;
    for(long i = 0; i <= xPanel; i++)
    {
        x = xMin + i*hx;
        errVal = abs(uk.values[i] -testSoln(x));
        if(errorMax < errVal) {errorMax = errVal;}
    }
    
    double residualNorm = evaluateResidualNorm(aCoeff, uk, f);
    //
    // Print out the results (using standard C i/o because I find it
    // easier to create nice output.
    //
    
    printf("XXXX   Fixed Timstep Relaxation Output XXXX \n\n");
    printf("Panel Count : %ld   \n", xPanel);
    printf("Timestep    : %10.5e \n",dt);
    printf("Difference between iterates : %10.5e \n",diffNorm);
    printf("Number of relaxation steps  : %ld \n",iter);
    printf("||Discrete - Exact||        : %10.5e \n",errorMax);
    printf("||   Residual     ||        : %10.5e \n",residualNorm);
    
    return 0;
}
