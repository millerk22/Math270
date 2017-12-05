#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <vector>
using namespace std;

//
// Math 270A : Assignment 9 Option II
//
// Test program for fixed step relaxation to a solution of the variable coefficient problem
//
// d/dx( a(x,y) du/dx )  +  d/dy( a(x,y) du/dy ) = f(x,y)  with Dirichlet boundary conditions
//
// for (x,y) in [xMin, xMax]x[yMin,yMax]
//
//
// Date: Mon 20 Nov 2017 08:36:39 AM PST
//

#include "GridFun2D.h"       // 2D grid function class
#include "VarRelaxOp2D.h"    // 2D Crank-Nicholson relaxation operator
#include "GNUplotUtility.h"

//
// This test problem for
//
// d/dx( a(x,y) du/dx )  +  d/dy( a(x,y) du/dy ) = f(x,y)  with Dirichlet boundary conditions
//
// is created by using the variable coefficient
//
// a(x,y) = (1 + ((x-xMin)/(xMax-xMin) + ((y-yMin)/(yMax-yMin))
//
// and obtaining the right hand side that corresponds to the solution
//
// u(x,y) = cos(2*pi*waveNumberX*(x-xMin)/(xMax-xMin)) * cos(2*pi*waveNumberY*(y-yMin)/(yMax-yMin))
//
// by inserting this solution into the equation and analytically differentiating.
//
//
class VarCoefficientTestProblem2D
{
public:
    
    VarCoefficientTestProblem2D(double xMin, double xMax, double waveNumberX, double yMin, double yMax, double waveNumberY)
    {
        this->xMin        = xMin;
        this->xMax        = xMax;
        this->waveNumberX = waveNumberX;
        this->yMin        = yMin;
        this->yMax        = yMax;
        this->waveNumberY = waveNumberY;
        this->pi          = 3.1415926535897932;
    }
    
    /// Returns the solution
    
    double operator()(double x,double y)
    {
        return  cos(2.0*pi*waveNumberX*((x-xMin)/(xMax-xMin)))*cos(2.0*pi*waveNumberY*((y-yMin)/(yMax-yMin)));
    }
    //
    // Returns f(x,y) = d/dx( a(x,y) du/dx )  +  d/dy( a(x,y) du/dy )
    //
    // where  a(x,y) =  (1.0 + (x - xMin)/(xMax-xMin) + (y - yMin)/(yMax-yMin))
    //
    // and  u(x,y) = cos(2*pi*waveNumberX*(x-xMin)/(xMax-xMin))*cos(2*pi*waveNumberY*(y-yMin)/(yMax-yMin))
    //
    //
    double evaluateRHS(double x, double y)
    {
        double Lx = (xMax-xMin);
        double Ly = (yMax-yMin);
        return    -(1.0 + (x - xMin)/Lx  + (y - yMin)/Ly)*
        (4.0*pi*pi*waveNumberX*waveNumberX*(1.0/(Lx*Lx)) + 4.0*pi*pi*waveNumberY*waveNumberX*(1.0/(Ly*Ly)) )
        *(cos(2.0*pi*waveNumberX*((x-xMin)/Lx))*cos(2.0*pi*waveNumberY*((y-yMin)/Ly)))
        
        - ((1.0/Lx)*2.0*pi*(waveNumberX/Lx)*sin(2.0*pi*waveNumberX*((x-xMin)/Lx))*cos(2.0*pi*waveNumberY*((y-yMin)/Ly)))
        - ((1.0/Ly)*2.0*pi*(waveNumberY/Ly)*sin(2.0*pi*waveNumberY*((y-yMin)/Ly))*cos(2.0*pi*waveNumberX*((x-xMin)/Lx)));
    }
    
    /// Returns the values of the variable coefficient
    
    double evaluateCoeff(double x,double y)
    {
        return (1.0 + (x - xMin)/(xMax-xMin) + (y - yMin)/(yMax-yMin));
    }
    
    
    double normInf()
    {
        return 1.0;
    }
    
    double pi;
    double xMin; double xMax; double waveNumberX;
    double yMin; double yMax; double waveNumberY;
};

/// Returns the size of the residual at interior points

double evaluateResidualNorm(const GridFun2D& a, const GridFun2D& u, const GridFun2D& f)
{
    double residualNorm = 0.0;
    double residualValue;
    
    double hx = u.hx;
    double hy = u.hy;
    
    for(long i = 1; i <  u.xPanel; i++)
    {
        for(long j = 1; j <  u.yPanel; j++)
        {
            residualValue  = (  0.5*(a.values(i+1,j)  + a.values(i,j))*u.values(i+1,j)
                              -   (0.5*a.values(i+1,j)   + a.values(i,j) + 0.5*a.values(i-1,j))*u.values(i,j)
                              +    0.5*(a.values(i,j)    + a.values(i-1,j))*u.values(i-1,j)) /(hx*hx)
            
            +(  0.5*(a.values(i,j+1)  + a.values(i,j))*u.values(i,j+1)
              -   (0.5*a.values(i,j+1)   + a.values(i,j) + 0.5*a.values(i,j-1))*u.values(i,j)
              +    0.5*(a.values(i,j)    + a.values(i,j-1))*u.values(i,j-1)) /(hy*hy)
            
            - f.values(i,j);
            
            if(abs(residualValue) > residualNorm){residualNorm = abs(residualValue);}
        }}
    return residualNorm;
}

int main()
{
    // Toggle between fixed and variable timestepping
    bool fixedTimeStep = true;
    // Set up test problem
    
    double waveNumberX =  1.0;  // test problem x-coordinate wave number
    double waveNumberY =  1.0;  // test problem x-coordinate wave number
    
    long N;
    cout << "NumPanels:" << endl;
    cin >> N;
    long   xPanel      =   N;  // X panel count
    long   yPanel      =   N;  // X panel count
    
    double dt          =   0.025;  // Relaxation timestep
    
    double tol        =  1.0e-6;  // Stopping tolerance
    
    double xMin        =  0.0;
    double xMax        =  1.0;
    double hx          = (xMax-xMin)/(double)xPanel;
    
    double yMin        =  0.0;
    double yMax        =  1.0;
    double hy          = (yMax-yMin)/(double)yPanel;
    
    
    // Echo input parameters
    
    cout << "XXXX  Douglas ADI 2D (Variable Coefficient) Program Start      XXXX " << endl;
    cout << "X Panel Count : " << xPanel << endl;
    cout << "Y Panel Count : " << xPanel << endl;
    cout << "Timestep      : " << dt << endl;
    
    
    // Instantiate the test problem solution
    
    VarCoefficientTestProblem2D testSoln(xMin, xMax, waveNumberX, yMin, yMax, waveNumberY);
    
    // Extract coefficients, right hand side values, and exact
    // solution from the test problem instance
    
    
    GridFun2D aCoeff(xPanel,xMin,xMax,yPanel,yMin,yMax);
    GridFun2D f(xPanel,xMin,xMax,yPanel,yMin,yMax);
    GridFun2D uExact(xPanel,xMin,xMax,yPanel,yMin,yMax);
    
    
    double x; double y;
    
    for(long i = 0; i <= xPanel; i++)
    {
        x = xMin + i*hx;
        for(long j = 0; j <= yPanel; j++)
        {
            y  = yMin + j*hy;
            aCoeff.values(i,j) = testSoln.evaluateCoeff(x,y);
            uExact.values(i,j) = testSoln(x,y);
            f.values(i,j)      = testSoln.evaluateRHS(x,y);
        }}
    
    GridFun2D uk(xPanel,xMin,xMax,yPanel,yMin,yMax);
    GridFun2D ukp1(xPanel,xMin,xMax,yPanel,yMin,yMax);
    GridFun2D uErr(xPanel,xMin,xMax,yPanel,yMin,yMax);
    
    //
    // Set initial iterate with interior values equal to zero
    // and perimeter (boundary) values set the exact solution.
    //
    // Set boundary values of f to zero.
    
    uk.setToValue(0.0);
    
    long i; long j;
    for(i = 0; i <= xPanel; i++)
    {
        j   = 0;
        x   =   xMin + i*hx;
        y   =   yMin + j*hy;
        uk.values(i,j) = testSoln(x,y);
        f.values(i,j)  = 0.0;
        
        j   = yPanel;
        x   = xMin + i*hx;
        y   = yMin + j*hy;
        uk.values(i,j) = testSoln(x,y);
        f.values(i,j)  = 0.0;
    }
    
    for(long j = 0; j <= yPanel; j++)
    {
        i = 0;
        x = xMin + i*hx;
        y = yMin + j*hy;
        
        uk.values(i,j) = testSoln(x,y);
        f.values(i,j)  = 0.0;
        
        i = xPanel;
        x = xMin + i*hx;
        y = yMin + j*hy;
        
        uk.values(i,j) = testSoln(x,y);
        f.values(i,j)  = 0.0;
    }
    
    
    // Output right hand side, initial data and exact solution
    
    GNUplotUtility::output(f,"f.dat");
    GNUplotUtility::output(aCoeff,"a.dat");
    GNUplotUtility::output(uk,"u0.dat");
    GNUplotUtility::output(uExact,"uExact.dat");
    
    
    VarRelaxOp2D varRelaxOp;
    varRelaxOp.initialize(dt, aCoeff, f);
    
    //
    // Relaxation loop
    //
    
    long   iterMax   = 1000;
    long   iter      = 0;
    double residNorm = 2.0*tol;
    
    //###################################################
    //         Fixed timestep
    //###################################################
    
    if(fixedTimeStep)
    {
        
        while((residNorm > tol)&&(iter < iterMax))
        {
            varRelaxOp.apply(uk,ukp1);
            uk  = ukp1;
            
            
            // Evaluate the residual at interior points
            
            residNorm = evaluateResidualNorm(aCoeff, uk, f);
            
            iter++;
            if(iter%10 == 0) cout << "Iteration : " << iter << endl;
        }
    }
    //###################################################
    //         Variable timestep
    //
    //  ToDo: Determine optimal parameters
    //
    //###################################################
    else
    {
        
        long sweepMax = 20;
        
        
        double h0  =  hx;
        h0  = (h0 < hy ) ? h0 : hy;
        
        long  panelCount = xPanel;
        panelCount = (xPanel > yPanel) ? xPanel : yPanel;
        
        double alphaBar   =  0.0;
        double* aData     = aCoeff.values.getDataPointer();
        for(long i = 0; i < aCoeff.values.getDataSize(); i++)
        {
            alphaBar += aData[i];
        }
        alphaBar = alphaBar/aCoeff.values.getDataSize();
        
        long levelMax      = (long)(log((double)panelCount)/log(2.0)) + 1;
        
        // Create a vector of operators for each timestep size
        
        vector < VarRelaxOp2D > varRelaxOp2Darray(levelMax+1);
        
        double  hj;
        double dtj;
        for(long j = 0; j <= levelMax; j++)
        {
            hj = h0*pow(2,j);
            dtj = (4.0*hj*hj)/(alphaBar*3.1415926535897932*3.1415926535897932);
            
            varRelaxOp2Darray[j].initialize(dtj, aCoeff, f);
        }
        
        long sweepCount = 0;
        residNorm = 100*tol;
        
        while((residNorm > tol)&&(sweepCount < sweepMax))
        {
            for(long j = 0; j <= levelMax; j++)
            {
                varRelaxOp2Darray[j].apply(uk,ukp1);
                uk=ukp1;
                iter++;
            }
            for(long j = levelMax; j >= 0; j--)
            {
                varRelaxOp2Darray[j].apply(uk,ukp1);
                uk=ukp1;
                iter++;
            }
            
            // Evaluate the residual at interior points
            
            residNorm = evaluateResidualNorm(aCoeff, uk, f);
            
            printf("Sweep Residual : %10.5e \n",residNorm);
            
            sweepCount++;
        }
        
        printf("Sweep      Count : %ld  \n",sweepCount);
        printf("Relaxation Count : %ld  \n",iter);
    }
    
    //
    // Evaluate the error with respect to the continuous solution
    //
    
    uErr =         uk;
    uErr -=    uExact;
    
    GNUplotUtility::output(uk,"uK.dat");
    GNUplotUtility::output(uErr,"uErr.dat");
    double solnNorm     = uExact.normInf();
    double rhsNorm      = f.normInf();
    double maxError     = uErr.normInf();
    double relError     = maxError/solnNorm;
    double residualNorm = evaluateResidualNorm(aCoeff, uk, f);
    //
    // Print out the results (using standard C i/o because I find it
    // easier to create nice output.
    //
    
    printf("XXXX   Fixed Timestep Relaxation Output XXXX \n\n");
    printf("X Panel Count : %ld   \n", xPanel);
    printf("Y Panel Count : %ld   \n", yPanel);
    printf("Timestep    : %10.5e \n",dt);\
    printf("Number of relaxation steps  : %ld \n",iter);
    printf("||   Residual     || : %10.5e \n",residNorm);
    printf("||   Residual     ||/ || F ||        : %10.5e \n",residualNorm/rhsNorm);
    printf("||  Relative Error|| : %10.5e \n", relError);
    printf("||Discrete - Exact||/||Exact||       : %10.5e \n",maxError);
    
    return 0;
}
