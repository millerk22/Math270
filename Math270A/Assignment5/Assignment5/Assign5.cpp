#include <cstdio>
#include <cmath>
using namespace std;

//
// Math 270A : Assignment 5
//
// Test program for variable step relaxation to a solution of
//
// alphaX * d^2 U/dx^2  alphaY * d^2 U/dy^2 = f(x,y)  with Dirichlet boundary conditions
//
// for (x,y) in [xMin, xMax]X[yMin,yMax] using the 2D version of Douglas ADI
//
//
// Date: Nov 3., 2017
//

#include "GridFun2D.h"    // 2D grid function class
#include "RelaxOp2D.h"    // 2D Douglas relaxation operator class

#include "GNUplotUtility.h"

class TestProblem2D
{
public:
    
    TestProblem2D(double alphaX, double waveNumberX,double xMin, double xMax,
                  double alphaY, double waveNumberY,double yMin, double yMax)
    {
        this->alphaX      = alphaX;
        this->alphaY      = alphaY;
        this->waveNumberX = waveNumberX;
        this->waveNumberY = waveNumberY;
        this->pi          = 3.1415926535897932;
    }
    
    double operator()(double x, double y)
    {
        double d2Xfact = (2.0*pi*waveNumberX)*(2.0*pi*waveNumberX);
        double d2Yfact = (2.0*pi*waveNumberY)*(2.0*pi*waveNumberY);
        return -(cos(2.0*pi*x*waveNumberX)*cos(2.0*pi*y*waveNumberY))/(alphaX*d2Xfact + alphaY*d2Yfact);
    }
    
    double normInf()
    {
        double d2Xfact = (2.0*pi*waveNumberX)*(2.0*pi*waveNumberX);
        double d2Yfact = (2.0*pi*waveNumberY)*(2.0*pi*waveNumberY);
        return 1.0/(alphaX*d2Xfact + alphaY*d2Yfact);
    }
    
    double pi;
    double alphaX;
    double alphaY;
    double waveNumberX;
    double waveNumberY;
};

//
// This routine computes the residual of the five-point discretization of Poisson's
// equation with Dirichlet boundary conditions at interior points of a rectangular domain.
//
// Input : Approximate solution u

//         Right hand side f which specifies the forcing function at interior grid points
//         (boundary values of f are ignored)
//
// Output : The residual obtained by evaluating
//
// r = f - (alpha_x*Delta_x + alpha_y*Delta_y) u   at interior points
//
// r = 0                                           at boundary points
//

void evaluatePoissonResidual(double alphaX, double alphaY, GridFun2D& u, GridFun2D& f, GridFun2D& residual)
{
    residual.initialize(f);   // Initialize residual with f
    residual.setToValue(0.0); // Set values to 0.0
    
    // Residual at interior points
    
    double hy = u.hy;
    double hx = u.hx;
    for(long i = 1; i <  u.xPanel; i++)
    {
        for(long j = 1; j <  u.yPanel; j++)
        {
            residual.values(i,j)  = f.values(i,j) - (alphaX*((u.values(i+1,j) - 2.0*u.values(i,j) + u.values(i-1,j))/(hx*hx))
                                                     + alphaY*((u.values(i,j-1) - 2.0*u.values(i,j) + u.values(i,j+1))/(hy*hy)));
        }
    }
}

int main()
{
    
    // Set up test problem
    
    double alphaX      =  2.0;  // Laplace operator-x prefactor
    double alphaY      =  2.0;  // Laplace operator-y prefactor
    double waveNumberX =  1.0;  // test problem x-coordinate wave number
    double waveNumberY =  1.0;  // test problem y-coordinate wave number
    
    // get the user input to specify the panel count
    long INPUT;
    cout << "Please input the panel count " << endl;
    cin >> INPUT;
    long   xPanel     =   INPUT;  // X panel count
    long   yPanel     =   INPUT;  // Y panel count
    double tol    =  1.0e-6;  // Stopping tolerance
    
    double xMin   =     0.0;
    double xMax   =     1.0;
    double hx     = (xMax-xMin)/(double)xPanel;
    
    double yMin   =     0.0;
    double yMax   =     1.0;
    double hy     = (yMax-yMin)/(double)yPanel;
    
    // List of files output
    
    vector<string> outputFileList;
    
    // Echo input parameters
    
    cout << "XXXX   Douglas ADI Program Start      XXXX " << endl;
    cout << "XXXX   Variable Time Stepping         XXXX " << endl;
    cout << "alpha_x : " << alphaX << endl;
    cout << "alpha_y : " << alphaY << endl;
    cout << "X Panel Count : " << xPanel << endl;
    cout << "Y Panel Count : " << yPanel << endl;
    cout << "Stopping Tolerance : " << tol << endl;
    //cout << "Timestep      : " << dt << endl;
    
    //
    // Instantiate the test problem solution
    //
    
    TestProblem2D     testSoln(alphaX,waveNumberX,xMin,xMax,
                               alphaY,waveNumberY,yMin,yMax);
    
    double solnNorm = testSoln.normInf();
    
    //
    // Instantiate 2D grid functions
    //
    
    GridFun2D f(xPanel,xMin,xMax,yPanel,yMin,yMax);
    GridFun2D uk(xPanel,xMin,xMax,yPanel,yMin,yMax);
    GridFun2D ukp1(xPanel,xMin,xMax,yPanel,yMin,yMax);
    GridFun2D uTmp(xPanel,xMin,xMax,yPanel,yMin,yMax);
    
    GridFun2D uExact(xPanel,xMin,xMax,yPanel,yMin,yMax);
    GridFun2D uSolnErr(xPanel,xMin,xMax,yPanel,yMin,yMax);
    GridFun2D residual(xPanel,xMin,xMax,yPanel,yMin,yMax);
    
    //
    // Construct the right hand side and the exact solution
    //
    
    double x; double y;
    double pi =  3.1415926535897932;
    for(long i = 0; i <= xPanel; i++)
    {
        x = xMin + i*hx;
        for(long j = 0; j <= yPanel; j++)
        {
            y =  yMin + j*hy;
            f.values(i,j)      = cos(2.0*pi*x*waveNumberX)*cos(2.0*pi*y*waveNumberY);
            uExact.values(i,j) = testSoln(x,y);
        }
    }
    
    
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
    outputFileList.push_back("f.dat      : Forcing function");
    
    GNUplotUtility::output(uk,"u0.dat");
    outputFileList.push_back("u0.dat     : Starting iterate");
    
    GNUplotUtility::output(uExact,"uExact.dat");
    outputFileList.push_back("uExact.dat : Exact solution");
    
    
    //
    /// Find all the dtj's so we don't have to calculate each time
    //
    
    long M = xPanel; // find M = max(xPanel, yPanel)
    if (yPanel > xPanel){
        M = yPanel;
    }
    int levels = int(log2(M)+ 1); // levels is the number of different dtj's we will do in a sweep
    vector<float> dTJ;
    vector<RelaxOp2D> RO2D_levels;
    dTJ.resize(levels);
    RO2D_levels.resize(levels);
    for (int i=0; i < levels; i++){
        dTJ[i] = 4.0*hx*hx*pow(2.0,i)/(pi*pi*alphaX); // note hx = hy = h0 and alphaX = alphaY = alpha
        RO2D_levels[i].initialize(dTJ[i], alphaX, alphaY, f);
    }
    
    
    //
    // Relaxation loop, with Variable Timestepping this time
    //
    
    long   iterMax  = 1000;
    long   iter     = 0;
    long iters_levels = 0;
    double residualNormInf = 2*tol;
    
    while((residualNormInf > tol)&&(iter < iterMax))
    {
        // forward sweep
        uTmp = uk;
        
        for(int j=0; j < levels; j++){
            RO2D_levels[j].apply(uTmp, ukp1);
            uTmp = ukp1;
        }
        
        // backward sweep
        for(int j=levels-1; j>= 0; j--){
            RO2D_levels[j].apply(uTmp, ukp1);
            uTmp = ukp1;
        }
        
        iters_levels += 2*levels;
        
        // Evaluate the residual and the relative norm of the residual
        evaluatePoissonResidual(alphaX,  alphaY, ukp1, f, residual);
        
        residualNormInf = residual.normInf();
        residualNormInf /= f.normInf();  // do we need this?
        
        
        // Update iterate
        
        uk  = ukp1;
        iter++;
        
        if(iter%10 == 0) cout << "Iteration : " << iter << endl;
    }
    
    
    
    // output the final residual's values to a dat file
    GNUplotUtility::output(residual,"residual.dat");
    outputFileList.push_back("Residual   : residual.dat");
    
    
    // Evaluate the error with respect to the continuous solution
    
    uSolnErr =         uk;
    uSolnErr -=    uExact;
    
    GNUplotUtility::output(uk,"uK.dat");
    outputFileList.push_back("uK.dat     : Computed solution");
    
    GNUplotUtility::output(uSolnErr,"uSolnErr.dat");
    outputFileList.push_back("uSolnErr.dat : Computed solution error");
    
    double absSolnError    = uSolnErr.normInf();
    double relSolnError    = absSolnError/solnNorm;
    
    //
    // Print out the results using standard C I/0
    //
    
    printf("\nXXXX   Variable Timestep Relaxation Output XXXX \n\n");
    printf("Number of levels per sweep         : %1d \n", levels);
    printf("Number of sweeps               : %ld \n",iter);
    printf("Total Number of iterations      : %1d \n", iters_levels);
    printf("Relative residual Size (InfNorm)   : %10.5e \n",residualNormInf);
    printf("Absolute solution error (InfNorm)  : %10.5e \n",absSolnError);
    printf("Relative solution error (InfNorm)  : %10.5e \n",relSolnError);
    
    // Print out list of files output
    
    printf("\nXXXX   Files Output  XXXX \n");
    
    for(size_t k = 0; k < outputFileList.size(); k++)
    {
        cout << outputFileList[k] << endl;
    }
    
    
    return 0;
}

