#include <cstdio>
#include <iostream>
#include <cmath>
using namespace std;
//
// LaplaceMP.cpp
//
// A test code demonstrating the computation of the two-dimensional discrete Laplace
// operator using different OMP constructs. There are 4 methods of computation
//
// 1) Standard nested for loop computation
//
// 2) Nested for loop computation with multi-threaded implementation of
//    inner loop.
//
// 3) X-Y sweep computation using the application of a one-dimensional discrete Laplace
//    operator.
//
// 4) OMP execution of the X-Y sweep computation. This uses an array of temporaries and
//    operators to isolate data can code for each thread.
//
//
// Dependencies:

// omp.h
// GridFun2D.h, ClockIt.h, D2operator.h
//
// You will need to specify -fopenmp on the compilation command line
//
// Math 270A - UCLA
// Tue 07 Nov 2017 08:04:52 AM PST
//

#ifdef _OPENMP
#include <omp.h>              // required for OpenMP
#endif

#include "ClockIt.h"
#include "GridFun2D.h"
#include "D2operator.h"

int main()
{

    int threadCount = -1;
    cout << "Enter in number of threads: ";
    cin >> threadCount;

    if(threadCount > omp_get_max_threads()){threadCount = omp_get_max_threads();}
    if(threadCount <= 0)                   {threadCount = omp_get_max_threads();}
    omp_set_num_threads(threadCount);

    printf("\n");
    printf("#############\n");
    printf("############# Using OpenMP With %d Threads\n",omp_get_max_threads());
    printf("#############\n");
    printf("\n");

    long               N;
    long repetitions = 1;

    cout << "Enter panel count N:  ";
    cin >> N;

    cout << "Enter in number of repetitions: ";
    cin >> repetitions;

    long   xPanel  =     N;
    long   yPanel  =     N;

    double xMin   =     0.0;
    double xMax   =     1.0;
    double yMin   =     0.0;
    double yMax   =     1.0;

    double alphaX      =  2.0;  // Laplace operator-x prefactor
    double alphaY      =  2.0;  // Laplace operator-y prefactor
    double waveNumberX =  1.0;  // test problem x-coordinate wave number
    double waveNumberY =  1.0;  // test problem y-coordinate wave number

    //
    // Two methods to evaluate
    //
    // alphaX * d^2 U/dx^2  alphaY * d^2 U/dy^2
    //
    // with Dirichlet boundary conditions on a rectangular domain.
    // The standard 2nd order 5 point finite difference stencil is used to approximate
    // the second derivatives.
    //
    //

    GridFun2D     u(xPanel,xMin,xMax,yPanel,yMin,yMax);     // grid function to be differentiated

    GridFun2D      d2u(xPanel,xMin,xMax,yPanel,yMin,yMax);  // standard loop result

    GridFun2D    d2u_A(xPanel,xMin,xMax,yPanel,yMin,yMax);  // OpenMP + standard loop result

    GridFun2D   d2u_B1(xPanel,xMin,xMax,yPanel,yMin,yMax);  // X-Y sweep result

    GridFun2D   d2u_B2(xPanel,xMin,xMax,yPanel,yMin,yMax);  // OpenMP  + X-Y sweep result


    d2u.setToValue(0.0);
    d2u_A.setToValue(0.0);
    d2u_B1.setToValue(0.0);
    d2u_B2.setToValue(0.0);

    double hx = u.hx;
    double hy = u.hy;

    //
    // Loop index variables. Declared outside for statements
    // to facilitate use of OpenMP pragmas
    //

    long i; long j; long k; long m;

    //
    // Create function to be differentiated
    //

    double x;
    double y;
    double pi =  3.1415926535897932;

    for(i = 0; i <= xPanel; i++)
    {
        x = xMin + i*hx;
        for(j = 0; j <= yPanel; j++)
        {
            y =  yMin + j*hy;
            u.values(i,j)      = cos(2.0*pi*x*waveNumberX)*cos(2.0*pi*y*waveNumberY);
        }
    }

    //
    // XXXXX Computation # 1 XXXXX
    // Nested loop evaluation
    //

    ClockIt clock_1;

    clock_1.start();

for(k = 0; k < repetitions; k++)
{
    for(i = 1; i < xPanel; i++)
    {
    for(j = 1; j < yPanel; j++)
    {
        d2u.values(i,j) =
        alphaX*((u.values(i+1,j) - 2.0*u.values(i,j) + u.values(i-1,j))/(hx*hx))
        +
        alphaY*((u.values(i,j+1) - 2.0*u.values(i,j) + u.values(i,j-1))/(hy*hy));
    }
    }

} // repetitions loop end

    clock_1.stop();
    double averageTimeStdLoop = clock_1.getMilliSecElapsedTime()/(double)repetitions;

    // XXXXX Computation # 2 XXXXX
    //
    // Nested loop evaluation with inner loop multi-threaded
    //
    // You must insert appropriate OpenMP pragmas to complete the programming
    //

    ClockIt clock_2;
    clock_2.start();

    for(k = 0; k < repetitions; k++)
    {
    for(i = 1; i < xPanel; i++)
    {
    for(j = 1; j < yPanel; j++)
    {
        d2u_A.values(i,j) =
        alphaX*((u.values(i+1,j) - 2.0*u.values(i,j) + u.values(i-1,j))/(hx*hx))
        +
        alphaY*((u.values(i,j+1) - 2.0*u.values(i,j) + u.values(i,j-1))/(hy*hy));
    }
    }

}  // repetitions loop end

    clock_2.stop();
    double averageTimeOMPloop = clock_2.getMilliSecElapsedTime()/(double)repetitions;

    //
    // Evaluate the difference between procedures to verify correct computation
    //

    d2u_A -= d2u;

    double nestedLoopDiff = d2u_A.normInf();


    printf("==============================================\n");
    printf("Nested Loop Results using %ld repetitions \n",repetitions);
    printf("==============================================\n");
    printf("Nested loop values difference  : %10.6e \n\n",nestedLoopDiff);
    printf("Nested loop time          (ms) : %10.6f \n",averageTimeStdLoop);
    printf("Nested loop with OMP time (ms) : %10.6f \n\n",averageTimeOMPloop);
    printf("Ratio                           : %10.6f \n",averageTimeStdLoop/averageTimeOMPloop);

    //
    // XXXXX Computation # 3 XXXXX
    //
    // Using x-y sweeps and employing D2operator class
    //

    vector<double> uXvalues(xPanel+1,0.0);
    vector<double> uYvalues(yPanel+1,0.0);

    vector<double> d2uXvalues(xPanel+1,0.0);
    vector<double> d2uYvalues(yPanel+1,0.0);

    D2operator d2Xoperator(alphaX,hx);
    D2operator d2Yoperator(alphaY,hy);


    ClockIt clock_3;

    clock_3.start();

for(k = 0; k < repetitions; k++)
{

    //
    // X-sweep : computation of the d^2/dx^2 for each j
    //

    for(j = 1; j < yPanel; j++)
    {

        for(i = 0; i <= xPanel; i++)
        {
        uXvalues[i] = u.values(i,j);
        }

        d2Xoperator.apply(uXvalues,d2uXvalues);

        for(i = 0; i <= xPanel; i++)
        {
        d2u_B1.values(i,j) = d2uXvalues[i];
        }

    }

    //
    // Y-sweep : computation of the d^2/dy^2 for each i
    //

    for(i = 1; i < xPanel; i++)
    {

        for(j = 0; j <= yPanel; j++)
        {
        uYvalues[j] = u.values(i,j);
        }

        d2Yoperator.apply(uYvalues,d2uYvalues);

        for(j = 0; j <= yPanel; j++)
        {
        d2u_B1.values(i,j) += d2uYvalues[j]; // Note the  += here!
        }

    }

}  // repetitions loop end

    clock_3.stop();

    double averageTimeXYloop = clock_3.getMilliSecElapsedTime()/(double)repetitions;

    //
    // Evaluate the difference between procedures to verify correct computation
    //

    d2u_B1 -= d2u;
    double xyLoopDiff = d2u_B1.normInf();

    //
    // XXXXX Computation # 4 XXXXX
    //
    // Computation using x-y sweeps employing D2operator and OpenMP constructs.
    // To keep data and execution code for each thread separate,
    // temporaries and operators are created for each thread.

    //
    // Capture thread count (determines size of arrays of temporaries and operators)
    //

    threadCount = omp_get_max_threads();

    //
    // Allocate and initialize arrays of temporaries of size threadCount
    //

    vector< vector<double> >     uXvaluesOMP(threadCount);
    vector< vector<double> >     uYvaluesOMP(threadCount);
    vector< vector<double> >   d2uXvaluesOMP(threadCount);
    vector< vector<double> >   d2uYvaluesOMP(threadCount);

    for(m = 0; m < threadCount; m++)
    {
        uXvaluesOMP[m].resize(xPanel+1,0.0);
        uYvaluesOMP[m].resize(yPanel+1,0.0);

        d2uXvaluesOMP[m].resize(xPanel+1,0.0);
        d2uYvaluesOMP[m].resize(yPanel+1,0.0);
    }

    //
    // Allocate and initialize arrays of difference operators of size threadCount
    //

    vector< D2operator > d2XoperatorOMP(threadCount);
    vector< D2operator > d2YoperatorOMP(threadCount);

    for(m = 0; m < threadCount; m++)
    {
        d2XoperatorOMP[m].initialize(alphaX,hx);
        d2YoperatorOMP[m].initialize(alphaY,hy);
    }

    //
    // variable to hold the theadIndex
    //

    ClockIt clock_4;
    clock_4.start();

for(k = 0; k < repetitions; k++)
{

//
// X-sweep : Computation of the d^2/dx^2 for each j
//

#pragma omp parallel for private(i,j,m) schedule(static)

    for(j = 1; j < yPanel; j++)
    {
        m = omp_get_thread_num();

        for(i = 0; i <= xPanel; i++)
        {
        uXvaluesOMP[m][i] = u.values(i,j);
        }

        d2XoperatorOMP[m].apply(uXvaluesOMP[m],d2uXvaluesOMP[m]);

        for(i = 0; i <= xPanel; i++)
        {
            d2u_B2.values(i,j) = d2uXvaluesOMP[m][i];
        }
    }

//
// Y-sweep : computation of the d^2/dy^2 for each i
//

#pragma omp parallel for private(i,j,m) schedule(static)
    for(i = 1; i < xPanel; i++)
    {
        m = omp_get_thread_num();

        for(j = 0; j <= yPanel; j++)
        {
        uYvaluesOMP[m][j] = u.values(i,j);
        }

        d2YoperatorOMP[m].apply(uYvaluesOMP[m],d2uYvaluesOMP[m]);

        for(j = 0; j <= yPanel; j++)
        {
        d2u_B2.values(i,j) += d2uYvaluesOMP[m][j];  // Note the  += here!
        }

    }


}  // repetitions loop end


    clock_4.stop();
    double averageTimeXYompLoop = clock_4.getMilliSecElapsedTime()/(double)repetitions;

    d2u_B2 -= d2u;

    double xyOMPloopDiff = d2u_B2.normInf();


    printf("==============================================\n");
    printf("XY Sweep Results using %ld repetitions \n",repetitions);
    printf("==============================================\n");

    printf("X-Y sweep     values difference  : %10.6e \n",xyLoopDiff);
    printf("OMP X-Y sweep values difference  : %10.6e \n\n",xyOMPloopDiff);
    printf("X-Y sweep     loop time (ms) : %10.6f \n",averageTimeXYloop);
    printf("OMP X-Y sweep loop time (ms) : %10.6f \n\n",averageTimeXYompLoop);
    printf("Ratio                        : %10.6f \n",averageTimeXYloop/averageTimeXYompLoop);

    return 0;

}

