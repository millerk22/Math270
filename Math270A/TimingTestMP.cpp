#include <cstdio>
#include <iostream>
#include <cmath>
using namespace std;
//
// TimingTestMP.cpp : A test code to time the execution of the  +=
// operations in the GridFun2D class and the evaluation of a product
// of cosines at each of the grid nodes.
//
// This program also times the execution of the evaluation of
// the product of cosines using OpenMP constructs.
//
//
// Math 270: Version  Tue 07 Nov 2017 08:04:52 AM PST
//
#include "ClockIt.h"
#include "GridFun2D.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int main()
{
    //
    // Specifying number of threads to use and checking that code was compiled with OpenMP enabled.
    //

    int threadCount = -1;


#ifdef _OPENMP
    cout << "Enter in number of threads: ";
    cin >> threadCount;
#else
    cout << "Program not compiled with OpenMP enabled. " << endl;
    cout << "Recompile with -fopenmp flag specified " << endl;
    exit(0);
#endif

    if(threadCount > omp_get_max_threads()){threadCount = omp_get_max_threads();}
    if(threadCount <= 0)                   {threadCount = omp_get_max_threads();}
    omp_set_num_threads(threadCount);

    printf("\n");
    printf("#############\n");
    printf("############# Using OpenMP With %d Threads\n",omp_get_max_threads());
    printf("#############\n");
    printf("\n");

    long   xPanel;
    long   yPanel;

    double xMin   =     0.0;
    double xMax   =     1.0;

    double yMin   =     0.0;
    double yMax   =     1.0;

    double waveNumberX =  1.0;  // test problem x-coordinate wave number
    double waveNumberY =  1.0;  // test problem y-coordinate wave number

    long N;
    long repetitions;

    cout << "Enter panel count N:  ";
    cin >> N;

    xPanel = N;
    yPanel = N;

    cout << "Enter in number of repetitions: ";
    cin >> repetitions;

    //
    // Allocate GridFun2D;
    //

    GridFun2D aFun(xPanel,xMin,xMax,yPanel,yMin,yMax);
    GridFun2D bFun(xPanel,xMin,xMax,yPanel,yMin,yMax);

    aFun.setToValue(0.0);
    bFun.setToValue(0.0);

    //
    // Time the += operator of the GridFun2D class.
    //

    ClockIt clock_1;

    clock_1.start();

    for(long k = 0; k < repetitions; k++)
    {
    	aFun += bFun;
    }

    clock_1.stop();

    double averageTimePlusEquals = clock_1.getMilliSecElapsedTime()/(double)repetitions;

    //
    // Time the evaluation of the product of cos functions at grid nodes
    //

    double x; double y;
    double pi =  3.1415926535897932;
    long i; long j;

    double hx = aFun.hx;
    double hy = aFun.hy;

    ClockIt clock_2;

    clock_2.start();

    for(long k = 0; k < repetitions; k++)
    {

    	for(i = 0; i <= xPanel; i++)
    	{
    		x = xMin + i*hx;
    		for(j = 0; j <= yPanel; j++)
    		{
    			y =  yMin + j*hy;
    			aFun.values(i,j)  = cos(2.0*pi*x*waveNumberX)*cos(2.0*pi*y*waveNumberY);
    		}
    }

    }  // repetition end

    clock_2.stop();

    double averageCosEvaluationTime = clock_2.getMilliSecElapsedTime()/(double)repetitions;

    //
    // Time the evaluation of the product of cos functions at grid nodes
    // with evaluation being performed with a multi-threaded loop.
    //


    ClockIt clock_3;

    clock_3.start();

    for(long k = 0; k < repetitions; k++)
    {

    for(i = 0; i <= xPanel; i++)
    {
        x = xMin + i*hx;
		#pragma omp parallel for default(shared) private(j,y) schedule(static)
        for(j = 0; j <= yPanel; j++)
        {
            y =  yMin + j*hy;
            bFun.values(i,j)  = cos(2.0*pi*x*waveNumberX)*cos(2.0*pi*y*waveNumberY);
        }

    }

    } // repetition end

    clock_3.stop();
    double averageOMPcosEvaluationTime = clock_3.getMilliSecElapsedTime()/(double)repetitions;

    //
    // Compute difference to verify accuracy
    //

    bFun -= aFun;
    double computationDifference = bFun.normInf();

    cout << endl;
    cout << "OMP vs. non-OMP evaluation value difference : " << computationDifference << endl;
    cout << endl;

    //
    // Print out the results and report the ratio
    //

    double totalTime = clock_1.getMilliSecElapsedTime() + clock_2.getMilliSecElapsedTime();
    double ratio     = averageCosEvaluationTime/averageTimePlusEquals;
    double ratio_2   = averageCosEvaluationTime/averageOMPcosEvaluationTime;

    printf("==============================================\n");
    printf("Results using %ld repetitions \n",repetitions);
    printf("==============================================\n");
    printf("Total Time (ms) : %10.6f \n",totalTime);
    printf("Cos computation difference : %10.5e \n",computationDifference);
    printf("Time  +=        operation (ms) : %10.6f \n",averageTimePlusEquals);
    printf("Time  cos evaluation      (ms) : %10.6f \n",averageCosEvaluationTime);
    printf("Time  cos OMP evaluation  (ms) : %10.6f \n\n",averageOMPcosEvaluationTime);


    printf("Ratio cos to +=     : %10.3f \n",ratio);
    printf("Ratio cos to OM cos : %10.3f \n\n",ratio_2);

}
