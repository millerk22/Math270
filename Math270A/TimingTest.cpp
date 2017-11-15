#include <cstdio>
#include <iostream>
#include <cmath>
using namespace std;
//
// TimingTest.cpp : A test code to time the execution of the +=
// operations in the GridFun2D class and the evaluation of a product
// of cosines at each of the grid nodes.
//
// Math 270A : Version Tue 07 Nov 2017 07:47:21 AM PST

#include "ClockIt.h"
#include "GridFun2D.h"

int main()
{
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

    }

    clock_2.stop();

    double averageCosEvaluationTime = clock_2.getMilliSecElapsedTime()/(double)repetitions;

    //
    // Print out the results and report the ratio
    //

    double totalTime = clock_1.getMilliSecElapsedTime() + clock_2.getMilliSecElapsedTime();
    double ratio     = averageCosEvaluationTime/averageTimePlusEquals;

    printf("==============================================\n");
    printf("Results using %ld repetitions \n",repetitions);
    printf("==============================================\n");
    printf("Total Time (ms) : %10.6f \n",totalTime);
    printf("Time  +=        operation (ms) : %10.6f \n",averageTimePlusEquals);
    printf("Time  cos evaluation      (ms) : %10.6f \n",averageCosEvaluationTime);

    printf("Ratio : %10.3f \n\n",ratio);

}
