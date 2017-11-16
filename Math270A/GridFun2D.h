//
//  GridFun2D.h
//  Assignment4
//
//  Created by Kevin Miller on 10/30/17.
//  Copyright © 2017 Kevin Miller. All rights reserved.
//
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "DoubleArray2D.h"


using namespace std;


#ifdef _OPENMP
#include<omp.h>
#endif


#ifndef _GridFun2D_
#define _GridFun2D_

class GridFun2D {
    
public:
    Math270A::DoubleArray2D values;
    double hx;
    double xMin;
    double xMax;
    long xPanel;
    double hy;
    double yMin;
    double yMax;
    long yPanel;
    
    // Null Constructor
    GridFun2D(){
        initialize();
    };
    
    // Copy Constructor
    GridFun2D(const GridFun2D& f){
        initialize(f);
    };
    
    // Constructor, calls the initializer
    GridFun2D(long xPanel, double xMin, double xMax, long yPanel, double yMin, double yMax){
        initialize(xPanel, xMin, xMax, yPanel, yMin, yMax);
    };
    
    // Null Initializer
    void initialize(){
        xPanel = 0;
        xMin = 0;
        xMax = 0;
        hx = 0;
        yPanel = 0;
        yMin = 0;
        yMax = 0;
        hy = 0;
        //values.clear();
    };
    
    // Copy Initializer
    void initialize(const GridFun2D& f){
        xPanel = f.xPanel;
        xMin = f.xMin;
        xMax = f.xMax;
        hx = f.hx;
        yPanel = f.yPanel;
        yMin = f.yMin;
        yMax = f.yMax;
        hy = f.hy;
        values = f.values;
    };
    
    // Initializer
    void initialize(long xPanel, double xMin, double xMax, long yPanel, double yMin, double yMax){
        this->xPanel = xPanel;
        this->xMin = xMin;
        this->xMax = xMax;
        this->yPanel = yPanel;
        this->yMin = yMin;
        this->yMax = yMax;
        //this->values.resize(xPanel + 1);
        this->values.initialize(xPanel+1, yPanel+1);
        this->hx = (xMax - xMin)/float(xPanel);
        this->hy = (yMax - yMin)/float(yPanel);
    };
    
    // Destructor
    ~GridFun2D(){};
    
    // define = operator for GridFun2D objects
    void operator=(const GridFun2D& g){
        initialize(g);
    };
    


    #ifndef _OPENMP
    // define += operator for GridFun2D objects
    void operator+=(const GridFun2D& g){
        for(long i=0; i < values.getIndex1Size(); i++){
            for(long j=0; j < values.getIndex2Size(); j++){
                values(i, j) += g.values(i, j);
            }
        }
    };
    
    #else
    // define += operator for GridFun2D objects (WITH OPENMP)
    void operator+=(const GridFun2D& g){
        for (long i=0; i<values.getIndex1Size(); i++){
#pragma omp parallel for default(shared) private(j) schedule(static)
            for(long j=0; j<values.getIndex2Size(); j++){
                values(i,j) += g.values(i,j);
            }
        }
    };
    #endif
    
    // define -= operator for GridFun2D objects
    void operator-=(const GridFun2D& g){
        for(long i=0; i < values.getIndex1Size(); i++){
            for(long j=0; j < values.getIndex2Size(); j++){
                values(i, j) -= g.values(i, j);
            }
        }
    };
    
    // define *= operator with constants for GridFun2D objects
    void operator*=(double alpha){
        for(long i=0; i < values.getIndex1Size(); i++){
            for(long j=0; j < values.getIndex2Size(); j++){
                values(i, j) *= alpha;
            }
        }
    };
    
    // define /= operator with constants for GridFun2D objects
    void operator/=(double alpha){
        for(long i=0; i < values.getIndex1Size(); i++){
            for(long j=0; j < values.getIndex2Size(); j++){
                values(i, j) /= alpha;
            }
        }
    };
    
    // setToValue
    void setToValue(double d){
        for(long i=0; i < values.getIndex1Size(); i++){
            for(long j=0; j < values.getIndex2Size(); j++){
                values(i, j) = d;
            }
        }
    };
    
    // normInf function -- to calculate the infinity norm of values
    double normInf(){
        double norm_val = 0.0;
        for(long i=0; i < values.getIndex1Size(); i++){
            for(long j=0; j < values.getIndex2Size(); j++){
                if(abs(values(i,j)) > norm_val){
                    norm_val = abs(values(i,j));
                }
            }
        }
        return norm_val;
    };
    
    // define << operator for outStream
    friend ostream& operator<<(ostream& outStream, const GridFun2D& v){
        outStream << "Grid Function 2D values:" << endl;
        outStream << "x\t\ty\t\tf(x):" << endl;
        double x_0 = v.xMin;
        double h_x = v.hx;
        double y_0 = v.yMin;
        double h_y = v.hy;
        for(long i=0; i < v.values.getIndex1Size(); i++){
            for(long j=0; j < v.values.getIndex2Size(); j++){
                outStream << setw(5) << x_0 + i*h_x << "\t" << y_0 + j*h_y << "\t" << v.values(i,j) << endl;
            }
            cout << endl;
        }
        return outStream;
    };
};



#endif /* _GridFun2D_ */
