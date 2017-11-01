//
//  GridFun1D.h
//  Assignment2
//
//  Created by Kevin Miller on 10/17/17.
//  Copyright © 2017 Kevin Miller. All rights reserved.
//

#include <iostream>   // For C++ input/output
#include <iomanip>
#include <vector>     // For STL vectors
#include <cstdio>     // For C input/output
#include <cmath>      // For math functions
#include <cstdlib>
using namespace std; // Using the "standard" (std) standard template library components

#ifndef _GridFun1D_
#define _GridFun1D_

class GridFun1D {
    
public:
    double hx;
    double xMin;
    double xMax;
    long xPanel;
    vector<double> values;
    
    // Null Constructor
    GridFun1D(){
        initialize();
    };
    
    // Copy Constructor
    GridFun1D(const GridFun1D& f){
        initialize(f);
    };
    
    // Constructor, calls the initializer
    GridFun1D(long xPanel, double xMin, double xMax){
        initialize(xPanel, xMin, xMax);
    };
    
    // Null Initializer
    void initialize(){
        xPanel = 0;
        xMin = 0;
        xMax = 0;
        values.clear();
    };
    
    // Copy Initializer
    void initialize(const GridFun1D& f){
        xPanel = f.xPanel;
        xMin = f.xMin;
        xMax = f.xMax;
        values = f.values;
    };
    
    // Initializer
    void initialize(long xPanel, double xMin, double xMax){
        this->xPanel = xPanel;
        this->xMin = xMin;
        this->xMax = xMax;
        this->values.resize(xPanel + 1);
        this->hx = (xMax - xMin)/float(xPanel);
    };
    
    // define = operator for GridFun1D objects
    void operator=(const GridFun1D& g){
        initialize(g);
    };
    
    // define += operator for GridFun1D objects
    void operator+=(const GridFun1D& g){
        for(long i=0; i < values.size(); i++){
            values[i] += g.values[i];
        }
    };
    
    // define -= operator for GridFun1D objects
    void operator-=(const GridFun1D& g){
        for(long i=0; i < values.size(); i++){
            values[i] -= g.values[i];
        }
    };
    
    // define *= operator with constants for GridFun1D objects
    void operator*=(double alpha){
        for(long i=0; i < values.size(); i++){
            values[i] *= alpha;
        }
    };
    
    // define /= operator with constants for GridFun1D objects
    void operator/=(double alpha){
        for(long i=0; i < values.size(); i++){
            values[i] /= alpha;
        }
    };
    
    // setToValue
    void setToValue(double d){
        for(long i=0; i < values.size(); i++){
            values[i] = d;
        }
    };
    
    // normInf function -- to calculate the infinity norm of values
    double normInf(){
        double norm_val = 0.0;
        for(long i = 0; i < values.size(); i++){
            if(abs(values[i]) > norm_val){
                norm_val = abs(values[i]);
            }
        }
        return norm_val;
    };
    
    // define << operator for outStream
    friend ostream& operator<<(ostream& outStream, const GridFun1D& v){
        outStream << "Grid Function 1D values:" << endl;
        outStream << "x\t\tf(x):" << endl;
        double x_0 = v.xMin;
        double h_x = v.hx;
        for(long i=0; i < v.values.size(); i++){
            outStream << setw(5) << x_0 + i*h_x << "\t" << v.values[i] << endl;
        }
        return outStream;
    };
    
    
};


#endif /* _GridFun1D_ */

