#ifndef _D2operator_
#define _D2operator_
//
// D2operator.h
//
// Class D2operator's apply member function uses the standard
// 3-point discrete approximation to create an approximation to 
//
//  alpha* d^2 U/ dx^2 
//
// The mesh size and parameter alpha must be specified using
// the constructor or initialize member functions before the
// apply(..) operator is used. 
//
// Math 270A UCLA -  Tue 07 Nov 2017 12:56:00 PM PST
//

#include <vector>
using namespace std;

class D2operator
{
public:
    D2operator()
    {
        alpha    = 0.0;
        h        = 0.0;
    }

    D2operator(const D2operator& S)
    {
        alpha = S.alpha;
        h     = S.h;
    }

    D2operator(double alpha, double h)
    {
        this->alpha    = alpha;
        this->h        = h;
    }


    void initialize(const D2operator& S)
    {
        alpha = S.alpha;
        h     = S.h;
    }

    void initialize(double alpha, double h)
    {
        this->alpha    = alpha;
        this->h        = h;
    }

    //
    // The apply member function evaluates the 3-point discrete
    // Laplacian with coefficient alpha at all "interior" points
    // and returns the result in d2u.
    //
    // The first and last element of d2u are set to 0 and u is
    // unchanged by this member function.
    //
    void apply(const vector<double>& u, vector<double>& d2u)
    {

    long panels = u.size() - 1;
    d2u[0]      = 0.0;

    for(long i = 1; i < panels; i++)
    {
        d2u[i] = alpha*(((u[i+1] - 2.0*u[i] + u[i-1]))/(h*h));
    }

    d2u[panels] = 0.0;
    }

    double alpha;
    double h;
};

#endif
