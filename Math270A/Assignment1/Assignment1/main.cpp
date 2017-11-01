#include <iostream>   // For C++ input/output
#include <vector>     // For STL vectors
#include <cstdio>     // For C input/output
#include <cmath>      // For math functions
using namespace std; // Using the "standard" (std) standard template library components


/**
 Instances of class TriOperator are linear operators whose apply member function
 multiplies a tri-diagonal matrix with an input vector.
 
 The tri-diagonal matrix associated with the class instance is
 specified by it's lower diagonal, diagonal and upper diagonal entries, e.g.
 
 n X n tridiagonal matrix T where n = systemSize
 
 
 T =  |  diag[0]  upDiag[0]                                  |
 |loDiag[0]    diag[1] upDiag[1]                        |
 |           loDiag[1]                                  |
 |                *       *       *                     |
 |                   *         *        *               |
 |                loDiag[n-3]   diag[n-2]  upDiag[n-2]  |
 |                             loDiag[n-2]   diag[n-1]  |
 
 
 
 Currently, STL (Standard Template Library) vectors are used to store
 the diagonals that define the non-zero values of T
 
 Version: Friday Sept 29 2017 01:50:45 PM PST
 */
class TriOperator
{
public:
    
    /// Null Constructor -- called when you declare an instance
    
    TriOperator()
    {initialize();};
    
    /// Copy constructor (creates duplicate of T)
    
    TriOperator(const TriOperator& T) // called when you declare an instance with an existing instance
    {initialize(T);};
    
    
    /**
     This constructor initializes the internal class data with entries
     of the n X n tridiagonal matrix T associated with the entries in
     the lower diagonal (loDiag), diagonal (diag), and upper diagonal (upDiag).
     */
    
    TriOperator(long systemSize, vector<double>& loDiag, vector<double>& diag,
                vector<double>& upDiag)
    {
        initialize(systemSize,loDiag,diag,upDiag);
    }
    
    // Default destructor
    
    virtual ~TriOperator(){};
    
    /// Null initializer, resizes the internal arrays to zero length
    
    void initialize()
    {
        systemSize = 0;
        loDiag.clear();
        diag.clear();
        upDiag.clear();
        
    };
    
    ///  Copy initializer. Duplicates the entries of T
    
    void initialize(const TriOperator& T)
    {
        systemSize = T.systemSize;
        loDiag     = T.loDiag;
        diag       = T.diag;
        upDiag     = T.upDiag;
    };
    
    
    /**
     This initializer initializes the internal class data with entries
     of the n X n tridiagonal matrix T associated with the entries in
     the lower diagonal (loDiag), diagonal (diag), and upper diagonal (upDiag).
     */
    
    void initialize(long systemSize, vector<double>& loDiag,
                    vector<double>& diag, vector<double>& upDiag)
    {
        this->systemSize = systemSize;
        this->loDiag     = loDiag;
        this->diag       = diag;
        this->upDiag     = upDiag;
    }
    
    /**
     The apply operator evaluates
     
     vOut = T*vIn
     
     where T is the tri-diagonal operator specified by the class constructor
     or initializer.
     
     Currently no bounds checking is performed.
     
     */
    
    
    void apply(vector<double>& vIn, vector<double>& vOut)
    {
        long i;
        
        if(systemSize == 1)
        {
            vOut[0] = diag[0]*vIn[0]; return;
        }
        
        vOut[0] = diag[0]*vIn[0] + upDiag[0]*vIn[1];
        
        for(i = 1; i < systemSize-1; i++)
        {
            vOut[i] = loDiag[i-1]*vIn[i-1] + diag[i]*vIn[i] + upDiag[i]*vIn[i+1];
        }
        
        i = systemSize-1;
        vOut[i] = loDiag[i-1]*vIn[i-1] + diag[i]*vIn[i];
    }
    
    // Internal data
    
    
    long        systemSize;
    vector<double>  loDiag;
    vector<double>    diag;
    vector<double>  upDiag;
    
};

// A specification of TriSolver class with "nothing in it". It's just
// a temporary specification with the absolutely minimum number
// of member functions implemented that will allow this
// program to comM_PIle and run.
//
// Since the apply member function of this
// version is just the identity, the errors in the
// computed solution of the linear system will be non-trivial.
//
// You will need to replace this class with your version of the
// TriSolver class as described in Assignment 1.
//

class TriSolver
{
public:
    
    // Null Constructor
    TriSolver()
    { initialize();};
    
    TriSolver(const TriSolver& T) // called when you declare an instance with an existing instance
    {initialize(T);};
    
    // Null initializer. Set the
    
    void initialize(){
        systemSize = 0;
        L_loDiag.clear();
        L_diag.clear();
        U_upDiag.clear();
    };
    
    ///  Copy initializer. Duplicates the entries of T
    
    void initialize(const TriSolver& T)
    {
        systemSize = T.systemSize;
        L_loDiag = T.L_loDiag;
        L_diag = T.L_diag;
        U_upDiag = T.U_upDiag;
    };
    
    
    // Constructor -- Here we solve for the diagonals of the LU factorization. Thus, when we apply this solver for different boundary conditions (different f vectors), we don't have to refactor into LU factorization.
    
    void initialize(long systemSize, vector<double>& loDiag, vector<double>& diag, vector<double>& upDiag)
    {
        this->systemSize = systemSize;
        
        // Compute LU factorization, as discussed in class
        
        this->L_loDiag = loDiag; // the subdiagonal of L is just the subdiagonal of the original matrix
        
        L_diag.resize(systemSize);
        U_upDiag.resize(systemSize-1);
        
        L_diag[0] = diag[0];
        
        U_upDiag[0] = upDiag[0]/L_diag[0];
        
        for (long i=1; i <= systemSize-2; i ++){
            L_diag[i] = diag[i] - loDiag[i-1]*U_upDiag[i-1]; // compute main diagonal of L
            U_upDiag[i] = upDiag[i]/L_diag[i]; // compute the super diagonal of U
        }
        
        L_diag[systemSize-1] = diag[systemSize-1] - loDiag[systemSize-2]*U_upDiag[systemSize-2]; // compute last entry in the main diagonal of L
        
        this->L_diag = L_diag; // assign the calculated diagonals to be parts of our Solver.
        this->U_upDiag = U_upDiag;
        
        return;
    }
    
    // Default destructor
    virtual ~TriSolver(){};
    
    
    // given the vector f, we calculate a solution u.
    
    void apply(long systemSize, vector<double>& f, vector<double>& u)
    {
        // Have the LU factorization A = LU computed in constructor
        
        // Step 1 -- Solve the system Lv = f
        vector<double> v; // vector we use as intermediate step of solving this system Lv = f
        v.resize(systemSize);
        
        v[0] = f[0] / L_diag[0];
        for(long i=1; i <= systemSize-1; i ++){
            v[i] = (f[i] - L_loDiag[i-1]*v[i-1])/L_diag[i];
        }
        
        // Step 2 -- Solve the system Uu = v
        u[systemSize-1] = v[systemSize-1];
        for(long i=systemSize-2; i >= 0; i--){
            u[i] = v[i] - U_upDiag[i]*u[i+1];
        }
        
        // delete v, clean up space
        v.clear();
        
        return;
    }
    
    // internal data
    long systemSize ;
    vector<double> L_loDiag;
    vector<double> L_diag;
    vector<double> U_upDiag;
};

//
// Local utility function to set up the coefficients for the
// discrete system.
//
void setUpSystem(double h, double alpha, double beta, long M,
                 vector<double>& loDiag, vector<double>& diag, vector<double>& upDiag)
{
    loDiag.resize(M);
    upDiag.resize(M);
    diag.resize(M+1);
    
    long i;
    
    // First equation is the identity
    
    i = 0;
    diag[0]   =  1.0;
    upDiag[0] =  0.0;
    
    // Interior equation coefficients associated with a standard
    // second order finite difference approximation
    
    for(long i = 1; i < M; i++)
    {
        loDiag[i-1] =  alpha/(h*h);
        upDiag[i]   =  alpha/(h*h);
        diag[i]     = -2.0*alpha/(h*h) + beta;
    }
    
    
    // Last equation is the identity
    
    i          =    M;
    diag[M]    =  1.0;
    loDiag[M-1] = 0.0;
}


// Exact solution for problem (a)

double Exact_a(double x)
{
    return 1.0 + x;
}

// Exact solution for problem (b), knowing that alpha = 2.0
double Exact_b(double x){
    return -cos(2.0*M_PI*x)/(8.0*M_PI*M_PI);
}

// Exact solution for problem (c), knowing that alpha = 2.0, and beta = -1.0
double Exact_c(double x){
    return -cos(2.0*M_PI*x)/(1.0 + 8.0*M_PI*M_PI);
}



int main()
{

    // Store the different values for each test we will be doing. Then, we will just loop through these different parameters and print out the results at the end.
    
    vector<double> leftVals = {1.0,-1.0/(8.0*M_PI*M_PI), -1.0/(1.0 + 8.0*M_PI*M_PI)};
    vector<double> rightVals = {2.0,-1.0/(8.0*M_PI*M_PI), -1.0/(1.0 + 8.0*M_PI*M_PI)};
    vector<long> numPanels = {10, 20, 40};
    vector<double> Betas = {0.0, 0.0, -1.0};
    
    vector<double> solError2_hist {0.0, 0.0, 0.0}; // vector to store the error vals for order of convergence calculation with 2 norm
    
    double alpha  =  2.0; // coefficients
    double xMin   =  0.0; // limits of the domain
    double xMax   =  1.0;
    
    string prob; // string variable to be used for printing
    
    for (int k=0; k <= 2; k++){ // Loop for the different problems
        double beta   =  Betas[k]; // coefficient beta
        
        for (int j=0; j <= 2; j++){ // Loop for the different mesh sizes (number of panels)
            TriOperator tOperator;  // instantiate a forward operator --- calls null constructor
            TriSolver     tSolver;  // instantiate a forward operator inverse operator
            
            long M = numPanels[j]; // constants dependent on our problem
            double h = (xMax-xMin)/(double)M;
            
            
            vector<double> loDiag;   // instantiating arrays for coefficients
            vector<double> upDiag;
            vector<double> diag;
            
            
            vector<double> f(M+1,0.0);      // declaring solution and right hand side vectors
            vector<double> u(M+1,0.0);      // and initializing to zero.
            vector<double> fStar(M+1,0.0);
            
            
            
            // Assign boundary values to first and last element of right hand side
            f[0] =  leftVals[k];
            f[M] =  rightVals[k];
            
            if (k == 0){ // if problem (a), then f = 0 for all interior pts
                for(long i = 1; i < M; i++)
                {
                    f[i] = 0.0;
                }
            }
            else{ // if not problem (a), then f = cos(2*M_PI*x)
                double x;
                for(long i=1; i < M; i++){
                    x = xMin + h*i;
                    f[i] = cos(2.0*M_PI*x);
                    
                }
            }
            
            
            // Setting up system
            
            setUpSystem(h, alpha, beta, M, loDiag, diag, upDiag);
            
            // Initializing the forward operator
            
            tOperator.initialize(M+1,loDiag,diag,upDiag);
            
            // Initializing the inverse operator
            
            tSolver.initialize(M+1,loDiag,diag,upDiag);
            
            // Applying inverse
            
            tSolver.apply(M+1,f,u); // Note that systemSize will be M+1
            
            // Applying forward operator to the result (to evaluate the residual)
            
            tOperator.apply(u,fStar);
            
            
            // Determining the size of the residual and the size of the error
            
            double aErrorMax = 0.0; // computed error, inf norm
            double fErrorMax = 0.0; // residual error, inf norm
            double aError2 = 0.0; // computed error, 2 norm
            double fError2 = 0.0; // residual error, 2 norm
            double errVal;
            double exact;
            double x;
            
            // based on k loop, we calculate the errors compared to exact soln as well as residual error
            if (k == 0){
                for (int i =0; i <= M; i++){
                    x = xMin + i*h;
                    exact = Exact_a(x); // calculate exact soln at the point x
                    
                    // calculate inf norms
                    errVal = fabs(u[i] - exact);
                    if(aErrorMax < errVal) {aErrorMax = errVal;}
                    errVal = fabs(fStar[i] - f[i]);
                    if(fErrorMax < errVal) {fErrorMax = errVal;}
                    
                    // calculate 2 norms
                    aError2 += (u[i] - exact)*(u[i] - exact);
                    fError2 += (fStar[i] - f[i])*(fStar[i] - f[i]);
                }
            }
            else if (k == 1){
                for (int i =0; i <= M; i++){
                    x = xMin + i*h;
                    exact = Exact_b(x); // calculate exact soln at the point x
                    
                    // calculate inf norms
                    errVal = fabs(u[i] - exact);
                    if(aErrorMax < errVal) {aErrorMax = errVal;}
                    errVal = fabs(fStar[i] - f[i]);
                    if(fErrorMax < errVal) {fErrorMax = errVal;}
                    
                    // calculate 2 norms
                    aError2 += (u[i] - exact)*(u[i] - exact);
                    fError2 += (fStar[i] - f[i])*(fStar[i] - f[i]);
                }
            }
            else{
                for (int i =0; i <= M; i++){
                    x = xMin + i*h;
                    exact = Exact_c(x); // calculate exact soln at the point x
                    
                    // calculate inf norms
                    errVal = fabs(u[i] - exact);
                    if(aErrorMax < errVal) {aErrorMax = errVal;}
                    errVal = fabs(fStar[i] - f[i]);
                    if(fErrorMax < errVal) {fErrorMax = errVal;}
                    
                    // calculate 2 norms
                    aError2 += (u[i] - exact)*(u[i] - exact);
                    fError2 += (fStar[i] - f[i])*(fStar[i] - f[i]);
                }
            }
            
            
            aError2 = sqrt(h)*sqrt(aError2); // Scaling for 2 norms
            fError2 = sqrt(h)*sqrt(fError2);
            
            
            cout.setf(ios::scientific); // Scientific notation
            cout.precision(15);         // 15 digits
            
            
            if (k == 0){
                prob = "(a)";
            }
            else if(k == 1){
                prob = "(b)";
            }
            else{
                prob = "(c)";
            }
            
            cout << "Problem " << prob << "  with panel size "<< M << ":" << endl << endl;
            cout << "Infinity norm:" << endl;
            cout << "    Residual error " << fErrorMax << endl;
            cout << "    Solution error " << aErrorMax << endl;
            cout << "2 norm:" << endl;
            cout << "    Residual error " << fError2 << endl;
            cout << "    Solution error " << aError2 << endl;
            cout << endl;
            
            // store the solution error in 2 norm for order of convergence calculation at end
            solError2_hist[j] = aError2;
        }
        cout << endl;
        
        // Order of convergence calculation
        
        double conv_ord1 = log(solError2_hist[0]/solError2_hist[1])/log(2.0);
        double conv_ord2 = log(solError2_hist[1]/solError2_hist[2])/log(2.0);
        cout << "Order of Convergence for " << prob << ":" << endl;
        cout << "\t"<< numPanels[0] << " --> "<< numPanels[1] << " panels : " << conv_ord1 << endl;
        cout << "\t"<< numPanels[1] << " --> "<< numPanels[2] << " panels : " << conv_ord2 << endl;
        cout << endl << "____________________________________________________________" << endl << endl;
        
        
        
    }
    
    return 0;
}
