// calculates the one-dimensional poisson equation
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <armadillo>
#include <time.h>

using namespace std;
using namespace arma;

double func (double, double);
double solution(double);
double error (double, unsigned int, double *);
void tridiag (double, unsigned int, double *);
void tridiag_poisson (double, unsigned int, double *); 
vec lu_armadillo (double, unsigned int);

int main()
{   
    clock_t time_tri, time_lu;                      // calculate time         //
    unsigned int n, i;                              // number of steps, index //
    double x_min, x_max;                            // x-boundaries           //
    double *u, u_min, u_max;                        // num. sol., bound.      //
    double h_step, err;                             // stepwidth, error       //  
    vec v;                                          // solution of lu-decomp. //
    
    n = 10;                                         // number of steps        //
    x_min = 0.;                                     // area of calculation    //
    x_max = 1.;                                     //                        //
    h_step = (x_max - x_min) / (n + 1);             // stepsize               // 

    u_min = 0.;                                     // dirichlet-boundaries   //
    u_max = 0.;                                     //                        //
    
    u = new double[n+2];                            // initialize array for   //
    
    u[0] = u_min;                                   // solution u             //
    u[n+1] = u_max;                                 //                        // 

    // tridiag solves the equation Au = f
    // tridiagonal poisson code
    time_tri = clock() ;
    tridiag_poisson (h_step, n, u);                 // own algorithm          //
    //tridiag (h_step, n, u);                       // thomas algorithm       //
    time_tri = clock() - time_tri;    

    // general lu decomposition code
    time_lu = clock();
    v = lu_armadillo(h_step, n);                    // lu decomposition       //
    time_lu = clock() - time_lu;

    err = error(h_step, n, u);                      // error-evalution        //

    // output
    cout << "--------------------------------------------" << endl;
    cout << "Maximum Error: " << err << endl;
    cout << "Calculation time / ms (Own Algorithm): " << time_tri << endl;
    cout << "Calculation time / ms (Armadillo LU Dec.): "  << time_lu << endl;
    cout << endl;
    cout << "--------------------------------------------" << endl;
    cout << "x" << "    " << "u" << "    " << "v" << endl;
    cout << h_step*0 << "    " << u[0] << "    " << v[n] << endl; 
    for (i = 1; i <= n + 1; i++)
    {
        cout << h_step*i << "    " << u[i] << "    " << v[i-1] << endl;
    }
    
    delete [] u;
    return 0;
}

void tridiag (double h_step, unsigned int n, double *u){
    /* Uses the Thomas-Algorithm, code from "Numerical Recipes,               *  
     * Third Edition", p. 56 f.                                               */
    double a,b,c;
    double btemp;
    unsigned int i; 
    vec temp(n);

    a = -1;
    c = a ;
    b = 2 ; 

    btemp = b;
    u[1] = func(h_step, h_step*1.)/btemp;
    for(i = 2 ; i <= n; i++){
        temp[i] = c/btemp;                                  // 1 FLOP
        btemp = b - a*temp[i];                              // 2 FLOPs
        u[i] = (func(h_step, h_step*i) - a*u[i-1])/btemp;   // 3 FLOPs
    } 

    for(i = n; i >= 1; i--){
        u[i] -= temp[i+1]*u[i+1];                           // 2 FLOPs
    }                                                       // -> 8N FLOPs
}

void tridiag_poisson (double h_step, unsigned int n, double *u){
    /* Just usable for poisson equation tridiagonal matrices. Reduces number  *
     * of FLOPs to 6N.                                                        */
    unsigned int i;
    vec ftemp(n+1);

    ftemp[1] = func(h_step, h_step*1);
    for (i = 2; i <= n; i++){
        ftemp[i] = ftemp[i-1] + i*func(h_step, h_step*i);   // 2 FLOPs
    }

    for (i = n; i >= 1; i--){
        u[i] = (ftemp[i] + u[i+1]*i)/(1.+i);                // 4 FLOPs
    }                                                       // -> 6 FLOPs
}

double error (double h_step, unsigned int n, double *u){
    /* Calculates the logarithmic error bewteen numerical and analytical      *
     * solution.                                                              */
    unsigned int i;
    double eps, eps_max, u_ana;
    eps_max = -100.0; 

    for (i = 1; i <= n; i++){
        u_ana = solution(h_step*i);
        eps = log10(fabs((u[i] - u_ana)/u_ana));
        if (eps > eps_max){
            eps_max = eps;
        }
    }
    return eps_max;
}

double func (double h_step, double x)
{
    // second derivative rhs-function 
    double f;

    f = h_step*h_step*100.*exp(-10.*x);
    return f; 
}

double solution (double x)
{
    // analytical solution 
    double f;
    f = 1.-(1.-exp(-10.))*x-exp(-10.*x);
    return f; 
}

vec lu_armadillo(double h_step, unsigned int n)
{
    /* lu-decomposition to solve the same Equation than the tridiag           *
     * algorithms but width a worse performance                               */
    unsigned int i;

    vec f(n);                                       // Ax = f                 //
    mat A(n,n);                                     //                        // 
    vec y, x;                                       // Ax = LUx = Ly = f      //
    mat L, U, P;                                    //                        //
    

    // define matrix A for as in poisson equation
    f(0) = func(h_step, h_step*(1));
    A(0,0) = 2.0;
    for (i = 1; i < n; i++){
        A(i,i) = 2.0;
        A(i,i-1) = -1.0;
        A(i-1,i) = -1.0;
        f(i) = func(h_step, h_step*(i+1));
    }

    lu(L, U, A);                                    // lu-decomposition       //
    solve(y, L, f);                                 // solve Ly = f           //
    solve(x, U, y);                                 // solve Ux = y           //
    return x; 
}
