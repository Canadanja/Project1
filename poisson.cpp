// calculates the one-dimensional poisson equation
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <armadillo>

using namespace std;
using namespace arma;

double func (double, double);
double solution(double);
double error (double, unsigned int, double *);
void tridiag (double, unsigned int, double *); 
void tridiag_poisson (double, unsigned int, double *); 



int main()
{   
    unsigned int n, i;                              
    double *x, x_min, x_max;                        
    double *u, u_min, u_max;                      
    double h_step, err;                               
    
    n = 10;                                         // number of steps        //
    x_min = 0.;                                     // area of calculation    //
    x_max = 1.;                                     //                        //
    x = new double[n+2];                            // initial value for var. // 

    u_min = 0.;                                     // dirichlet-boundaries   //
    u_max = 0.;                                     //                        //
    h_step = (x_max - x_min) / (n + 1);             // stepsize               // 
    
    u = new double[n+2];                            // initialize array for   //
    
    u[0] = u_min;                                   // solution u             //
    u[n+1] = u_max;                                 //                        // 
    
    // tridiag solves the equation Au = f
    tridiag_poisson (h_step, n, u); 
    
    err = error(h_step, n, u);

    // output
    cout << "Maximum Error: " << err << endl;
    cout << endl;
    cout << "---------------------------------" << endl;
    cout << "x" << "    " << "u" << endl;
    for (i = 0; i <= n + 1; i++)
    {
        cout << h_step*i << "    " << u[i] << endl;
    }
    
    delete [] u;
    return 0;
}


//TODO: Change the call of func and or func itself because of unnecessary FLOPs
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
        temp[i] = c/btemp;
        btemp = b - a*temp[i];
        u[i] = (func(h_step, h_step*i) - a*u[i-1])/btemp;
    } 

    for(i = n; i >= 1; i--){
        u[i] -= temp[i+1]*u[i+1];
    }
}

//TODO: Change the call of func and or func itself because of unnecessary FLOPs
void tridiag_poisson (double h_step, unsigned int n, double *u){
    /* Just usable for poisson equation tridiagonal matrices. Reduces number  *
     * of FLOPs to 6N.                                                        */
    unsigned int i;
    vec ftemp(n+1);

    ftemp[1] = func(h_step, h_step*1);
    for (i = 2; i <= n; i++){
        ftemp[i] = ftemp[i-1] + i*func(h_step, h_step*i);
    }

    for (i = n; i >= 1; i--){
        u[i] = (ftemp[i] + u[i+1]*i)/(1.+i);
    }
}

double error (double h_step, unsigned int n, double *u){
    unsigned int i;
    double eps, eps_max, u_ana;
    eps_max = -100.0; // not a nice solution 

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
    double f;
    f = 1.-(1.-exp(-10.))*x-exp(-10.*x);
    return f; 
}

