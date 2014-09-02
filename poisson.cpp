// calculates the one-dimensional poisson equation
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

double func (double, double);
void tridiag (double, int, double *); 


int main()
{   
    unsigned int n, i;                              
    double x, x_min, x_max;                        
    double *u, u_min, u_max;                      
    double h_step;                               
    
    n = 50;                                         // number of steps        //
    x_min = 0.;                                     // area of calculation    //
    x_max = 1.;                                     //                        //
    x = x_min;                                      // initial value for var. // 
    u_min = 0;                                      // dirichlet-boundaries   //
    u_max = 0;                                      //                        //
    h_step = (x_max - x_min) / (n + 1);             // stepsize               // 
    
    u = new double[n];                              // initialize array for   //
    u[0] = u_min;                                   // solution u             //
    u[n] = u_max;                                   //                        // 
    
    // tridiag solves the equation Au = f
    tridiag (h_step, n, u); 
    
    // output
    for (i = 0; i <= n + 1; i++)
    {
        cout << h_step*i << "    " << u[i] << endl;
    }
    
    delete [] u;
    return 0;
}

double func (double h_step, double x)
{
    // second derivative rhs-function 
    double f;

    f = h_step*h_step*100.*exp(-10.*x);
    return f; 
}

void tridiag (double h_step, int n, double *u){
    /* Uses the thomas algorithm. Code from "Numerical Recipes,               *  
     * Third Edition", p. 56 f.                                               */
    double a,b,c;
    double btemp;
    unsigned int i; 
    vec temp(n);

    a = -1;
    c = a ;
    b = 2 ; 

    btemp = b;
    u[1] = func(h_step, h_step*1)/btemp;
    for(i = 2 ; i <= n; i++){
        temp[i] = c/btemp;
        btemp = b - a*temp[i];
        u[i] = (func(h_step, h_step*i) - a*u[i-1])/btemp;
    } 

    for(i = n-1; i >= 1; i--){
        u[i] -= temp[i+1]*u[i+1];
    }


}

