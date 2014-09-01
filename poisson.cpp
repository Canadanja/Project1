// calculates the first derivative of a function
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

double func (double x);
void tridiag (double h_step, int n, double * u); 


int main()
{   
    unsigned int n, i;                              // number of steps        //
    double x, x_min, x_max;                         // range of variable      // 
    double *u, u_min, u_max;                        // derichlet boundaries   //
    double h_step;                                  // stepsize               //
    
    n = 50;                                         // number of steps        //
    x_min = 0.;                                     //                        //
    x_max = 1.;
    x = x_min; 
    u_min = 0;                                      // derichlet-boundaries   //
    u_max = 0;
    h_step = (x_max - x_min) / (n + 1);
    
    u = new double[n + 2];
    u[0] = u_min;
    u[n+1] = u_max;
    
    tridiag (h_step, n, u); 
    
    for (i = 0; i <= n + 1; i++)
    {
    cout << u[i] << endl;
    }
    delete [] u;
    return 0;
}

double func (double x)
{
    double f;

    f = 100.*exp(-10.*x);
    return f; 
}

void tridiag (double h_step, int n, double *u){
    double a,b,c;
    double btemp;
    unsigned int i; 
    vec temp(n);

    a = -1;
    c = a ;
    b = 2 ; 

    btemp = b;
    u[i] = func(h_step*1)/btemp;
    for(i = 2 ; i <= n; i++){
        temp[i] = c/btemp;
        btemp = b - a*temp[i];
        u[i] = (func(h_step*i) - a*u[i-1])/btemp;
    } 

    for(i = n-1; i >= 1; i--){
        u[i] -= temp[i+1]*u[i+1];
    }


}

