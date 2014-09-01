// calculates the first derivative of a function
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <armadillo>

using namespace std;
using namespace arma;

//void derivative1 (int, double, double, double *, double *);
//void derivative2 (int, double, double, double *, double *);
//void output(double *, double *, double, int)
double func (double);
void sec_derivative (


int main()
{   
    unsigned int n;                                 // number of steps        //
    double x, x_min, x_max;                         // range of variable      // 
    double *u, u_min, u_max;                        // derichlet boundaries   //
    double h_step;                                  // stepsize               //
    
    n = 50;
    x_min = 0.;
    x_max = 1.;
    x = x_min; 
    u_min = 0;                                      // derichlet-boundaries   //
    u_max = 0;
    h_step = (x_max - x_min) / (n + 1);
    
    u = new double[n + 2];
    u[0] = u_min;
    u[n+1] = u_max;
    
    sec_derivative (h_step, x, 
    // Dimension
    //int n = atoi(argv[1]);
    //mat A(n,p);
    //C = zeros<mat> (n,m);
    //int _numRowA = A.n_rows;
    
    
    cout << tridiag[0] << func(20.) << endl;
    delete [] u;
    return 0;
}

//void derivative1 (int number_of_steps, double x, double initial_step, double \
//        *h_step, double *computed_derivative1)
//{
//}

double func (double x)
{
    double f;

    f = 100.*exp(-10.*x);
    return f; 
}

void sec_derivative (int number_of_steps, double x, double initial_step, \
        double *h_step, double *computed_derivative)
{
    vec tridiag(2);
    tridiag[0] = -1.;
    tridiag[1] = 2;

     
}

void tridiag (int n, double u, double initial_step, \
        double *h_step, double *computed_derivative){
    double btemp;
    vec temp(number_of_steps);

    btemp = b[1];
    u[i] = f[1]/btemp;
    for(i = 2 ; i <= n; i++){
        temp[i] = c[i-1]/temp;
        btemp = b[i] - a[i]*temp[i];
        u[i] = (f[i] - a[i]*u[i-1])/temp;
    } 

    for(i = n-1; i >= 1; i--){
        u[i] -= temp[i+1]*u[i+1];
    }


}

/*void derivative2 (int number_of_steps, double x, double initial_step, double \
        *h_step, double *computed_derivative); */

