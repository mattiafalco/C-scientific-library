//header utili
#include <iostream>
#include <cmath>
#include <iomanip> //per abilitare le cifre decimali
#include <cstdlib>
#include <fstream>
#include <ctime>

#define TRUE 1
#define FALSE 0

using namespace std;

//prototipi integratori
double MidQuad(double(*Func)(double ),double , double , int );
double TrapezoidQuad(double(*Func)(double ),double , double , int );
double SimpsonQuad(double(*Func)(double ),double , double , int );
double GaussQuad(double(*Func)(double ),double , double ,int ,int);


//prototipi funzioni root_finder
double Bisection(double (*Func)(double ),double , double , double );
double FalsePos(double (*Func)(double ),double , double , double );
double Secant(double (*Func)(double ),double , double , double );
double Newton(double (*Func)(double ), double (*Deriv)(double ),double , double , double );
void Braketing( double (*Func)(double ), double , double , int , double *, double *, int& );


//prototipi funzioni ode_solver
void EulerStep(double t, double *Y, void (*Ydot)(double , double *, double *), double dt, int neq);
void RKStepMid(double t, double *Y, void (*Ydot)(double , double *, double *), double dt, int neq);
void RKStepTrap(double t, double *Y, void (*Ydot)(double , double *, double *), double dt, int neq);
void RKStep4th(double t, double *Y, void (*Ydot)(double , double *, double *), double dt, int neq);
void Position_Verlet(double *x, double *v, int np, double dt, void (*acceleration)(double *r,double *a));
void Symplectic4(double *x, double *v, int np, double dt, void (*acceleration)(double *r,double *a));


//prototipi matrix_reduction
void Matrix_times_Vector(double **M, double *v, double *S, int n);
void Matrix_product(double **M1, double **M2, double **S, int n);
void Print_Matrix(double **M, int n, int m);
void Swap_Rows(double **A, double *b, int i1, int i2, int n);
void Gauss_elim(double **A, double *b, int n);
void Print_Vector(double *v, int n);
void Tridiag_Solver(double *a, double *b, double *c, double *r, double *x, int n);
