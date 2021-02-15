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

//prototipi funzioni
void EulerStep(double t, double *Y, void (*Ydot)(double , double *, double *), double dt, int neq);
void RKStepMid(double t, double *Y, void (*Ydot)(double , double *, double *), double dt, int neq);
void RKStepTrap(double t, double *Y, void (*Ydot)(double , double *, double *), double dt, int neq);
void RKStep4th(double t, double *Y, void (*Ydot)(double , double *, double *), double dt, int neq);
void Position_Verlet(double *x, double *v, int np, double dt, void (*acceleration)(double *r,double *a));
