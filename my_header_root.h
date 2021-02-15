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

//prototipi di funzioni
double Bisection(double (*Func)(double ),double , double , double );
double FalsePos(double (*Func)(double ),double , double , double );
double Secant(double (*Func)(double ),double , double , double );
double Newton(double (*Func)(double ), double (*Deriv)(double ),double , double , double );
void Braketing( double (*Func)(double ), double , double , int , double *, double *, int& );
