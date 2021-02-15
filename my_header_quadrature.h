#include <iostream>
#include <cmath>
#include <iomanip> //per abilitare le cifre decimali
#include <cstdlib>
#include <fstream>
#include <ctime>

using namespace std;


double MidQuad(double(*Func)(double ),double , double , int );
double TrapezoidQuad(double(*Func)(double ),double , double , int );
double SimpsonQuad(double(*Func)(double ),double , double , int );
double GaussQuad(double(*Func)(double ),double , double ,int ,int);
