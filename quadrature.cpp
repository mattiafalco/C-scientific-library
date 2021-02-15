#include "my_header_quadrature.h"



///////////////////////////////////////////////////////////////////
double MidQuad(double(*Func)(double x),double xa, double xb, int N)
{
  double integral=0.0;
  double h = (xb-xa)/(double)N; //base dell'intevallo
  double xi=xa+h/2;
  
  for(int i=0;i<N;i++)
  {
    xi=xa+h/2+i*h;
    integral+=Func(xi);
  }
  return integral*h;
}
///////////////////////////////////////////////////////////////////
double TrapezoidQuad(double(*Func)(double x),double xa, double xb, int N)
{
  double integral=(Func(xa)+Func(xb))/2;
  double h = (xb-xa)/(double)N; //base dell'intevallo
  double xi=xa;
  
  for(int i=1;i<N;i++)
  {
    xi=xa+i*h;
    integral+=Func(xi);
  }
  return integral*h;
}
///////////////////////////////////////////////////////////////////
double SimpsonQuad(double(*Func)(double x),double xa, double xb, int N)
{
  double integral=Func(xa)+Func(xb);
  double h = (xb-xa)/(double)N; //base dell'intevallo
  double xi=xa;
  int w=4;
  
  for(int i=1;i<N;i++)
  {
    xi=xa+i*h;
    integral+=Func(xi)*w;
    w=6-w;
  }
  return integral*h/3;
}

///////////////////////////////////////////////////////////////////
double GaussQuad(double (*Func)(double x),double xa, double xb, int N, int ngauss)
//////////////////
//////////////////
{
  
  double p[16];
  double w[16];
  double x[16];
  
  //punti gaussiani
  if(ngauss==2)
  {
    p[0]=sqrt(1./3.);   w[0]=1.;
    p[1]=-sqrt(1./3.);  w[1]=1.;
    
  }
  else if(ngauss==3)
  {
    p[0]=0;              w[0]=8./9;
    p[1]=sqrt(3./5);     w[1]=5./9;
    p[2]=-sqrt(3./5);    w[2]=5./9;
    
  }
  else if(ngauss==4)
  {
    p[0]=0.339981043584856;    w[0]=0.888888888888889;
    p[1]=-p[0];                w[1]=w[0];
    p[2]=0.861136311594053;    w[2]=0.555555555555556;
    p[3]=-p[4];                w[3]=w[2];
  
  }
  else if(ngauss==5)
  {
    p[0]=0.000000000000000;   w[0]=0.568888888888889;
    p[1]=0.538469310105683;   w[1]=0.478628670499366;
    p[2]=-p[1];               w[2]=w[1];
    p[3]=0.906179845938664;   w[3]=0.236926885056189;
    p[4]=-p[3];               w[4]=w[3];
    
  }
  else
  {
    cout << "GausQuad() not defined for ngauss = " << ngauss << endl;
    exit(1);
  }
    
  double sum=0.0, sumk=0.0;
  double h = (xb-xa)/(double)N; //base dell'intevallo
  
  for(int i=0; i<N; i++)
  {
    sumk=0.0;  //somma parziale = 0
    
    for(int j=0; j<ngauss; j++){x[j]=0.5*h*p[j]+xa+i*h+0.5*h;}
    
    for(int j=0; j<ngauss; j++)
    {
      sumk+=w[j]*Func(x[j]);
    }
    
    sum+=sumk;
  }
  return 0.5*h*sum;
}

