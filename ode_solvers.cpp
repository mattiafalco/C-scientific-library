#include "my_prototype.h" //le virgolette vogliono dire di cercare l'header qui nella directory

////////////////////////////////////////////////////////
void EulerStep(double t, double *Y, void (*Ydot)(double , double *, double *), double dt, int neq)
//
//la funzione integra le condizioni iniziali Y(t)  dato il passo dt, le funzioni derivate,
//il numero di equazioni e le condizioni iniziali Y
//
{

  int n;
  double rhs[neq];

  //La funzione Ydot valuta assegna a rhs i valori del right hand side a partire da t e Y
  Ydot(t,Y,rhs);
  for (n = 0; n<neq ; n++){
    Y[n] += dt*rhs[n];
  }
}
 //////////////////////////////////////////////////////
void RKStepMid(double t, double *Y, void (*Ydot)(double , double *, double *), double dt, int neq)
//
//la funzione integra le condizioni iniziali Y(t)  dato il passo dt, le funzioni derivate,
//il numero di equazioni e le condizioni iniziali Y
//
{
  double Ys[neq],k1[neq],k2[neq];

  Ydot(t,Y,k1);

  for(int i=0; i<neq; i++){
    Ys[i] = Y[i] +0.5*dt*k1[i];
  }

  Ydot(t+0.5*dt,Ys,k2);

  for(int i=0; i<neq; i++){
    Y[i] += dt*k2[i];
  }


}
//////////////////////////////////////////////////////
void RKStepTrap(double t, double *Y, void (*Ydot)(double , double *, double *), double dt, int neq)
//
//la funzione integra le condizioni iniziali Y(t)  dato il passo dt, le funzioni derivate,
//il numero di equazioni e le condizioni iniziali Y
//
{
  double Ys[neq],k1[neq],k2[neq];

  Ydot(t,Y,k1);

  for(int i=0; i<neq; i++){
    Ys[i] = Y[i] + dt*k1[i];
  }

  Ydot(t+dt,Ys,k2);

  for(int i=0; i<neq; i++){
    Ys[i] += 0.5*dt*(k2[i]+k1[i]);
  }


}

//////////////////////////////////////////////////////
void RKStep4th(double t, double *Y, void (*Ydot)(double , double *, double *), double dt, int neq)
//
//la funzione integra le condizioni iniziali Y(t)  dato il passo dt, le funzioni derivate,
//il numero di equazioni e le condizioni iniziali Y
//
{
  double Ys1[neq],Ys2[neq],Ys3[neq],k1[neq],k2[neq],k3[neq],k4[neq];
  //iniz k1
  Ydot(t,Y,k1);

  for(int i=0; i<neq; i++){
    Ys1[i] = Y[i] + 0.5*dt*k1[i];
  }
  //iniz k2
  Ydot(t+0.5*dt,Ys1,k2);

  for(int i=0; i<neq; i++){
    Ys2[i] = Y[i] + 0.5*dt*k2[i];
  }
  //iniz k3
  Ydot(t+0.5*dt,Ys2,k3);

  for(int i=0; i<neq; i++){
    Ys3[i] = Y[i] + dt*k3[i];
  }
  //iniz k4
  Ydot(t+dt,Ys3,k4);
  for(int i=0; i<neq; i++){
    Y[i] += dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
  }
}

///////////////////////////////////////////////////////////
void Position_Verlet(double *x, double *v, int np, double dt, void (*acceleration)(double *r,double *a))
//
//metodo simpletticco, integra l'equazione dati gli array di posizioe e velocità
//la funzione acceleration prende i valori di x e riempie a[] con le accelerazioni
//
{
  double x1[np],a[np];

  //calcolo x1 dopo mezzo step
  for(int i=0; i<np; i++){
    x1[i] = x[i] + 0.5*dt*v[i];
  }

  //assegno ad a i valori di accelerazione di x1
  acceleration(x1,a);

  //calcolo v dopo 1 passo
  for(int i=0; i<np; i++){
    v[i] += dt*a[i];
  }

  //calcolo x dopo 1 passo
  for(int i=0; i<np; i++){
    x[i] = x1[i] + 0.5*dt*v[i];
  }
}

/////////////////////////////////////////////
void Symplectic4(double *x, double *v, int np, double dt, void (*acceleration)(double *r,double *a))
//
//metodo simpletticco, integra l'equazione dati gli array di posizioe e velocità
//la funzione acceleration prende i valori di x e riempie a[] con le accelerazioni
//vedi Wikipedia : https://en.wikipedia.org/wiki/Symplectic_integrator
{
  double gamma = 1/(2-pow(2,1/3));

  double a[np],v1[np],x1[np],v2[np],x2[np],v3[np],x3[np];

  //calcolo x1
  for(int i=0; i<np; i++){
    x1[i] = x[i] + 0.5*gamma*dt*v[i];
  }

  //calcolo a(x1)
  acceleration(x1,a);

  //calcolo v1
  for(int i=0; i<np; i++){
    v1[i] = v[i] + gamma*dt*a[i];
  }

  //calcolo x2
  for(int i=0; i<np; i++){
    x2[i] = x1[i] + (1-gamma)*0.5*dt*v1[i];
  }

  //calcolo a(x2)
  acceleration(x2,a);

  //calcolo v2
  for(int i=0; i<np; i++){
    v2[i] = v1[i] + (1-2*gamma)*dt*a[i];
  }

  //calcolo x3
  for(int i=0; i<np; i++){
    x3[i] = x2[i] + (1-gamma)*dt*v2[i]*0.5;
  }

  //calcolo a(x3)
  acceleration(x3,a);

  //calcolo v3
  for(int i=0; i<np; i++){
    v3[i] = v2[i] + (gamma)*dt*a[i];
  }

  //calcolo x
  for(int i=0; i<np; i++){
    x[i] = x3[i] + gamma*dt*v3[i]*0.5;
  }

  //calcolo v;
  for(int i=0; i<np; i++){
    v[i] = v3[i] ;
  }


}
