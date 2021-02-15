#include "my_header_root.h" //le virgolette vogliono dire di cercare l'header qui nella directory

#define DEBUG FALSE

////////////////////////////////////////////////////////////////////
double Bisection(double (*Func)(double x),double inf, double sup, double tol)
/////////////
//funziona solo per uno zero
/////////////
{
  
  double xm=0.,fm=0.;
  int k=1;
  double fa=Func(inf);
  
  while(1){
    
    xm = 0.5 * (sup + inf);
    fm=Func(xm);
   
    #if DEBUG == TRUE
    cout << "Bisection(): k = " << k << "; [a,b] = [" << inf << "," << sup << "; xm = " << xm << "; dx = " << sup-inf << "; fm = " << fm << endl;
    #endif
    
    k++;
    
    //metodo di bisezione non ridefinisco fa perché quando dovrebbe essere ridefinito il segno è uguale
    if(fm*fa<0){sup=xm;}
    else{inf=xm;}
    
    if(fabs(sup-inf)<tol){
      break;
    }
  }
  
  return xm;
  
}

////////////////////////////////////////////////////////////////////
double FalsePos(double (*Func)(double x),double inf, double sup, double tol)
/////////////
//funziona solo per uno zero
/////////////
{
  double fa=Func(inf), fb=Func(sup);
  double p=0.,q=0, xm=0., fm=0,del=0;
  int k=1;
  
  while(1){
    
    //calcolo parametri retta
    p=(fb-fa)/(sup-inf);//coeff ang.
    q=fa-p*inf;//termine noto
    
    xm=-q/p;
    fm=Func(xm);
    
     #if DEBUG == TRUE
    cout << "FalsePos(): k = " << k << "; [a,b] = [" << inf << "," << sup << "; xm = " << xm << "; |del| = " << del << "; fm = " << fm << endl;
    #endif
    
    k++;
    
    if(fm*fa<0){
      del=fabs(xm-sup);
      sup=xm;
      fb=fm;
    }
    else{
      del=fabs(xm-inf);
      inf=xm;
      fa=fm;
    }
    
    
    if(del<tol){
      break;
    }
    if(k>100){
      cout << "Error!! : iteration > 100 " << endl;
      break;} //controllo: questo metodo può non convergere
    
  }
  
  if(k>100){
    return 0;
  }
  else{
  return xm;
  }
  
}


////////////////////////////////////////////////////////////////////
double Secant(double (*Func)(double x),double inf, double sup, double tol)
/////////////
//funziona solo per uno zero
/////////////
{
  double fa=Func(inf), fb=Func(sup);
  double xm=0., fm=0, dx=sup-inf;
  int k=1;
  
  while(1){
    
    dx=fb*(sup-inf)/(fb-fa);
    
    #if DEBUG == TRUE
    cout << "Secant(): k = " << k << "; xa = " << inf << "; xb = " << sup << "; dx = " << dx << endl;
    #endif
    
    
    //vedi regola per metodo delle secanti
    inf=sup;
    fa=fb;
    sup-=dx;
    fb=Func(sup);
    
    k++;
    
    if(fabs(dx)<tol){
      break;
    }
    if(k>100){
      cout << "Error!! : iteration > 100 " << endl;
      break;} //controllo: questo metodo può non convergere
    
    
  }
  
  if(k>100){
    return 0;
  }
  else{
    return sup;
  }
  
}


////////////////////////////////////////////////////////////////////
double Newton(double (*Func)(double x), double (*Deriv)(double y),double inf, double sup, double tol)
/////////////
//funziona solo per uno zero
/////////////
{
  double  fb=Func(sup), db = Deriv(sup);
  double dx=sup-inf;
  int k=1;
  
  while(1){
    
    dx=fb/db;
    
    #if DEBUG == TRUE
    cout << "Newton(): k = " << k << "; xa = " << inf << "; xb = " << sup << "; dx = " << dx << endl;
    #endif
    
    //vedi regola per metodo di Newton
    sup-=dx;
    fb=Func(sup);
    db=Deriv(sup);
    
    k++;
    
    if(fabs(dx)<tol){
      break;
    }
    if(k>1000){
      cout << "Error!! : iteration > 1000 " << endl;
      break;} //controllo: questo metodo può non convergere
    
    
  }
  
  if(k>1000){
    return 0;
  }
  else{
    return sup;
  }
  
}
/////////////////////////////////////////////////
void Braketing( double (*Func)(double x), double xa, double xb, int n, double *xL, double *xR, int& Nr)
//
//
{
  double dx = fabs(xb-xa)/n, x=xa;
  int i=0;
  int j=0;
  
  while(1){
    
    xL[i]=x;
    xR[i]=x+dx;
    
    if( Func(x) * Func(x+dx) < 0){
      i++;
    }
    
    //passo all'intervallo dopo
    x+=dx;
    j++;
    
    if( j>=n ) {break;}
    
  }
  
  Nr=i;
}
