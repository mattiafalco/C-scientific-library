#include "my_prototype.h" //le virgolette vogliono dire di cercare l'header qui
// nella directory
//////////////////////////////
//
//Project : three body problem
//
//Sun->1
//Earth->2
//Jupiter->3
//
//ultima modifica 11/01/20

#define SESSION 2

#define SUN TRUE
#define GIRI FALSE
#define MTIMES 1000

#define NP 3
#define MH 1   //Earth
#define MS 332946 //Sun
#define MJ 317*MTIMES  //Jupiter
#define TIME 2.8e9  //time


void acceleration(double *r,double *a);
void acc_analitic(double *r,double *a);
void Solution(double t, double *X);
double Conversion(double t);


int main() {

  #if SESSION==1

  string fname1 = "analitic.dat";
  ofstream fdata1;    // declare Output stream class to operate on files
  fdata1.open("analitic.dat"); // open output file
  fdata1 << setiosflags(ios::scientific)<<setprecision(14);

  string fname2 = "error.dat";
  ofstream fdata2;    // declare Output stream class to operate on files
  fdata2.open("error.dat"); // open output file
  fdata2 << setiosflags(ios::scientific)<<setprecision(14);

  int giri_max = 40;

  //analitic
  double X_an[2];

  //initial condition position
  double x1 = 1., y1 = 0;
  double r = sqrt(x1*x1+y1*y1);

  double a = 1.;
  //initial condition velocity
  double vx1 = 0., vy1 = sqrt(a/r);


  //array position and velocity 2nd order
  double *x = new double[2];


  double h = 0.1;
  double temp=vx1;
  int giri=0;
  double t=0.;
  double err = 0.;
  double err4 = 0.;
  double E = 0.;

  //ciclo per dimezzare h
  for(int i=0; i<10; i++){
  //filling array
  x[0] = x1; x[1] = y1;
  v[0] = vx1; v[1] = vy1;

  //initial condition
  temp=vx1;
  giri=0;
  t=0.;

    while(1){

      //intergration with Position_Verlet
      Position_Verlet(x,v,2,h,acc_analitic);
      t+=h;

      //analitic solution at time t
      Solution(t,X_an);

      r = sqrt(x[0]*x[0] + x[1]*x[1]);

      E = 0.5*(v[0]*v[0]+v[1]*v[1]) - 1/r;

      //check Earth
      if(temp*v[0] < 0){giri++;}
      temp = v[0];

      //save data
      fdata1 << t << " ";
      fdata1 << x[0] << " " << x[1] << " " << X_an[0] << " " << X_an[1] << " "
      << err << " " << E << " " ;
      fdata1 << endl;

      if(giri > 2*giri_max){break;}

    }

    //error
    err = sqrt((x[0]-X_an[0])*(x[0]-X_an[0]) + (x[1]-X_an[1])*(x[1]-X_an[1]));

    cout << " numero giri Terra : " << (double)giri/2.0 << endl;

    fdata2 << h << " " << err << " " << endl;
    fdata1 << endl <<endl;
    h/=2;
  }

  fdata1.close();
  fdata2.close();
  #endif

  #if SESSION==2

  cout << setiosflags(ios::scientific)<<setprecision(14); //abilita tutte
  //le cifre decimali

  string fname1 = "sun.dat";
  ofstream fdata1;    // declare Output stream class to operate on files
  fdata1.open("sun.dat"); // open output file
  fdata1 << setiosflags(ios::scientific)<<setprecision(14);

  string fname2 = "earth.dat";
  ofstream fdata2;    // declare Output stream class to operate on files
  fdata2.open("earth.dat"); // open output file
  fdata2 << setiosflags(ios::scientific)<<setprecision(14);

  string fname3 = "jupiter.dat";
  ofstream fdata3;    // declare Output stream class to operate on files
  fdata3.open("jupiter.dat"); // open output file
  fdata3 << setiosflags(ios::scientific)<<setprecision(14);

  int giri_max = 40;
  double anni = 20;

  //initial condition position
  double x1 = 0, y1 = 0;
  double x2 = 1, y2 = 0;
  double x3 = 5.2, y3 = 0;

  //parametri a usare per impostare orbite generiche
  double r12 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
  double r13 = sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3));
  double r23 = sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3));

  double a1 = 1., a2 = 1., a3 = 1.;
  //initial condition velocity
  double vx1 = 0., vy1= 0;
  double vx2 = 0., vy2 = 557.779; //sqrt(a1*MS/r12 + a2*MJ/r23);
  double vx3 = 0., vy3 = 244.465; //sqrt(a1*MS/r13 + a3*MH/r23);

  //array position and velocity
  double *x = new double[2*NP];
  double *v = new double[2*NP];


  //filling array
  x[0] = x1; x[1] = y1;
  x[2] = x2; x[3] = y2;
  x[4] = x3; x[5] = y3;

  v[0] = vx1; v[1] = vy1;
  v[2] = vx2; v[3] = vy2;
  v[4] = vx3; v[5] = vy3;

  double h = 0.00001;
  double temp=vx2;
  int giri=0;
  double t=0.;

  while(1){

    //intrgration with Position_Verlet
    Position_Verlet(x,v,2*NP,h,acceleration);
    t+=h;

    //check Earth
    if(temp*v[2] < 0){giri++;}
    temp = v[2];

    //save data
    fdata1 << t << " ";
    fdata1 << x[0] << " " << x[1] << " ";
    fdata1 << endl;

    fdata2 << t << " ";
    fdata2 << x[2] << " " << x[3] << " " << giri << " " ;
    fdata2 << endl;

    fdata3 << t << " ";
    fdata3 << x[4] << " " << x[5] << " ";
    fdata3 << endl;

    #if GIRI==TRUE
    if(giri > 2*giri_max){break;}
    #endif
    #if GIRI==FALSE
    if(Conversion(t) > anni ){break;}
    #endif

  }

  cout << " numero giri Terra : " << (double)giri/2.0 << endl;
  cout << " tempo trascorso : " << Conversion(t) << " anni " << endl;

  fdata1.close();
  fdata2.close();
  fdata3.close();

  #endif

  return 0;
}

/////////////////////////////////////////////////////
void Solution(double t, double *X)
//
//
//this function gives the analitic solution of
//the equation for the two body problem
//
{

    X[0] = cos(t);
    X[1] = sin(t);

}
/////////////////////////////////////////////////////
void acc_analitic(double *r,double *a)
//
//acceleration function
//r = {x,y}
//
{
  double x1,y1;

  x1 = r[0]; y1 = r[1];

  double r1 = sqrt(x1*x1+y1*y1);

  //x
  a[0] = -x1/(r1*r1*r1) ;
  //y
  a[1] = -y1/(r1*r1*r1) ;


}


/////////////////////////////////////////////////////
void acceleration(double *r,double *a)
//
//acceleration function
//r = {x1,y1,x2,y2,x3,y3}
//
{
  double x1,x2,x3,y1,y2,y3;

  x1 = r[0]; y1 = r[1];
  x2 = r[2]; y2 = r[3];
  x3 = r[4]; y3 = r[5];

  double r12 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
  double r13 = sqrt((x1-x3)*(x1-x3) + (y1-y3)*(y1-y3));
  double r23 = sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3));

  #if SUN==FALSE
  /////////
  //Sun
  ////////
  //x1
  a[0] = 0.;
  //y1
  a[1] = 0.;
  #endif


  #if SUN==TRUE
  /////////
  //Sun
  ////////
  //x1
  a[0] = -MH*(x1-x2)/(r12*r12*r12) - MJ*(x1-x3)/(r13*r13*r13);
  //y1
  a[1] = -MH*(y1-y2)/(r12*r12*r12) - MJ*(y1-y3)/(r13*r13*r13);
  #endif

  /////////
  //Earth
  ////////
  //x2
  a[2] = -MS*(x2-x1)/(r12*r12*r12) - MJ*(x2-x3)/(r23*r23*r23);
  //y2
  a[3] = -MS*(y2-y1)/(r12*r12*r12) - MJ*(y2-y3)/(r23*r23*r23);

  /////////
  //Jupiter
  ////////
  //x3
  a[4] = -MS*(x3-x1)/(r13*r13*r13) - MH*(x3-x2)/(r23*r23*r23);
  //y3
  a[5] = -MS*(y3-y1)/(r13*r13*r13) - MH*(y3-y2)/(r23*r23*r23);



}

////////////////////////////////////////////////////////////////
double Conversion(double t){
/////////
//
// converte il tempo adimensionale in anni
//
////////

  double yr = t*TIME/(60.0*60.0*24.0*365.0);

  return yr;
}
