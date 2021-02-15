#include "my_header_root.h" //le virgolette vogliono dire di cercare l'header qui nella directory

//////////////////////////////////////
//
//
//////////////////////////////////////
void Matrix_times_Vector(double **M, double *v, double *S, int n){
  
  for(int i=0; i<n; i++){
    
    S[i] = 0;
    
    for(int j=0; j<n; j++){
      
      S[i] += M[i][j]*v[j];
    }
    
  }
  
  
}

//////////////////////////////////////
//
//
//////////////////////////////////////
void Matrix_product(double **M1, double **M2, double **S, int n){
  
  for(int i=0; i<n; i++){
    
    for(int j=0; j<n; j++){
      
      S[i][j] = 0;
      
      for(int k=0; k<n; k++){
        S[i][j] += M1[i][k]*M2[k][j];
      }
    }
  }
}

//////////////////////////////////////
//
//
//////////////////////////////////////
void Print_Matrix(double **M, int n, int m){
  
  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      
      cout << setw(10) << right << M[i][j] << " ";
    }
    cout << endl;
  }
  
}
//////////////////////////////////////
//
//
//////////////////////////////////////
void Print_Vector(double *v, int n){
  
  for(int i=0; i<n; i++){
    cout << setw(10) << right << v[i] << " ";
  }
  cout << endl;
  
}
//////////////////////////////////////
//
//
//////////////////////////////////////
void Swap_Rows(double **A, double *b, int i1, int i2, int n){
  
  double temp;
  for(int j=0; j<n; j++){
    temp = A[i1][j];
    A[i1][j] = A[i2][j];
    A[i2][j] = temp;
  }
  
  temp = b[i1];
  b[i1] = b[i2];
  b[i2] = temp;
}

void Gauss_elim(double **A, double *b, int n){
  
  double m;
  int jmax;
  double tmp;
  
  for(int k=0; k<n-1; k++){
    
    //partial pivoting
    tmp = A[k][k];
    jmax = k;
    for(int j=k+1; j<n; j++){
      if(fabs(A[j][k])>fabs(tmp)){
        tmp = A[j][k];
        jmax = j;
      }
    }
    
    Swap_Rows(A,b,k,jmax,n);
    
    for(int i=k+1; i<n; i++){
      m = A[i][k]/A[k][k];
      for(int j=k+1; j<n; j++){
        A[i][j] -= m*A[k][j];
      }
      A[i][k] = 0.0;
      b[i] -= m*b[k];
      
    }
    
  }
  
}
//////////////////////////////////////
//
//
//////////////////////////////////////
void Tridiag_Solver(double *a, double *b, double *c, double *r, double *x, int n){
  
  double *h = new double[n];
  double *p = new double[n];
  double den;
  
  h[0] = c[0]/b[0], p[0] = r[0]/b[0];
  
  for(int i=1; i<n; i++){
    den = b[i]-a[i]*h[i-1];
    h[i] = c[i]/den;
    p[i] = (r[i]-a[i]*p[i-1])/den;
  }
  
  x[n-1] = p[n-1];
  
  for(int i=n-2; i>=0; i--){
    x[i] = p[i] - h[i]*x[i+1];
  }
  
  delete[] p;
  delete[] h;
}
