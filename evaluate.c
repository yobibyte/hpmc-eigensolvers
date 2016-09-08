// TODO read csv with maximum precision
// TODO check for memleaks
// TODO add computing the eigenvectors also
// TODO properly load N_PROBLEMS

#include <cblas.h>
#include <stdio.h>
#include <lapacke.h>
#include <stdlib.h>
#include <string.h>

#define MATRIX_LAYOUT 1
#define MODE 'V'
#define RANGE 'A'
#define N_PROBLEMS 5

struct Eigenproblem {
  int p_size;
  double *eigenvalues;
  // before reduction to symmetric tridiagonal
  double *matrix;
  // diagonal after reduction to symm tridiagonal
  double *D;
  // subdiagonal after reduction to symm tridiagonal
  double *E;
  // orthogonal matrix for solvers
  double *Q;
};

void eye(int dim, double *mat) {
  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      if(i==j) { 
        mat[i*dim+j%dim] = 1.0;
      } else {
        mat[i*dim+j%dim] = 0.0;
      }
    }
  }
}

void print_array(int len, double *array) {
  for(int i=0; i<len;i++){
    printf("%f ", array[i]);
  }
  printf("\n");
}

void load_problems(char *filename, struct Eigenproblem *problems) {
  FILE *f = fopen("data.csv", "r");
  char *curr_line;
  size_t len = 0;
  int n_problems;
  
  getline(&curr_line, &len, f);
  sscanf(curr_line, "%d", &n_problems);

  char *curr_buf;
  char *eigenvalues_str;
  char *curr_matrix_str;

  for(int i=0; i<n_problems;i++) {
    char *end_str = NULL;
    
    getline(&curr_line, &len, f);
    curr_buf        = strtok_r(curr_line, ";", &end_str);
    eigenvalues_str = strtok_r(NULL,      ";", &end_str);
    curr_matrix_str = strtok_r(NULL,      ";", &end_str); 
    
    printf("Reading eigenproblem #%d with dim %s\n", i, curr_buf);
    
    int cpdim;
    sscanf(curr_buf, "%d", &cpdim); 
    int cplen = cpdim*cpdim;
    //int cplen = cpdim*(cpdim+1)/2;
    double curr_eigenvals[cpdim];
    double curr_matrix[cplen];
    
    curr_buf = strtok_r(eigenvalues_str, ",", &end_str);
    sscanf(curr_buf, "%lf", &curr_eigenvals[0]);
    for(int j=1;j<cpdim;j++){
      curr_buf = strtok_r(NULL, ",", &end_str);
      sscanf(curr_buf, "%lf", &curr_eigenvals[j]);
    }    
  
    curr_buf = strtok_r(curr_matrix_str, ",", &end_str);
    sscanf(curr_buf, "%lf", &curr_matrix[0]);
    for(int j=1;j<cplen;j++){
      curr_buf= strtok_r(NULL, ",", &end_str);
      sscanf(curr_buf, "%lf", &curr_matrix[j]);
    }    

    struct Eigenproblem p;
    p.p_size = cpdim;  
    p.matrix = curr_matrix;    
    p.eigenvalues = curr_eigenvals;
    p.D = malloc(sizeof(double) * cpdim);
    p.E = malloc(sizeof(double) * (cpdim-1));
    p.Q = malloc(sizeof(double) * cplen);
    eye(cpdim, p.Q);
    //p.Q = calloc(cplen, sizeof(double));
    double TAU[cpdim-1];
    LAPACKE_dsytrd(LAPACK_ROW_MAJOR, 'U', cpdim, p.matrix, cpdim, p.D, p.E, TAU);

    double V[cpdim];
    for(int j=cpdim-1;j>0;j--){
      for(int k=0;k<cpdim;k++) {
        if (k < j-1){
          V[k] = curr_matrix[k*cpdim+j];
        } else if(k == j-1) {
          V[k] = 1;
        } else {
          V[k] = 0;
        }
      }
      double EYE[cplen];
      eye(j+1, EYE);
      double *Q_curr = calloc(cplen, sizeof(double));
      cblas_dger(CblasRowMajor, cpdim, cpdim, -TAU[j], V, 1, V, 1, Q_curr, cpdim );
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cpdim, cpdim, cpdim, 1.0, p.Q, cpdim, Q_curr, cpdim, 1.0, p.Q, cpdim);
      free(Q_curr);
    }
    
    problems[i] = p;
    
  ;}
  printf("Read %d eigenproblems into the memory.\n", n_problems);

  fclose(f);
  if (curr_line) {
    free(curr_line); 
  }
}


int main(void) {
  struct Eigenproblem problems[N_PROBLEMS];
  load_problems("data.csv", problems);
  
  for (int i=0;i<N_PROBLEMS;i++) {
    struct Eigenproblem p = problems[i];
    double *Z = p.Q;
    // DSTEQR
    lapack_int info = LAPACKE_dsteqr(LAPACK_ROW_MAJOR, MODE, p.p_size, p.D, p.E, Z, p.p_size); 
    /*
    int VL, VU, IL, IU;
    lapack_int M[p.p_size];
    double W[p.p_size];
    double ABSTOL = 0.001;
    int ISUPPZ[p.p_size*2];
    
    // DSTEVX
    int *ifail;
    lapack_int info = LAPACKE_dstevx(LAPACK_ROW_MAJOR, MODE, RANGE, p.p_size, p.D, p.E, VL, VU, IL, IU, ABSTOL, M, W, Z, p.p_size, ifail);
    */
   
    /* 
    //DSTEMR
    lapack_logical *tryrac = malloc(sizeof(lapack_logical));
    lapack_int info = LAPACKE_dstemr(LAPACK_ROW_MAJOR, MODE, RANGE, p.p_size, p.D, p.E, VL, VU, IL, IU, M, W, Z, p.p_size, p.p_size, ISUPPZ, tryrac); 
    
    printf("INFO : %d; ", info); 
    print_array(p.p_size, W);
    
    free(tryrac);
    */
    
    printf("INFO : %d;\n", info); 
  }

  for(int i=0;i<N_PROBLEMS;i++) {
    free(problems[i].D);
  }
  return 0;
}
