// TODO read csv with maximum precision
// TODO check for memleaks
// TODO add computing the eigenvectors also
// TODO properly load N_PROBLEMS

#include <stdio.h>
//#include <lapack.h>
#include <lapacke.h>
#include <stdlib.h>
#include <string.h>

#define MATRIX_LAYOUT 1
#define MODE 'N'
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
};

void load_problems(char *filename, struct Eigenproblem *problems) {
  FILE *f = fopen("data.csv", "r");
  char *curr_line;
  size_t len = 0;
  size_t read;

  int n_problems;
  
  read = getline(&curr_line, &len, f);
  sscanf(curr_line, "%d", &n_problems);

  char *curr_buf;
  char *eigenvalues_str;
  char *curr_matrix_str;
  char *curr_eigenval;
  char **cmdarray;

  for(int i=0; i<n_problems;i++) {
    read = getline(&curr_line, &len, f);
    char *end_str = NULL;
    curr_buf        = strtok_r(curr_line, ";", &end_str);
    eigenvalues_str = strtok_r(NULL,      ";", &end_str);
    curr_matrix_str = strtok_r(NULL,      ";", &end_str); 
    
    printf("Reading eigenproblem #%d with dim %s\n", i, curr_buf);
    
    int cpdim;
    sscanf(curr_buf, "%d", &cpdim); 
 
    char *eigenvalues_rem;
    curr_eigenval = strtok_r(eigenvalues_str, ",", &eigenvalues_rem);
    double curr_eigenvals[cpdim];
    sscanf(curr_eigenval, "%lf", &curr_eigenvals[0]);
    //printf("%f\t", curr_eigenvals[0]); 
    for(int j=1;j<cpdim;j++){
      curr_eigenval     = strtok_r(NULL, ",", &eigenvalues_rem);
      sscanf(curr_eigenval, "%lf", &curr_eigenvals[j]);
      //printf("%f\t", curr_eigenvals[j]);
    }    
    //printf("\n");
  
    char *curr_elem;
    char *cur_mat_rem;
    curr_elem = strtok_r(curr_matrix_str, ",", &cur_mat_rem);
    int cplen = (int) (cpdim*(cpdim+1)/2);
    
    double curr_matrix[cplen];
    sscanf(curr_elem, "%lf", &curr_matrix[0]);
    //printf("%lf ", curr_matrix[0]); 
    for(int j=1;j<cplen;j++){
      curr_elem= strtok_r(NULL, ",", &cur_mat_rem);
      sscanf(curr_elem, "%lf", &curr_matrix[j]);
      //printf("%lf ", curr_matrix[j]);
    }    
    //printf("\n");

        
    struct Eigenproblem p;
    p.p_size = cpdim;  
    p.matrix = curr_matrix;    
    p.eigenvalues = curr_eigenvals;


    p.D = malloc(sizeof(double) * cpdim);
    p.E = malloc(sizeof(double) * (cpdim-1));
    //double D[cpdim];
    //double E[cpdim-1];
    double TAU[cpdim-1];
    //lapack_int LAPACKE_dsytrd(int matrix_layout, char uplo, lapack_int n, double* a,
    //lapack_int lda, double* d, double* e, double* tau )
    
    LAPACKE_dsytrd(LAPACK_ROW_MAJOR, 'U', cpdim, p.matrix, cpdim, p.D, p.E, TAU);
    /*
    for (int l=0;l<cplen;l++) {
      printf("%f   ", p.D[l]);
    }  
    printf("\n");
    */
    problems[i] = p;
  }
  printf("Read %d eigenproblems into the memory.\n", n_problems);

  fclose(f);
  if (curr_line) {
    free(curr_line); 
  }
}

void print_array(int len, double *array) {
  for(int i=0; i<len;i++){
    printf("%f ", array[i]);
  }
  printf("\n");
}

int main(void) {
  struct Eigenproblem problems[N_PROBLEMS];
  load_problems("data.csv", problems);
  
  double *Z;
  
  for (int i=0;i<N_PROBLEMS;i++) {
    struct Eigenproblem p = problems[i];
   
    /*
    // DSTEQR
    //lapack_int info = LAPACKE_dsteqr(LAPACK_ROW_MAJOR, MODE, p.p_size, p.D, p.E, NULL, p.p_size); 
    */

    /*
    // DSTEVX
    lapack_int *M = malloc(sizeof(double)*p.p_size);
    int VL, VU, IL, IU;
    double *W = malloc(sizeof(double)*p.p_size);
    double ABSTOL = 0.001;
    int *ifail;
    lapack_int info = LAPACKE_dstevx(LAPACK_ROW_MAJOR, MODE, RANGE, p.p_size, p.D, p.E, VL, VU, IL, IU, ABSTOL, M, W, Z, p.p_size, ifail);
    */
    
    //DSTEMR
    lapack_int *M = malloc(sizeof(double)*p.p_size);
    int VL, VU, IL, IU;
    double *W = malloc(sizeof(double)*p.p_size);
    double ABSTOL = 0.001;
    int ISUPPZ[p.p_size*2];
    lapack_logical *tryrac = malloc(sizeof(lapack_logical));
    lapack_int info = LAPACKE_dstemr(LAPACK_ROW_MAJOR, MODE, RANGE, p.p_size, p.D, p.E, VL, VU, IL, IU, M, W, Z, p.p_size, p.p_size, ISUPPZ, tryrac); 
    
    printf("INFO : %d; ", info); 
    print_array(p.p_size, W);
  }
}
