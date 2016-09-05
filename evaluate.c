// TODO read csv with maximum precision
// TODO check for memleaks
// TODO add computing the eigenvectors also


#include <stdio.h>
//#include <lapack.h>
#include <lapacke.h>
#include <stdlib.h>
#include <string.h>

#define MATRIX_LAYOUT 1
#define MODE 'N'
struct Eigenproblem {
  int p_size;
  double *eigenvalues;
  double *matrix;
};

struct Eigenproblem_ref {
  int p_size;
  double *eigenvalues;
  double *D;
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
  char *cmdarray[512];

  //struct Eigenproblem problems[n_problems];
  
  for(int i=0; i<n_problems;i++) {
    
    read = getline(&curr_line, &len, f);
    char *end_str = NULL;
    curr_buf        = strtok_r(curr_line, ";", &end_str);
    eigenvalues_str = strtok_r(NULL,      ";", &end_str);
    curr_matrix_str = strtok_r(NULL,      ";", &end_str); 
    
    printf("Reading eigenproblem #%d with dim %s\n", i, curr_buf);
    
    int cpsize;
    sscanf(curr_buf, "%d", &cpsize); 
 
    char *eigenvalues_rem;
    curr_eigenval = strtok_r(eigenvalues_str, ",", &eigenvalues_rem);
    double curr_eigenvals[cpsize];
    sscanf(curr_eigenval, "%lf", &curr_eigenvals[0]);
    //printf("%f\t", curr_eigenvals[0]); 
    for(int j=1;j<cpsize;j++){
      curr_eigenval     = strtok_r(NULL, ",", &eigenvalues_rem);
      sscanf(curr_eigenval, "%lf", &curr_eigenvals[j]);
      //printf("%f\t", curr_eigenvals[j]);
    }    
    //printf("\n");
  
    char *curr_elem;
    char *cur_mat_rem;
    curr_elem = strtok_r(curr_matrix_str, ",", &cur_mat_rem);
    double curr_matrix[cpsize*cpsize];
    sscanf(curr_elem, "%lf", &curr_matrix[0]);
    //printf("%lf ", curr_matrix[0]); 
    for(int j=1;j<cpsize*cpsize;j++){
      curr_elem= strtok_r(NULL, ",", &cur_mat_rem);
      sscanf(curr_elem, "%lf", &curr_matrix[j]);
      //printf("%lf ", curr_matrix[j]);
    }    
    //printf("\n");

    
    struct Eigenproblem p;
    p.p_size = cpsize;  
    p.matrix = curr_matrix;    
    p.eigenvalues = curr_eigenvals;  
  
    problems[i] = p;
  }
  printf("Read %d eigenproblems into the memory.\n", n_problems);

  fclose(f);
  if (curr_line) {
    free(curr_line); 
  }
}

//LAPACKE_dsteqr(MATRIX_LAYOUT, MODE, p_size, D, E, Z, ldz);
//LAPACKE_dstevx(MATRIX_LAYOUT, MODE, char range, p_size, D, E, VL, VU, IL, IU, double abstol, M, W, Z, ldz, ifail);
//LAPACKE_dsemr(MATRIX_LAYOUT, MODE, char range, p_size, D, E, VL, VU, IL, IU, M, W, Z, ldz, nzc, isuppz, tryrac);


void solve_dstemr() {
  // MR3 
}

void solve_dstevx(struct Eigenproblem *p, double *result){
  // BX+II  
}

void solve_dsteqr(int n, double *D, double *E, double *result) {
  // QR
  lapack_int ldz = 1;
  double *Z;
  LAPACKE_dsteqr(LAPACK_ROW_MAJOR, MODE, n, D, E, Z, ldz); 
}

int main(void) {
  struct Eigenproblem problems[5]; //TODO CHANGE TO VAR
  load_problems("data.csv", problems);
//
  double D[3] = {1,2,3};
  double E[2] = {1,3};
  double *Z;
  int ldz = 3;
  lapack_int t = LAPACKE_dsteqr(MATRIX_LAYOUT, 'N', 3, D, E, Z, ldz); 
}
