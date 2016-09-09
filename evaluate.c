// TODO read csv with maximum precision
// TODO check for memleaks
// TODO add computing the eigenvectors also
// TODO properly load N_PROBLEMS

#include <cblas.h>
#include <stdio.h>
#include <lapacke.h>
#include <stdlib.h>
#include <string.h>

#define DEBUG 0

#define MATRIX_LAYOUT 1
#define MODE 'V'
#define RANGE 'A'
#define N_PROBLEMS 1000
#define ABSTOL 0.000001

#define get_ticks(var) {\
      unsigned int __a, __d;\
      asm volatile("rdtsc" : "=a" (__a), "=d" (__d));\
      var = ((unsigned long) __a) | (((unsigned long) __d) << 32); \
   } while(0)

struct Eigenproblem {
  int p_size;
  double *eigenvalues;
  // before reduction to symmetric tridiagonal
  // double *matrix;
  // diagonal after reduction to symm tridiagonal
  double *D;
  // subdiagonal after reduction to symm tridiagonal
  double *E;
  // orthogonal matrix for solvers
  double *Q;
};

void construct_eigenproblem(struct Eigenproblem *p, int size, double *eigenvalues, double *D, double *E, double *Q) {
  p->eigenvalues = malloc(sizeof(double)*size);
  p->D = malloc(sizeof(double)*size);
  p->E = malloc(sizeof(double)*(size-1));
  p->Q = malloc(sizeof(double)*size*size);

  p->p_size = size;
  memcpy(p->eigenvalues, eigenvalues, sizeof(double)*size);
  memcpy(p->D, D, sizeof(double)*size);
  memcpy(p->E, E, sizeof(double)*(size-1));
  memcpy(p->Q, Q, sizeof(double)*size*size);
}

void destroy_eigenproblem(struct Eigenproblem *p) {
  free(p->eigenvalues);
  free(p->D);
  free(p->E);
  free(p->Q);
}

double get_relative_accuracy(double *real_values, double *computed_values, int len) {
  // Returns mean relative accuracy across all eigenvalues
  // notes on relative accuracy for eigensolvers:
  // http://gauss.uc3m.es/web/personal_web/molera/talks/cedya05charla.pdf
  double res = 0;
  for(int i=0; i<len;res+=abs(computed_values[i]-real_values[i])/real_values[i], i++);
  return 1.0 - res/((double) len);
} 

double get_absolute_accuracy(double *real_values, double *computed_values, int len) {
  // Returns mean absolute accuracy across all eigenvalues
  double res = 0;
  for(int i=0; i<len;res+=abs(computed_values[i]-real_values[i]), i++) {
  }
  return 1.0 - res/((double) len);
}
  
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
    
    if(DEBUG) {
      printf("Reading eigenproblem #%d with dim %s\n", i, curr_buf);
    }
    
    int cpdim;
    sscanf(curr_buf, "%d", &cpdim); 
    int cplen = cpdim*cpdim;
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

    double Q[cplen];
    double D[cpdim];
    double E[cpdim-1];
    double TAU[cpdim-1];
    double V[cpdim];
    
    eye(cpdim, Q);
    //print_array(N_PROBLEMS, curr_eigenvals);
    LAPACKE_dsytrd(LAPACK_ROW_MAJOR, 'U', cpdim, curr_matrix, cpdim, D, E, TAU);
  

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
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, cpdim, cpdim, cpdim, 1.0, Q, cpdim, Q_curr, cpdim, 1.0, Q, cpdim);
      free(Q_curr);
    }
   
    construct_eigenproblem(&problems[i], cpdim, curr_eigenvals, D, E, Q);
    
  }
  printf("Read %d eigenproblems into the memory.\n", n_problems);

  fclose(f);
  if (curr_line) {
    free(curr_line); 
  }
}

double get_mean(double *array, int len) {
  double res=0.0;
  for(int i=0; i<len; res+=array[i],i++);
  return res/((double) len);
}

void test_dsteqr() {
  int VL, VU, IL, IU;
  struct Eigenproblem problems[N_PROBLEMS];
  load_problems("data.csv", problems);

  double accuracies[N_PROBLEMS];

  for (int i=0;i<N_PROBLEMS;i++) {
    struct Eigenproblem p = problems[i];
    lapack_int info = LAPACKE_dsteqr(LAPACK_ROW_MAJOR, MODE, p.p_size, p.D, p.E, p.Q, p.p_size); 
    
    if(info != 0) {
      printf("Eigenproblem #%d was not solved correctly!\n", i);  
    }
    
    accuracies[i] = get_absolute_accuracy(p.eigenvalues, p.D, p.p_size);
    destroy_eigenproblem(&p);
  }
  printf("Mean accuracy of DSTEQR is %f\n", get_mean(accuracies, N_PROBLEMS));    
}

void test_dstevx() {
  int VL, VU, IL, IU;
  struct Eigenproblem problems[N_PROBLEMS];
  load_problems("data.csv", problems);
  
  double accuracies[N_PROBLEMS];
  for (int i=0;i<N_PROBLEMS;i++) {
    struct Eigenproblem p = problems[i];
    lapack_int M;
    int psize = p.p_size;
    int plen = psize*psize;
    double E[psize-1];
    double D[psize];
    double Z[plen];
    memcpy(E, p.E, sizeof(double)*(psize-1)); 
    memcpy(D, p.D, sizeof(double)*psize); 
    memcpy(Z, p.Q, sizeof(double)*plen);
    double W[psize];
    lapack_int ifail[psize];
    lapack_int info = LAPACKE_dstevx(LAPACK_ROW_MAJOR, MODE, RANGE, psize, D, E, VL, VU, IL, IU, ABSTOL, &M, W, Z, psize, ifail);
    
    if(info != 0) {
      printf("Eigenproblem #%d was not solved correctly!\n", i);  
    }
    
    accuracies[i] = get_absolute_accuracy(p.eigenvalues, W, psize);
    destroy_eigenproblem(&p);
  }
  printf("Mean accuracy of DSTEVX is %f\n", get_mean(accuracies, N_PROBLEMS));    
}

// E here should be N dimensional!
void test_dstemr() {
  int VL, VU, IL, IU;
  struct Eigenproblem problems[N_PROBLEMS];
  load_problems("data.csv", problems);
  double accuracies[N_PROBLEMS];
  for (int i=0;i<N_PROBLEMS;i++) {
    struct Eigenproblem p = problems[i];
    lapack_int M = 0;
    int psize = p.p_size;
    int plen = psize*psize;
    double E[psize];
    double D[psize];
    double Z[plen];
    memcpy(E, p.E, sizeof(double)*(psize-1)); 
    memcpy(D, p.D, sizeof(double)*psize); 
    memcpy(Z, p.Q, sizeof(double)*plen); 
    double W[psize];
    int ifail = 0;
    lapack_int ISUPPZ[psize*2];
    lapack_logical tryrac = (lapack_logical) 0;
    lapack_int info = LAPACKE_dstemr(LAPACK_ROW_MAJOR, MODE, RANGE, psize, D, E, VL, VU, IL, IU, &M, W, Z, psize, psize, ISUPPZ, &tryrac); 

    if(info != 0) {
      printf("Eigenproblem #%d was not solved correctly!\n", i);  
    }

    accuracies[i] = get_absolute_accuracy(p.eigenvalues, W, psize);
    destroy_eigenproblem(&p);   
  }
  printf("Mean accuracy of DSTEMR is %f\n", get_mean(accuracies, N_PROBLEMS));    
}


int main(void) {
  test_dsteqr();
  test_dstevx(); 
  test_dstemr(); 
  return 0;
}
