// 5) Small eigenproblems                                                       
// You are given thousands of small tridiagonal eigenproblems (n<60)           
// Compare the accuracy and speed of three eigensolvers: BX+II, QR, MR3.
// Compare different matrix types:
//   * Random eigenvalues, uniform distribution (0, 1)                            
//   * Uniform eigenvalue distribution (page 119 of http://arxiv.org/pdf/1401.4950v1.pdf
// Instrument the code to count flops. (use PAPI)
// Study flops vs accuracy, for different accuracy levels. (use dlarrv for emr)
// BX+II: DSTEVX, QR: DSTEQR, MR3: DSTEMR  
// 
// 
// ################################################################################
// TODO 3d flps, accuracy, problem size                                           #
// TODO how flops vs accuracy for different accuracy lvls? Only one func has rtol #
// TODO compile with -O0?                                                         # 
// TODO play with TRYRAC for dstemr                                               #
// TODO process correctly if we fail to solve                                     #
// ################################################################################
//
//
//
// example code from here: https://icl.cs.utk.edu/projects/papi/wiki/PAPITopics:Getting_Started#The_Source_Code

#include <papi.h>
#include <cblas.h>
#include <stdio.h>
#include <lapacke.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define DEBUG 0

#define MATRIX_LAYOUT 1
#define MODE 'V'
#define RANGE 'A'
#define N_PROBLEMS 10000
#define ABSTOL 0.0000000000000001

#define DATA_PREFIX "data/"
#define RES_PREFIX "res/"

#define L3_SIZE 4096*1024 // run lscpu to know

struct Eigenproblem {
  int p_size;
  double *eigenvalues;
  // before reduction to symmetric tridiagonal
  // double *matrix;
  // diagonal after reduction to symm tridiagonal
  double *D;
  // subdiagonal after reduction to symm tridiagonal
  double *E;
};

void construct_eigenproblem(struct Eigenproblem *p, int size, double *eigenvalues, double *D, double *E) {
  p->eigenvalues = malloc(sizeof(double)*size);
  p->D = malloc(sizeof(double)*size);
  p->E = malloc(sizeof(double)*(size-1));

  p->p_size = size;
  memcpy(p->eigenvalues, eigenvalues, sizeof(double)*size);
  memcpy(p->D, D, sizeof(double)*size);
  memcpy(p->E, E, sizeof(double)*(size-1));
}

void destroy_eigenproblem(struct Eigenproblem *p) {
  free(p->eigenvalues);
  free(p->D);
  free(p->E);
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
    printf("%lf ", array[i]);
  }
  printf("\n");
}


double get_accuracy(int dim, double *V){
  // |V'V-I|
  double res[dim*dim];
  eye(dim, res);
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, dim, dim, dim, 1.0, V, dim, V, dim, -1.0, res, dim);
  double max = -1;
  for(int i=0;i<dim*dim;i++) {
    if (fabs(res[i]) > max){
      max = fabs(res[i]);
    }
  }
  return max;
}
 
void compile_accuracy_speed_flops(double *accuracy, double *speed, long long *flops, char res[][256]) {
  for (int i = 0; i< N_PROBLEMS; i++) {
    snprintf(res[i], sizeof res[i], "%.20lf;%.20lf;%lld\n", accuracy[i], speed[i], flops[i]);
  }
}


void write_results(char* filename, char results[][256], int nb_res, char *exp_type, char *method, char *accuracy_type){
  char path[256];
  snprintf(path, sizeof path, "%s%s/%s_%s_%s", RES_PREFIX, exp_type, method, accuracy_type, filename);
  FILE *f = fopen(path, "w");
  for(int i=0;i<nb_res;i++){
    fprintf(f, "%s", results[i]);
  }
  fclose(f);
}

void load_problems(char *filename, struct Eigenproblem *problems) {
  char path[256];
  snprintf(path, sizeof path, "%s%s", DATA_PREFIX, filename);
  FILE *f = fopen(path, "r");
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

    double D[cpdim];
    double E[cpdim-1];
    double TAU[cpdim-1];
    double V[cpdim];
    
    //print_array(N_PROBLEMS, curr_eigenvals);
    LAPACKE_dsytrd(LAPACK_ROW_MAJOR, 'U', cpdim, curr_matrix, cpdim, D, E, TAU);
    construct_eigenproblem(&problems[i], cpdim, curr_eigenvals, D, E);
    
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

struct Timing {
  float real_time, proc_time, mflops;
  long long flpins;
};

static void test_fail(char *file, int line, char *call, int retval){
    printf("%s\tFAILED\nLine # %d\n", file, line);
    if ( retval == PAPI_ESYS ) {
        char buf[128];
        memset( buf, '\0', sizeof(buf) );
        sprintf(buf, "System error in %s:", call );
        perror(buf);
    }
    else if ( retval > 0 ) {
        printf("Error calculating: %s\n", call );
    }
    else {
        printf("Error in %s: %s\n", call, PAPI_strerror(retval) );
    }
    printf("\n");
    exit(1);
}

void call_PAPI(struct Timing *t){
  int retval;
  if(retval=PAPI_flops(&t->real_time, &t->proc_time, &t->flpins, &t->mflops)<PAPI_OK){
    test_fail(__FILE__, __LINE__, "PAPI_flops", retval);  
  }
}

void test_dsteqr(char *filename, char *exp_type, char *accuracy_type) {
  printf("Performing tests for DSTEQR\n");
  int VL, VU, IL, IU;
  struct Eigenproblem problems[N_PROBLEMS];
  load_problems(filename, problems);
  double accuracies[N_PROBLEMS];
  double real_time[N_PROBLEMS];
  long long flops[N_PROBLEMS];
  
  for (int i=0;i<N_PROBLEMS;i++) {
    //garbage initialization
    int n_elem = L3_SIZE/sizeof(double);
    double garbage[n_elem];
    for(int g=0; g<n_elem;g++){
      garbage[g]+=0.3;
    }
    
    struct Eigenproblem p = problems[i];
    double Z[p.p_size*p.p_size];
    eye(p.p_size, Z);
   
    //START TIMING 
    struct Timing t;
    call_PAPI(&t); 
    lapack_int info = LAPACKE_dsteqr(LAPACK_ROW_MAJOR, 'I', p.p_size, p.D, p.E, Z, p.p_size); 
    call_PAPI(&t);
    flops[i] = t.flpins;
    real_time[i] = t.real_time;
    PAPI_shutdown();
    //END TIMING

    accuracies[i] = get_accuracy(p.p_size, Z);
    if(info != 0) {
      printf("Eigenproblem #%d was not solved correctly!\n", i);  
      //accuracies[i] = accuracies[i-1];
    }
    destroy_eigenproblem(&p);
  }
  printf("Mean accuracy of DSTEQR is %.20lf\n", get_mean(accuracies, N_PROBLEMS));    
  char result[N_PROBLEMS][256];
  compile_accuracy_speed_flops(accuracies, real_time, flops, result);
  write_results(filename , result, N_PROBLEMS, exp_type, "dsteqr", accuracy_type);
  printf("DONE.\n");
}

void test_dstevx(char *filename, char *exp_type, double tolerance, char *accuracy_type) {
  printf("Performing tests for DSTEVX\n");
  int VL, VU, IL, IU;
  struct Eigenproblem problems[N_PROBLEMS];
  load_problems(filename, problems);
  double accuracies[N_PROBLEMS];
  double real_time[N_PROBLEMS];
  long long flops[N_PROBLEMS];
    
  
  for (int i=0;i<N_PROBLEMS;i++) {
    //garbage initialization
    int n_elem = L3_SIZE/sizeof(double);
    double garbage[n_elem];
    for(int g=0; g<n_elem;g++){
      garbage[g]+=0.3;
    }

    struct Eigenproblem p = problems[i];
    lapack_int M;
    int psize = p.p_size;
    int plen = psize*psize;
    double E[psize-1];
    double D[psize];
    double Z[plen];
    memcpy(E, p.E, sizeof(double)*(psize-1)); 
    memcpy(D, p.D, sizeof(double)*psize); 
    double W[psize];
    lapack_int ifail[psize];
    
    //START TIMING
    struct Timing t;
    call_PAPI(&t); 
    lapack_int info = LAPACKE_dstevx(LAPACK_ROW_MAJOR, MODE, RANGE, psize, D, E, VL, VU, IL, IU, tolerance, &M, W, Z, psize, ifail);
    call_PAPI(&t);
    flops[i] = t.flpins;
    real_time[i] = t.real_time;
    PAPI_shutdown();
    //END TIMING

    accuracies[i] = get_accuracy(p.p_size, Z);
    if(info != 0) {
      printf("Eigenproblem #%d was not solved correctly!\n", i);  
      //accuracies[i] = accuracies[i-1];
    }

    destroy_eigenproblem(&p);
  }
  printf("Mean accuracy of DSTEVX is %.20lf\n", get_mean(accuracies, N_PROBLEMS));    
  char result[N_PROBLEMS][256];
  compile_accuracy_speed_flops(accuracies, real_time, flops, result);
  write_results(filename , result, N_PROBLEMS, exp_type, "dstevx", accuracy_type);
  printf("DONE.\n");
}

// E here should be N dimensional!
void test_dstemr(char *filename, char *exp_type, lapack_logical tryrac, char *accuracy_type) {
  printf("Performing tests for DSTEMR\n");
  int VL, VU, IL, IU;
  struct Eigenproblem problems[N_PROBLEMS];
  load_problems(filename, problems);
  double accuracies[N_PROBLEMS];
  double real_time[N_PROBLEMS];
  long long flops[N_PROBLEMS];
  
  for (int i=0;i<N_PROBLEMS;i++) {
    //garbage initialization
    int n_elem = L3_SIZE/sizeof(double);
    double garbage[n_elem];
    for(int g=0; g<n_elem;g++){
      garbage[g]+=0.3;
    }
    
    struct Eigenproblem p = problems[i];
    lapack_int M = 0;
    int psize = p.p_size;
    int plen = psize*psize;
    double E[psize];
    double D[psize];
    double Z[plen];
    memcpy(E, p.E, sizeof(double)*(psize-1)); 
    memcpy(D, p.D, sizeof(double)*psize); 
    double W[psize];
    int ifail = 0;
    lapack_int ISUPPZ[psize*2];
    
    //START TIMING 
    struct Timing t;
    call_PAPI(&t); 
    lapack_int info = LAPACKE_dstemr(LAPACK_ROW_MAJOR, MODE, RANGE, psize, D, E, VL, VU, IL, IU, &M, W, Z, psize, psize, ISUPPZ, &tryrac); 
    call_PAPI(&t);
    flops[i] = t.flpins;
    real_time[i] = t.real_time;
    PAPI_shutdown();
    //END TIMING

    accuracies[i] = get_accuracy(p.p_size, Z);
    if(info != 0) {
      printf("Eigenproblem #%d was not solved correctly!\n", i);  
      //accuracies[i] = accuracies[i-1];
    }

    accuracies[i] = get_accuracy(p.p_size, Z);
    destroy_eigenproblem(&p);   
  }
  printf("Mean accuracy of DSTEMR is %.20lf\n", get_mean(accuracies, N_PROBLEMS));    
  char result[N_PROBLEMS][256];
  compile_accuracy_speed_flops(accuracies, real_time, flops, result);
  write_results(filename , result, N_PROBLEMS, exp_type, "dstemr", accuracy_type);
  printf("DONE.\n");
}


int main(int argc, char **argv) {
  char path[256];
  if (strcmp(argv[1], "speed-vs-accuracy") == 0) {
    test_dsteqr(argv[2], "speed-vs-accuracy", "");
    test_dstevx(argv[2], "speed-vs-accuracy", ABSTOL, "");
    test_dstemr(argv[2], "speed-vs-accuracy", (lapack_logical) 0, "");
  } else if(strcmp((argv[1]), "flops-given-accuracy") == 0) {
    // NO TEST FOR DSTEQR since it has no nolerance parameter
    for(int i = 1; i < 15; i++){
      char accstr[32];
      snprintf(accstr, sizeof(accstr), "%.20lf", 1.0/pow(10.0, i));
      test_dstevx(argv[2], "flops-given-accuracy", 1.0/pow(10.0,i), accstr);
    }
    //DSTEMR test using tryrac: TRUE if high accuracy, FALSE otherwise
    test_dstemr(argv[2], "flops-given-accuracy", (lapack_logical) 0, "low");
    test_dstemr(argv[2], "flops-given-accuracy", (lapack_logical) 1, "high");
  } else {
    printf("Unknown test type parameter\n");
  }
  return 0;
}
