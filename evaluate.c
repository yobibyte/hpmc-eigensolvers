//TODO read csv with maximum precision
// TODO check for memleaks

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct Eigenproblem {
  int p_size;
  double *eigenvalues;
  double *matrix;
};

int main(void) {

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
  
  for(int i=0; i<n_problems;i++) {
    printf("Reading eigenproblem #%d with dim %s\n", i, curr_buf);
    
    read = getline(&curr_line, &len, f);
    char *end_str = NULL;
    curr_buf        = strtok_r(curr_line, ";", &end_str);
    eigenvalues_str = strtok_r(NULL,      ";", &end_str);
    curr_matrix_str = strtok_r(NULL,      ";", &end_str); 
    
    int cpsize;
    sscanf(curr_buf, "%d", &cpsize); 
 
    char *eigenvalues_rem;
    curr_eigenval = strtok_r(eigenvalues_str, ",", &eigenvalues_rem);
    double curr_eigenvals[cpsize];
    sscanf(curr_eigenval, "%lf", &curr_eigenvals[0]);
    printf("%f\t", curr_eigenvals[0]); 
    for(int j=1;j<cpsize;j++){
      curr_eigenval     = strtok_r(NULL, ",", &eigenvalues_rem);
      sscanf(curr_eigenval, "%lf", &curr_eigenvals[j]);
      printf("%f\t", curr_eigenvals[j]);
    }    
    printf("\n");
  
    char *curr_elem;
    char *cur_mat_rem;
    curr_elem = strtok_r(curr_matrix_str, ",", &cur_mat_rem);
    double curr_matrix[cpsize*cpsize];
    sscanf(curr_elem, "%lf", &curr_matrix[0]);
    printf("%lf ", curr_matrix[0]); 
    for(int j=1;j<cpsize*cpsize;j++){
      curr_elem= strtok_r(NULL, ",", &cur_mat_rem);
      sscanf(curr_elem, "%lf", &curr_matrix[j]);
      printf("%lf ", curr_matrix[j]);
    }    
    printf("\n");

    
    struct Eigenproblem p;
    p.p_size = cpsize;  
    p.matrix = curr_matrix;    
    p.eigenvalues = curr_eigenvals;  
  
  }
  
  printf("Read %d eigenproblems into the memory.\n", n_problems);

  fclose(f);
  if (curr_line) {
    free(curr_line); 
  } 
}
