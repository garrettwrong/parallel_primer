#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int neps=100;

/* This preprocessed file is located in 'input_data' from the git root.
   If using the provided Makefile from this directory, it should be found. */
#include "monthly_sunspots_data.h"

double* all_pairs_distances(double* X, int n, int d){
  int i, j, k;
  double xi, xj, tmp;
  
  double* D = (double*)calloc(n*n, sizeof(double));
  
  for(i=0; i<n; i++){
    for(j=0; j<n; j++){
      /* compute norm in dim d,
	 this is sort of stupid for d=1, alas.
       */
      tmp = 0;
      for(k=0; k<d; k++){
	xi = X[d*i + k];
	xj = X[d*j + k];
	/* square */
	tmp += (xi-xj) * (xi-xj);
      }
      D[i*n + j] = sqrtf(tmp);
    }
  }

  return D;
}


double* generate_epsilons(double* D, int n){
  int i;
  double minD, maxD, val;
  double step;

  /* allocate space for epsilon array */
  double* epsilons = (double*)calloc(neps, sizeof(double));  

  /* compute min max of D*/
  minD = D[0];
  maxD = D[0];   
  for(i=1; i<n*n; i++){
    val = D[i];    

    if(val < minD){
      minD = val;
    }

    if(val > maxD){
      maxD = val;
    }
  }

  /* compute a step size */
  step = (maxD - minD) / (neps+1);
  /* but skip the first, since it would be count of 0...*/
  for(i=0; i<neps; i++){
    epsilons[i] = minD + step*(i+1);
  }

  return epsilons;
}

int* correlation_integrals(double* D, int n, double* epsilons){

  int i, m, cnt;
  double eps;
  
  int* C = (int*)calloc(neps, sizeof(int));
  
  for(i=0; i<neps; i++){
    cnt = 0;
    eps = epsilons[i];
    /* loop through D, counting if closer than eps */
    for(m=0; m<n*n; m++){
      if(D[m] < eps) cnt++;
    }
    /* assign */
    C[i] = cnt;
  }

  return C;
}

void write_file(double* epsilons, int* C){
  int n;

  FILE* fh = fopen("correlation_integrals.dat", "w");
  
  for(n=0; n<neps; n++){
    fprintf(fh, "%f %d\n", epsilons[n], C[n]);
  }

  fclose(fh);
}


int main(int argc, char** argv){
  double* D;
  double* epsilons;
  int* C;

  /* monthly_sunspots and monthly_sunspots_n are provided by the pre baked header at TOF. */
  D = all_pairs_distances(monthly_sunspots, monthly_sunspots_n, 1);

  epsilons = generate_epsilons(D, monthly_sunspots_n);

  C = correlation_integrals(D, monthly_sunspots_n, epsilons);

  write_file(epsilons, C);

  free(D);
  free(epsilons);
  free(C);
  
  return 0;
};
