#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

const int numstep=100;
const int neps=numstep-1;

/* This preprocessed file is located in 'input_data' from the git root.
   If using the provided Makefile from this directory, it should be found. */
#include "monthly_sunspots_data.h"

/* a larger data set */
/* #include "random_2d_data.h" */

double* all_pairs_distances(double* X, int n, int d){
  int i, j, k;
  double xi, xj, tmp;

  double* D = (double*)calloc(n*n, sizeof(double));

  /* There are several ways to parallize this with omp.
     I have written it this way so you may easily remove
     some variables from the private list to generate
     a common error... try removing j, xj ...

     Also, while this section is O(N^2), it doesn't
     yield a great speedup alone.  We need to
     parallelize another quadratic section below.

     There are more intermediate techniques, like
     loop "collapsing" etc, but this is meant to be simple.
  */
#pragma omp parallel private(i, xi, j, xj, k)
#pragma omp for reduction(+:tmp)
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
  int i, j;
  double minD, maxD, val;
  double step;

  /* allocate space for epsilon array */
  double* epsilons = (double*)calloc(neps, sizeof(double));

  /* compute min max of D*/
  /* We expect non diagonal elements to be non zero,
     though this is not "robust".
     We'll also use the symmetry here.
  */
  minD = D[1];
  maxD = D[1];
#pragma omp parallel for private(i, j, val)
  for(i=0; i<n; i++){
    for(j=0; j<i; j++){
      val = D[i*n + j];

      if(val < minD){
	minD = val;
      }

      if(val > maxD){
	maxD = val;
      }
    }
  }

  /* compute a step size */
  step = (maxD - minD) / neps;
  /* printf("step %f\n", step); */
  /* remember, we'll skip the first step, so (i+1), since it would be count of 0...*/
  for(i=0; i<neps; i++){
    epsilons[i] = minD + step*(i+1);
  }

  return epsilons;
}


int* correlation_integrals(double* D, int n, double* epsilons){

  int i, m, cnt;
  double eps;

  int* C = (int*)calloc(neps, sizeof(int));

  /* This section makes better use of shared data,
   and yields a much more dramatic speed increase.*/
#pragma omp parallel private(i)
#pragma omp for reduction(+:cnt)
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


double estimate_dimension(double* epsilons, int* C){
  int i;
  double X[neps];
  double Y[neps];
  double xhat, yhat, num, den;
  double slope;
  /* double inter; */

  write_file(epsilons, C);

  xhat = 0.;
  yhat = 0.;

  for(i=0; i<neps; i++){
    X[i] = logf(epsilons[i]);
    Y[i] = logf(C[i]);

    xhat += X[i];
    yhat += Y[i];
  }

  xhat /= neps;
  yhat /= neps;
  /* printf("xhat yhat %f %f %d\n", xhat, yhat, neps); */

  num = 0.;
  den = 0.;
  for(i=0; i<neps; i++){
    num += X[i] * Y[i];
    den += X[i] * X[i];
  }

  num -= (neps) * xhat * yhat;
  den -= (neps) * xhat * xhat;

  slope = num / den;
  /*inter = yhat - slope * xhat; */

  return slope;

}


int main(int argc, char** argv){
  double* D;
  double* epsilons;
  int* C;

  /* input_data and input_data_n are provided by the pre baked header at TOF. */
  D = all_pairs_distances(input_data, input_data_n, input_data_dim);

  epsilons = generate_epsilons(D, input_data_n);

  C = correlation_integrals(D, input_data_n, epsilons);

  printf("Estimated Correlation Dimension: %f\n",
	 estimate_dimension(epsilons, C));

  free(D);
  free(epsilons);
  free(C);

  return 0;
};
