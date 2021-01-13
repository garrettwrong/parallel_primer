#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

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
  int nelem, ubnd;

  /* MPI setup*/
  int mpi_pid=0;
  int numprocs=1;
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_pid);


  /* Note we'll "cheat" and make this a little larger to avoid boundary issues.
     This is fine so long as we ignore pad values in our loops.*/
  nelem = (n + numprocs - 1) / numprocs;
  double* D = (double*)calloc(n*n + numprocs - 1, sizeof(double));
  ubnd = nelem*(mpi_pid+1);
  if(n < ubnd){
    ubnd = n;
  }
  /* an alternative to this would be to use strides of numprocs.. */

  for(i=nelem*mpi_pid; i<ubnd; i++){
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
  for(i=0; i<n; i++){
    for(j=0; j<i; j++){
      val = D[i*n +j];

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
  int nelem, ubnd;

  /* MPI setup*/
  int mpi_pid=0;
  int numprocs=1;
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_pid);

  int* C = (int*)calloc(neps + numprocs - 1, sizeof(int));

  nelem = (neps + numprocs - 1) / numprocs;
  ubnd = (mpi_pid+1)*nelem;
  if(neps<ubnd){
    ubnd = neps;
  }

  for(i=mpi_pid*nelem; i<ubnd; i++){
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

  /* Since we don't have a real limit situation here,
     we'll truncate the tail of this dataset. */
  int n = (int)(0.5*neps);

  xhat = 0.;
  yhat = 0.;
  for(i=0; i<n; i++){
    X[i] = logf(epsilons[i]);
    Y[i] = logf(C[i]);

    xhat += X[i];
    yhat += Y[i];
  }

  xhat /= n;
  yhat /= n;
  /* printf("xhat yhat %f %f %d\n", xhat, yhat, neps); */

  num = 0.;
  den = 0.;
  for(i=0; i<n; i++){
    num += X[i] * Y[i];
    den += X[i] * X[i];
  }

  num -= n * xhat * yhat;
  den -= n * xhat * xhat;

  slope = num / den;
  /*inter = yhat - slope * xhat; */

  return slope;

}


int main(int argc, char** argv){
  double* D;
  double* epsilons;
  int* bufC;
  int* C;

  /*
    Note that it is pretty common for folks to use #ifdef
    to switch MPI code in an out, but this complicates things,
    perhaps too much for a basic example.
  */

  /* MPI Initialization */
  int mpi_pid, numprocs, nelem;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_pid);
  bufC = (int*)calloc(neps+numprocs-1, sizeof(int));

  /* input_data and input_data_n are provided by the pre baked header at TOF. */

  /* Each process will do roughly 1/numprocs of the computation D. */
  D = all_pairs_distances(input_data, input_data_n, input_data_dim);

  /* Then we'll share the results for each piece of D. */
  nelem = (input_data_n + numprocs - 1) / numprocs;
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, D,
		nelem*input_data_n, MPI_DOUBLE, MPI_COMM_WORLD);

  /* This is function is cheap, all processes can do it;
     not worth additional code for a broadcast. */
  epsilons = generate_epsilons(D, input_data_n);

  /* Each process will do roughly 1/numprocs computation of C. */
  C = correlation_integrals(D, input_data_n, epsilons);

  /* Then we'll gather the results for C on mpi_pid=0 */
  nelem = (neps + numprocs - 1) / numprocs;
  /* Note the send address is tied to process here*/
  MPI_Gather(&C[mpi_pid*nelem], nelem, MPI_INT, bufC, nelem, MPI_INT, 0, MPI_COMM_WORLD);
  /* don't forget to copy root's values to buffer */
  for(int i=0; i<nelem; i++){
    bufC[i] = C[i];
  }

  /* Only the main process needs to compute this. */
  if(mpi_pid == 0){
    printf("Estimated Correlation Dimension: %f\n",
	   estimate_dimension(epsilons, bufC));
  }

  free(D);
  free(epsilons);
  free(C);
  free(bufC);

  /* MPI shutdown */
  MPI_Finalize();

  return 0;

};
