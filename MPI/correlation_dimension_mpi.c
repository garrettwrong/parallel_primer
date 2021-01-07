#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

const int numstep=100;
const int neps=numstep-1;

/* This preprocessed file is located in 'input_data' from the git root.
   If using the provided Makefile from this directory, it should be found. */
#include "monthly_sunspots_data.h"

double* all_pairs_distances(double* X, int n, int d, int pid, int numprocs){
  int i, j, k;
  double xi, xj, tmp;
  int nelem, ubnd;

  /* Note we'll "cheat" and make this a little larger to avoid boundary issues.
     This is fine so long as we ignore pad values in our loops.*/
  nelem = (n + numprocs - 1) / numprocs;
  double* D = (double*)calloc(n*n + numprocs - 1, sizeof(double));
  ubnd = nelem*(pid+1);
  if(n < ubnd){
    ubnd = n;
  }
  /* an alternative to this would be to use strides of numprocs.. */

  for(i=nelem*pid; i<ubnd; i++){
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
  step = (maxD - minD) / neps;
  /* printf("step %f\n", step); */
  /* remember, we'll skip the first step, so (i+1), since it would be count of 0...*/
  for(i=0; i<neps; i++){
    epsilons[i] = minD + step*(i+1);
  }

  return epsilons;
}


int* correlation_integrals(double* D, int n, double* epsilons, int pid, int numprocs){

  int i, m, cnt;
  double eps;
  int nelem, ubnd;

  int* C = (int*)calloc(neps + numprocs - 1, sizeof(int));

  nelem = (neps + numprocs - 1) / numprocs;
  ubnd = (pid+1)*nelem;
  if(neps<ubnd){
    ubnd = neps;
  }

  for(i=pid*nelem; i<ubnd; i++){
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
  int bufC[neps+32];
  int* C;


  /* MPI Initialization */
  int mpi_pid, numprocs, nelem;
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_pid);


  /* monthly_sunspots and monthly_sunspots_n are provided by the pre baked header at TOF. */

  /* Each process will do roughly 1/numprocs of the computation D. */
  D = all_pairs_distances(monthly_sunspots, monthly_sunspots_n, 1, mpi_pid, numprocs);

  /* Then we'll share the results for each piece of D. */
  nelem = (monthly_sunspots_n*monthly_sunspots_n + numprocs - 1) / numprocs;
  MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, D, nelem, MPI_DOUBLE, MPI_COMM_WORLD);

  /* This is function is cheap, all processes can do it;
     not worth additional code for abroadcast. */
  epsilons = generate_epsilons(D, monthly_sunspots_n);

  /* Each process will do roughly 1/numprocs computation of C. */
  C = correlation_integrals(D, monthly_sunspots_n, epsilons, mpi_pid, numprocs);

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

  /* MPI shutdown */
  MPI_Finalize();

  return 0;

};
