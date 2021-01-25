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

int get_chunk_size(int n, int nchunks){
  /* Get the number of elements required
     in a chunk to cover n by nchunks. */
  return (n + nchunks - 1) / nchunks;
}

double* all_pairs_distances(double* X, int n, int d){
  int i, j, k, global_i;
  double xi, xj, tmp;
  int chunk_size;

  /* MPI setup*/
  int mpi_pid;
  int numprocs;
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_pid);


  /* Note we'll "cheat" and make this a little larger to avoid boundary issues.
     This is fine so long as we ignore pad values in our loops. */
  chunk_size = get_chunk_size(n, numprocs);

  /* Each of the chunk_size will contribute n distances,
     so in total this rank is responsible for chunk_size*n. */
  double* D = (double*)calloc(chunk_size*n, sizeof(double));

  /* We'll compute the exact elements we're responsible for now */
  for(i=0; i<chunk_size; i++){
    /* I'll check boundary condition,
       which covers the common case where the number or processors
       does not evenly divide n.
       An alternative to this would be to use strides of numprocs,
       with different index arith/checks.
    */

    /*  Note it is absoultely ciritical we understand when
	to use local i index into a chunk,
	versus when we care about the global_i index,
	such as accessing X.
    */
    global_i = mpi_pid * chunk_size + i;
    if( global_i >= n){
      /* Out of bounds (we're past n on the last chunk). */
      break;
    }

    for(j=0; j<n; j++){
      /* Compute norm in dim d;
	 this is sort of stupid for d=1, alas.
       */
      tmp = 0;
      for(k=0; k<d; k++){
	xi = X[d*global_i + k];
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

  /* MPI setup*/
  int mpi_pid=0;
  int numprocs=1;
  int Dn;
  double global_minD, global_maxD;
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_pid);

  /* compute min max of D*/

  /* We expect non diagonal elements to be non zero,
     though this is not "robust".
     Here, to find non diagonal, we can't rely on D[1],
     since D is a chunk which might represent the row where
     this element is a diagonal.
  */


  int chunk_size = get_chunk_size(n, numprocs);
  minD = D[mpi_pid*chunk_size + 1];
  maxD = D[mpi_pid*chunk_size + 1];

  for(i=0; i<chunk_size; i++){
    /*     We'll also use the symmetry here as before,
	   but the triangular index calculation needs to
	   account for D only representing a chunk of a larger
	   matrix that is _distributed_ amoung the ranks.
    */
    Dn = mpi_pid * chunk_size + i;
    if(Dn>n){
      /* Should look familiar, tail end out of bounds wrt n. */
      break;
    }

    for(j=0; j<Dn; j++){
      /* You may note this is not ideal because
	 ranks have asymetric work loads.

	 That is some might have a small section of triangle,
	 while others have something on the order of chunk_size * n.

	 In practice you might combine the longest & shortest,
	 second longest & second shortest, and so on, workloads.
	 One might also get acceptable balance using the strided
	 approach.

	 This would complicate the program further though...

	 In profiling this would probably reveal itself as
	 a synchonization or rank time/utilization imbalance.
	 Then you might make a judgement call.
       */
      val = D[i*n +j];

      if(val < minD){
	minD = val;
      }

      if(val > maxD){
	maxD = val;
      }
    }
  }

  /*
     At this point, all processes have local minD and maxD for their chunk.
     We'll need to aggregate them so each process can compute the epsilon bins.
  */
  MPI_Allreduce(&minD, &global_minD, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&maxD, &global_maxD, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  /* compute a step size */
  step = (global_maxD - global_minD) / neps;
  printf("local min max %f %f \t global min max %f %f\n", minD, maxD, global_minD, global_maxD);
  printf("step %f\n", step);
  /* remember, we'll skip the first step, so (i+1), since it would be count of 0...*/
  for(i=0; i<neps; i++){
    epsilons[i] = global_minD + step*(i+1);
  }

  /* All processes should have their own (identical),
     copy of epsilons. */
  return epsilons;
}


int* correlation_integrals(double* D, int n, double* epsilons){

  int i, j, m, cnt;
  double eps;

  /* MPI setup*/
  int mpi_pid=0;
  int numprocs=1;
  int chunk_size;
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_pid);

  /* Each process will count their contribution to each bins' sum,
     and again, we'll exchange to aggregate the result.

     That is because we can't simply partition the binning of C,
     since each process only has a chunk of D.  To avoid communicating the
     entire D matrix, at any point, we'll instead calculate and communicate
     partial counts in C in each process.

     This is a common refactor that comes up in MPI programs.
  */
  int* C = (int*)calloc(neps, sizeof(int));

  chunk_size = get_chunk_size(n, numprocs);

  for(i=0; i<neps; i++){
    cnt = 0;
    eps = epsilons[i];
    /* loop through our chunk of D, counting if closer than eps */
    for(m=0; m<chunk_size; m++){
      /* if we are still in bounds.. */
      if(mpi_pid * chunk_size + m < n){
	for(j=0; j<n; j++){
	  if(D[m*n + j] < eps) cnt++;
	}
      }
    }
    /* assign */
    C[i] = cnt;
  }

  /*
    Cool, so now each process has their own counts stored locally in C.

     We now need to sum C. Unlike before, we only need this on the
     main (root, 0) process.
   */
  int* global_C = (int*)calloc(neps, sizeof(int));

  MPI_Reduce(C, global_C, neps, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  free(C);

  return global_C;
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
  double* X = (double*)calloc(neps, sizeof(double));
  double* Y = (double*)calloc(neps, sizeof(double));
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

  free(X);
  free(Y);

  return slope;

}


int main(int argc, char** argv){
  double* D;
  double* epsilons;
  int* C;

  /*
    Note that it is pretty common for folks to use #ifdef
    to switch MPI code in an out, but this complicates things,
    perhaps too much for a basic example.
  */

  /* MPI Initialization */
  int mpi_pid, numprocs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_pid);

  /* input_data and input_data_n are provided by the pre baked header at TOF. */

  /* Each process will do roughly 1/numprocs of the computation D.
     The exact chunk size is displayed.
  */
  if(mpi_pid == 0){
    printf("Partitioning %d elements into %d ranks, " \
	   "each responsible for %d elements.\n",
	   input_data_n,
	   numprocs,
	   get_chunk_size(input_data_n, numprocs));
  }
  D = all_pairs_distances(input_data, input_data_n, input_data_dim);

  /* All processes need the result of the following function.
     It is modified such that each process
     computes their chunk of D known from above.

     Then they will exchange information required,
     namely the off diagonal min and max of D.

     Note throughout all of these processes, we never construct the entire
     n*n D matrix on any single process.  This distribution of memory is
     one of the core strengths of MPI approaches.
  */
  epsilons = generate_epsilons(D, input_data_n);

  /* Each process will contribute a portion of C
     by scanning their chunk of D.

     The results are then aggregated (reduction sum),
     forming the global C result.
   */
  C = correlation_integrals(D, input_data_n, epsilons);

  /* Only the main process needs to compute this. */
  if(mpi_pid == 0){
    printf("Estimated Correlation Dimension: %f\n",
	   estimate_dimension(epsilons, C));
  }

  free(D);
  free(epsilons);
  free(C);

  /* MPI shutdown */
  MPI_Finalize();

  return 0;

};
