#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int numstep=100;
const int neps=numstep-1;

/* This preprocessed file is located in 'input_data' from the git root.
   If using the provided Makefile from this directory, it should be found. */
#include "monthly_sunspots_data.h"

/* a larger data set */
/* #include "random_2d_data.h" */

__global__ void all_pairs_distances_kernel(double* X, int n, int d, double* D){
  int i, k;
  double xi, xj, tmp;

  const int tidx = blockDim.x * blockIdx.x + threadIdx.x;
  if( tidx >= n){
    return;
  }

  for(i=0; i<n; i++){

    /* compute norm in dim d,
       this is sort of stupid for d=1, alas.
    */
    tmp = 0;
    for(k=0; k<d; k++){
      xi = X[d*i + k];
      xj = X[d*tidx + k];
      /* square */
      tmp += (xi-xj) * (xi-xj);
    }
    D[i*n + tidx] = sqrtf(tmp);
  }
}



double* all_pairs_distances(double* X, int n, int d){
  double* D;
  cudaMalloc(&D, n * n * sizeof(double));

  double* X_dev;
  cudaMalloc(&X_dev, n * d * sizeof(double));
  cudaMemcpy(X_dev, X, n * d * sizeof(double), cudaMemcpyHostToDevice);

  int blocksz = 1024;
  int nblocks = (n + blocksz -1 ) / blocksz;

  all_pairs_distances_kernel<<<nblocks, blocksz>>>(X_dev, n, d, D);

  return D;
}

__global__ void min_max_kernel(double* D, int n, double* buf){
  int i;
  double minD, maxD, val;

  const int tidx = blockDim.x * blockIdx.x + threadIdx.x;
  if(tidx >= n){
    return;
  }

  /* compute min max of D*/
  /* We expect non diagonal elements to be non zero,
     though this is not "robust".
     I think using symmetry here would actually hurt performance...
     do you know why?
  */
  minD = D[1];
  maxD = D[1];

  for(i=0; i<n; i++){
    /* skip diagonals */
    if(tidx == i){
      continue;
    }
    val = D[i*n + tidx];

    if(val<minD){
      minD = val;
    }

    if(val>maxD){
      maxD = val;
    }
  }
  buf[tidx] = minD;
  buf[n + tidx] = maxD;

}

__global__ void generate_epsilons_kernel(double* D, int n, double* epsilons_dev, double minD, double maxD){
  double step;

  const int tidx = blockDim.x * blockIdx.x + threadIdx.x;
  if(tidx >= neps){
    return;
  }

  /* compute a step size */
  step = (maxD - minD) / neps;
  /* printf("step %f\n", step); */
  /* remember, we'll skip the first step, so (i+1), since it would be count of 0...*/
  epsilons_dev[tidx] = minD + step*(tidx+1);
}

double* generate_epsilons(double* D, int n){
  int i;
  double minD, maxD;

  /* compute min max of D*/
  double* buf = (double*)calloc(2*n, sizeof(double));
  double* buf_dev;
  cudaMalloc(&buf_dev, 2*n*sizeof(double));

  int blocksz = 1024;
  int nblocks = (n + blocksz -1 ) / blocksz;

  min_max_kernel<<<blocksz, nblocks>>>(D, n, buf_dev);
  cudaMemcpy(buf, buf_dev, 2*n*sizeof(double), cudaMemcpyDeviceToHost);
  minD = buf[0];
  maxD = buf[n];
  for(i=1; i<n; i++){
    if(buf[i]<minD){
      minD = buf[i];
    }

    if(buf[n+i]>maxD){
      maxD = buf[n+i];
    }
  }
  /* printf("minxD %f maxD %f\n", minD, maxD); */

  /* allocate space for epsilon array */
  double* epsilons_dev;
  cudaMalloc(&epsilons_dev, neps*sizeof(double));

  blocksz = 1024;
  nblocks = (neps + blocksz -1 ) / blocksz;

  generate_epsilons_kernel<<<blocksz, nblocks>>>(D, n, epsilons_dev, minD, maxD);

  cudaFree(buf_dev);
  free(buf);

  return epsilons_dev;
}

__global__ void correlation_integrals_kernel(double* D, int n, double* epsilons, int* C){
  int i, cnt;
  double eps;

  const int tidx = blockDim.x * blockIdx.x + threadIdx.x;
  const int tidy = blockDim.y * blockIdx.y + threadIdx.y;
  if( tidx >= n){
    return;
  }
  if( tidy >= neps){
    return;
  }

  eps = epsilons[tidy];
    /* loop through D, counting if closer than eps */
  cnt = 0;
  for(i=0; i<n; i++){
    if(D[i*n + tidx] < eps) {
      cnt++;
    }
    /* assign */
  }

  /* Each thread writes its partial count to a scratch buffer, */
  C[neps*tidx + tidy] = cnt;

}

__global__ void correlation_integrals_redux(int* C, int n){
  int i;

  const int tidx = blockDim.x * blockIdx.x + threadIdx.x;
  if( tidx >= neps){
    return;
  }
  /* then we'll sum that scratch buffer. */
  for(i=1; i<n; i++){
    C[tidx] += C[neps*i + tidx];
  }
}


int* correlation_integrals(double* D, int n, double* epsilons){

  /* we'll need to make some space and transfer the epsilons */
  double* epsilons_dev;
  cudaMalloc(&epsilons_dev, neps*sizeof(double));
  cudaMemcpy(epsilons_dev, epsilons, neps*sizeof(double), cudaMemcpyHostToDevice);

  /* note C_dev is a much larger buffer for reduction sum scratch space*/
  int* C = (int*)calloc(neps, sizeof(int));
  int* C_dev;
  cudaMalloc(&C_dev, neps*n*sizeof(int));

  dim3 blocksz(128,128);
  dim3 nblocks((n + blocksz.x - 1 ) / blocksz.x,
	       (neps + blocksz.y - 1 ) / blocksz.y);

  correlation_integrals_kernel<<<blocksz, nblocks>>>(D, n, epsilons_dev, C_dev);
  correlation_integrals_redux<<<blocksz.y, nblocks.y>>>(C_dev, n);

  /* copy results to host*/
  cudaMemcpy(C, C_dev, neps*sizeof(int), cudaMemcpyDeviceToHost);

  cudaFree(C_dev);

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


double estimate_dimension(double* epsilons_dev, int* C){
  int i;
  double* X = (double*)calloc(neps, sizeof(double));
  double* Y = (double*)calloc(neps, sizeof(double));
  double xhat, yhat, num, den;
  double slope;
  /* double inter; */

  double* epsilons = (double*)calloc(neps, sizeof(double));
  cudaMemcpy(epsilons, epsilons_dev, neps*sizeof(double), cudaMemcpyDeviceToHost);

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
  /* printf("xhat yhat %f %f %d\n", xhat, yhat, n); */

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
  double* D_dev;
  double* epsilons_dev;
  int* C;

  /* input_data and input_data_n are provided by the pre baked header at TOF. */
  D_dev = all_pairs_distances(input_data, input_data_n, input_data_dim);

  epsilons_dev = generate_epsilons(D_dev, input_data_n);

  C = correlation_integrals(D_dev, input_data_n, epsilons_dev);

  printf("Estimated Correlation Dimension: %f\n",
	 estimate_dimension(epsilons_dev, C));

  cudaFree(D_dev);
  cudaFree(epsilons_dev);
  free(C);

  return 0;
};
