# Parallel Programming Primer Examples

## About

For these example programs we will compute "Correlation Dimension".  This implementation is trivially O(N^2) steps and easily subdivided in various ways.  Fractional dimension measure is sort of a fun idea.  It is a cute stepping stone into more formal non linear time series analysis. Something we can understand and code in a few minutes... mostly.

The high level steps to the algorithm are simple:

Given a set of N points in R^d,

1) Compute all pairs of distances D between two points n_i, n_j.  We will use Euclidean distance for appropriate input data dimension (l2 norm).

2) For a range of epsilons, we count how many distances in D are less than each epsilon. This is usually called the "correlation sum", up to normalization constants.

3) Fit a line to the log-log scatter plot of counts vs epsilon. The slope of this line is the estimate of Correlation Dimension.

Because we don't have infinite points, the line is only an estimate;
fitting it requires interpretation.
For this toy example, we'll simply truncate the tail effect before fitting.
The python code should generate a plot saved as png if you are curious.

This is not how it is done in practice,
but seems to make an example that is
a little more intersting than dot product.

## Navigating this Repo

Each implementation is located in it's own directory.
The serial version was prototyped in Python.
It is intentionally not optimized beyond using numpy.
There is _unoptimized python version that uses native python loops (and is about 100x slower, yikes).
This version was mostly translated to C, again not optimized, just a basic port.
The C variant was then minimally adapted to other approaches.

Currently:

* Python
* OpenMP (gomp)
* MPI (OpenMPI)
* CUDA
* CUPY

There is an example environment that I have used on TigerCPU and TigerGPU.
I'd encourage you to look at the contents before sourcing.  You may prefer
or require something different.

`source tiger.environ`

## License

Released as GPLv3.

This code was written as supporting material for a short Parallel Computing Primer workshop,
as part Princeton University's larger "Winter Break 2021 Research Computing Bootcamp".

Copyright Garrett Wright 2021.
