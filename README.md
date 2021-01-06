For our example program we will compute "Correlation Dimension".  This implementation is trivially O(N^2) steps and easily subdivided in various ways.  Fractional dimension measure is sort of a fun idea.  It is a cute stepping stone into more formal non linear time series analysis. Something we can understand and code in a few minutes...

The high level steps to the algorithm:

Given a set of N points in R^d,

1) Compute all pairs of distances between two points n_i, n_j.  We will use Euclidean distance for dimension d (l2 norm).

2) For a range of epsilons, count how many distances are less than epsilon.

3) Fit a line log-log scatter plot of counts vs epsilon. The slope of this line is the estimate of Correlation Dimension.

Each implementation is located in it's own directory.  The serial version was prototyped in python.  The algorithm is split into two programs.  One which computes steps 1 and 2, then outputs a simple text file. The second program plots and line fits for us.

Steps 1 and 2 were translated to C and extended to various implementations.

