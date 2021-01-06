#/usr/bin/env python

import os
import numpy as np

import matplotlib.pyplot as plt

def all_pairs_distances(X):
    n = len(X)
    D = np.zeros((n,n), dtype=np.float64)

    # Compute pair wise distances.
    #   Yes, this is intentionally not optimized much... it is an example for non programmers.
    #   For those inclined, you can easily rewrite this as array opts and only compute upper/lower triangle.
    for i, xi in enumerate(X):
        for j, xj in enumerate(X):
            D[i, j] = np.linalg.norm(xi - xj)

    return D

def correlation_integrals(D, epsilons=None, fname='correlation_integrals.dat'):

    if epsilons is None:
        # Lets generate a reasonable epsilon range.
        minD = np.min(D)
        maxD = np.max(D)
        # Just ignore minD, for log(0) reasons
        epsilons = np.linspace(minD, maxD, num=100)[1:]

    C = np.zeros(len(epsilons), dtype=np.int)

    for i, epsilon in enumerate(epsilons):
        C[i] = np.sum(D<epsilon)

    # Optionally Write out to a text file
    if fname is not None:
        with open(fname, 'w') as fh:
            for i, epsilon in enumerate(epsilons):
                # We'll use the same formatting later in C
                fh.write('%f %d\n' % (epsilon, C[i]))

    return epsilons, C


def estimate_dimension(epsilons, C):

    # Fit a line
    X = np.log(epsilons)
    Y = np.log(C)
    poly = np.poly1d(np.polyfit(X, Y, 1))
    slope = poly[1]
    
    # Generate log-log plot, saving as a file
    plt.plot(X, Y)           # Correlation Integrals
    plt.plot(X, poly(X))     # Fit line
    plt.savefig('correlation_dimension.png')
    return slope


def estimate_dimension_from_file(fn):
    with open(fname, 'r') as fh:
        lines = fh.readlines()

    m = len(lines)

    epsilons = np.zeros(m, dtype=np.float64)
    C = np.zeros(m, dtype=np.float64)

    for i, line in enumerate(lines):
        epsilons[i], C[i] = line.strip().split()

    return estimate_dimension(epsilons, C)


def correlation_dimension(X):

    # Step 1
    D = all_pairs_distances(X)

    # Step 2
    epsilons, C = correlation_integrals(D)

    # Step 3
    dim = estimate_dimension(epsilons, C)
    return dim


if __name__ == "__main__":

    # Load some example data
    dirp = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dirp, '../input_data/monthly-sunspots.csv'), 'r') as fh:
        lines = fh.readlines()

    # Take second column, make into an array, skipping header
    X = np.array([l.strip().split(',')[1] for l in lines[1:]], dtype=np.float64)

    print('Estimated Correlation Dimension: %f' % correlation_dimension(X))