#!/usr/bin/env python

import numpy as np

# We'll take column of data we want and translate that into a static C array,
#   which we will save in a header file.

n = 10000
dim = 2
dats = np.random.rand(n, dim)

np.save('random_2d_data.npy', dats)

with open('random_2d_data.h', 'w') as fh:

    fh.write('const int input_data_n = %d;\n\n' % n)
    fh.write('const int input_data_dim = %d;\n\n' % dim)
    fh.write('double input_data[%d] = {\n' % (n*dim))
    
    for i, line in enumerate(dats):
        s = '\t' + ', '.join([str(x) for x in line])
        if i<n-1:
            s += ',\n'

        fh.write(s)

    fh.write('\n};\n')



