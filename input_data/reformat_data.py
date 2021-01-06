#!/usr/bin/env python


# We'll take column of data we want and translate that into a static C array,
#   which we will save in a header file.

with open('monthly-sunspots.csv', 'r') as fh:
    lines = fh.readlines()

# skip header line
lines = lines[1:]
n = len(lines)

with open('monthly_sunspots_data.h', 'w') as fh:

    fh.write('const int monthly_sunspots_n = %d;\n\n' % n)
    fh.write('double monthly_sunspots[%d] = {\n' % n)
    
    for i, line in enumerate(lines):
        s = '\t' + line.strip().split(',')[1]
        if i<n-1:
            s += ',\n'

        fh.write(s)

    fh.write('\n};\n')

