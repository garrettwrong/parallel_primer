CC=gcc
OMP_FLAGS=-fopenmp
CFLAGS=-O2 -Wall -lm $(OMP_FLAGS) -I../input_data

correlation_dimension_omp.x: correlation_dimension_omp.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f correlation_dimension_omp.x

.PHONY: clean
