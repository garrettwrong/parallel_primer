CC=mpicc
MPI_FLAGS=
CFLAGS=-O2 -Wall $(MPI_FLAGS) -lm -I../input_data

correlation_dimension_mpi.x: correlation_dimension_mpi.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f correlation_dimension_mpi.x

.PHONY: clean
