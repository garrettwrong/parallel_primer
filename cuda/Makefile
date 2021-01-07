CC=nvcc
CFLAGS=-O2 -lm -I../input_data --compiler-options -Wall

correlation_dimension_cuda.x: correlation_dimension_cuda.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f correlation_dimension_cuda.x

.PHONY: clean
