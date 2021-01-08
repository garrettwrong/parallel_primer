CC=nvcc
CFLAGS=-O2 -lm -I../input_data --compiler-options -Wall -arch=sm_35

correlation_dimension_cuda.x: correlation_dimension_cuda.cu
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f correlation_dimension_cuda.x

.PHONY: clean
