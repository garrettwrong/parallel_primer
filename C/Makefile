CC=gcc
CFLAGS=-std=c99 -O2 -Wall -lm -I../input_data

correlation_dimension.x: correlation_dimension.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -f correlation_dimension.x

.PHONY: clean
