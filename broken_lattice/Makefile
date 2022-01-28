CFLAGS=-Wall -Wextra -Wpedantic -g -std=gnu99 -lfftw3_threads -lfftw3 -lm -lpthread -lgsl -lgslcblas
.DEFAULT:= all
.PHONY: all

OBJECTS=main.o util.o field_gen.o power_spectrum.o 
all: main

main: $(OBJECTS)
	 cc $(CFLAGS) $(OBJECTS) -o main 

clean:
	rm *.o main
