
#Define C Compiler

CC=gcc -march=native -O3 -funroll-loops -ftree-vectorize

# Compiler flags:
# -g Add debuging information to the executable file
# -Wall turn on most of compiler warnings
CFLAGS=-Wall -pedantic

# OpenMP Library path
LFLAGS=-L/usr/lib/gcc/x86_64-linux-gnu/12/include

# Define libraries to be used
LIBS= -lm  

#Link openmp to gcc compiler
OMPLINK=-fopenmp 

#Object for compilation
OBJ_OMP=modelling_2D_omp.o  modelling_utils.o

# Name of executable file
MAIN_OMP=modelling_omp

all:  $(MAIN_OMP)
	@echo Compiling executable $(MAIN_OMP)


$(MAIN_OMP): $(OBJ_OMP)
	$(CC) $(OMPLINK) -o $@ $^ $(LIBS) $(LFLAGS)

%.o: %.c 
	$(CC) -c -o $@ $< $(CFLAGS)


modelling_2D_omp.o: modelling_2D_omp.c
	$(CC) $(OMPLINK) -c -o $@ $< $(CFLAGS)


clean:
	rm -f *.o $(MAIN_OMP) 
