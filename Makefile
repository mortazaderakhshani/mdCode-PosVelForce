

CXX = g++
CC = gcc
LAPACK =/Users/mortaza/lib 
OPTS = -O3 -ftree-vectorize

mc: code2.c stringlib.c stringlib.h  
	$(CC) -c  code2.c stringlib.c $(OPTS) 
	$(CC) code2.o stringlib.o $(OPTS) -o code2.x






