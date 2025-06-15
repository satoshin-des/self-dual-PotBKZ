#CC	= g++
CC = mpic++
SRC	= src/SelfDualPotBKZ.cpp
CFLAGS	= -shared -fPIC -Ofast -fopenmp -mtune=native -march=native -mfpmath=both -o libSDPotBKZ.so
LDFLAGS = -lntl

all:
	${CC} ${CFLAGS} ${SRC} ${LDFLAGS}
