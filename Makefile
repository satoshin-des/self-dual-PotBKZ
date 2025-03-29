CC	= g++
SRC	= src/SelfDualPotBKZ.cpp
CFLAGS	= -shared -fPIC -O3 -fopenmp -mtune=native -march=native -mfpmath=both -o libSDPotBKZ.so
LDFLAGS = -lntl

all:
	${CC} ${CFLAGS} ${SRC} ${LDFLAGS}
