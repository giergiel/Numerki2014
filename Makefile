CC=g++
#delete -DLINUX if system doesnt provide sys/types.h and unistd.h POSIX OS API
FLAGS=-O2 -Wall -std=gnu++1y -DLINUX
objects := $(patsubst %.cpp,%.out,$(wildcard ./*.cpp))

all: CW1.out $(objects)

%.out: %.cpp
	${CC} -o $@ ${FLAGS} $< 

CW1.out: CW1.cpp
	${CC} ${FLAGS} -o CW1.out CW1.cpp


.PHONY: all
