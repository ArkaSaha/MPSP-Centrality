

CC = g++ 
CCFLAGS = -std=gnu++17 -O2

topk: topk.cpp
	${CC} ${CCFLAGS} topk.cpp  -o topk

mpsp: mpsp.cpp
	${CC} ${CCFLAGS} mpsp.cpp  -o mpsp

all: topk mpsp

