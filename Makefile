

CC = g++ 
CCFLAGS = -Wall -Wextra -pedantic -std=gnu++17 -O2

topk: topk.cpp
	${CC} ${CCFLAGS} topk.cpp  -o topk

mpsp: mpsp.cpp
	${CC} ${CCFLAGS} mpsp.cpp  -o mpsp

all: topk mpsp

clean:
	rm topk mpsp
	

