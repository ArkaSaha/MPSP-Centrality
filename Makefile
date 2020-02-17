
CC = g++ 
CCFLAGS = -Wall -Wextra -pedantic -std=gnu++17 -O2
CCFLAGS_BAD_PRACTICE = -Wno-unused-parameter

experiments: experiments.cpp topk.cpp topk.h
	$(CC) $(CCFLAGS) $(CCFLAGS_BAD_PRACTICE) experiments.cpp topk.cpp -o experiments

topk: topk.cpp
	$(CC) $(CCFLAGS) topk.cpp  -o topk

mpsp: mpsp.cpp
	$(CC) $(CCFLAGS) mpsp.cpp  -o mpsp

all: experiments topk mpsp

clean:
	rm topk mpsp experiments
	

