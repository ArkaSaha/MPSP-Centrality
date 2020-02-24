
CC = g++ 
CCFLAGS = -Wall -Wextra -pedantic -std=gnu++17 -O2
CCFLAGS_BAD_PRACTICE = -Wno-unused-parameter

experiments: experiments.cpp topk.cpp 
	$(CC) $(CCFLAGS) $(CCFLAGS_BAD_PRACTICE) experiments.cpp topk.cpp statistics.h -o experiments

topk: topk.cpp topk.h
	$(CC) $(CCFLAGS) statistics.h topk.cpp  -o topk

mpsp: mpsp.cpp
	$(CC) $(CCFLAGS) mpsp.cpp  -o mpsp

all: experiments topk mpsp

clean:
	rm topk mpsp experiments
	

