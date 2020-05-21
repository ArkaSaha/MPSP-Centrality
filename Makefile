
CC = g++ 
CCFLAGS = -Wall -Wextra -pedantic -std=gnu++17 -O3
CCFLAGS_BAD_PRACTICE = -Wno-unused-parameter

experiments: experiments.cpp topk.cpp 
	$(CC) $(CCFLAGS) $(CCFLAGS_BAD_PRACTICE) experiments.cpp topk.cpp -o experiments


test: test.cpp topk.cpp 
	$(CC) $(CCFLAGS) test.cpp topk.cpp -o test

topk: topk.cpp 
	$(CC) $(CCFLAGS) topk.cpp  -o topk

mpsp: mpsp.cpp topk.cpp
	$(CC) $(CCFLAGS) mpsp.cpp  topk.cpp -o mpsp

mpsp_baseline: mpsp_baseline.cpp topk.cpp
	$(CC) $(CCFLAGS) mpsp_baseline.cpp  topk.cpp -o mpsp_baseline

all: experiments topk mpsp mpsp_baseline

clean:
	rm -f topk mpsp mpsp_baseline experiments test
	
