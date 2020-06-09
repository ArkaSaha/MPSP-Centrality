
CC = g++ 
CCFLAGS_BAD_PRACTICE = -Wno-unused-parameter -Wno-unused-variable
CCFLAGS = -Wall -Wextra -pedantic -std=c++17 -O3 $(CCFLAGS_BAD_PRACTICE)

experiments: experiments.cpp topk.cpp 
	$(CC) $(CCFLAGS)  experiments.cpp topk.cpp -o experiments

test: test.cpp topk.cpp 
	$(CC) $(CCFLAGS) test.cpp topk.cpp -o test

mpsp: mpsp.cpp topk.cpp
	$(CC) $(CCFLAGS) mpsp.cpp topk.cpp -o mpsp

mpsp_baseline: mpsp_baseline.cpp topk.cpp
	$(CC) $(CCFLAGS) mpsp_baseline.cpp topk.cpp -o mpsp_baseline

all: experiments mpsp mpsp_baseline test

clean:
	rm -f mpsp mpsp_baseline experiments test
	
