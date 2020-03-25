
CC = g++ 
CCFLAGS = -Wall -Wextra -pedantic -std=c++17 -O3
CCFLAGS_BAD_PRACTICE = -Wno-unused-parameter

experiments: experiments.cpp topk.cpp 
	$(CC) $(CCFLAGS) $(CCFLAGS_BAD_PRACTICE) experiments.cpp topk.cpp -o experiments -lstdc++fs

test: test.cpp topk.cpp 
	$(CC) $(CCFLAGS) test.cpp topk.cpp -o test

topk: topk.cpp 
	$(CC) $(CCFLAGS) topk.cpp  -o topk

mpsp: mpsp.cpp topk.cpp
	$(CC) $(CCFLAGS) mpsp.cpp  topk.cpp -o mpsp

all: experiments topk mpsp

clean:
	rm -f topk mpsp experiments test
	
