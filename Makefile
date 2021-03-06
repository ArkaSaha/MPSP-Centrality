
CC = g++
CCFLAGS_BAD_PRACTICE = -Wno-unused-parameter -Wno-unused-variable
CCFLAGS = -Wall -Wextra -pedantic -pthread -std=c++17 -O3 $(CCFLAGS_BAD_PRACTICE)

mpsp: mpsp.cpp topk.cpp
	$(CC) $(CCFLAGS) mpsp.cpp topk.cpp -o mpsp

wise: wise_experiment.cpp topk.cpp
	$(CC) $(CCFLAGS) wise_experiment.cpp topk.cpp -o wise

dasfaa: dasfaa.cpp
	$(CC) $(CCFLAGS) dasfaa.cpp -o dasfaa

single: single.cpp
	$(CC) $(CCFLAGS) single.cpp -o single

all: wise mpsp dasfaa single

clean:
	rm -f mpsp dasfaa wise single
