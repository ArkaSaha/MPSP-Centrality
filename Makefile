
CC = g++
CCFLAGS_BAD_PRACTICE = -Wno-unused-parameter -Wno-unused-variable
CCFLAGS = -Wall -Wextra -pedantic -std=c++17 -O3 $(CCFLAGS_BAD_PRACTICE)

wise: wise_experiment.cpp topk.cpp
	$(CC) $(CCFLAGS)  wise_experiment.cpp topk.cpp -o wise

mpsp: mpsp.cpp topk.cpp
	$(CC) $(CCFLAGS) mpsp.cpp topk.cpp -o mpsp

mpsp_baseline: mpsp_baseline.cpp topk.cpp
	$(CC) $(CCFLAGS) mpsp_baseline.cpp topk.cpp -o mpsp_baseline

all: wise mpsp mpsp_baseline

clean:
	rm -f mpsp mpsp_baseline wise
