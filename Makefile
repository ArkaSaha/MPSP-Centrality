
CC = g++
CCFLAGS_BAD_PRACTICE = -Wno-unused-parameter -Wno-unused-variable
CCFLAGS = -Wall -Wextra -pedantic -pthread -std=c++17 -O3 $(CCFLAGS_BAD_PRACTICE)

wise: wise_experiment.cpp topk.cpp
	$(CC) $(CCFLAGS) wise_experiment.cpp topk.cpp -o wise

mpsp: mpsp.cpp topk.cpp
	$(CC) $(CCFLAGS) mpsp.cpp topk.cpp -o mpsp

dasfaa: dasfaa.cpp
	$(CC) $(CCFLAGS) dasfaa.cpp -o dasfaa

all: wise mpsp dasfaa

clean:
	rm -f mpsp dasfaa wise
