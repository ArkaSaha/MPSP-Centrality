#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "You must enter exactly 2 command line arguments: first BA or ER and second is the number of seconds";
    exit 1;
fi

type=$1
seconds=$2

nrnodes=(10000 20000 50000 100000 500000 1000000 5000000 10000000)
for nodes in ${nrnodes[@]}; do
	let edges=$nodes*2;
	if [ "$1" == "BA" ]; then
		let edges=$edges-3;
	fi
	./experiments $1_${nodes}_${edges} ${seconds} &
done

