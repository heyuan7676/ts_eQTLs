#!/bin/bash


ml matlab
for tau in 1000 500 100 20 10 5 1
do 
	for seed in {1..10}
	do
		echo "tau${tau}_seed${seed}"
		matlab -r "run_nsfa_mydata($tau, $seed);exit(0)"
	done
done
