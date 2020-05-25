#!/bin/bash


ml matlab
tau="$1"
seed="$2"

echo "tau${tau}_seed${seed}"
matlab -r "run_nsfa_mydata($tau, $seed);exit(0)"
