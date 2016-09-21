#!/bin/bash

GREEN='\033[0;32m'
NC='\033[0m'

make

# clean old experiment results
rm -f res/speed_vs_accuracy/*
rm -f res/flops_given_accuracy/*

#for i in $( find 'data' -name '*csv' ); do
for i in $( ls data | grep .csv ); do
  echo -e ${GREEN}Run evaluation for $i${NC}
  ./evaluate.o speed_vs_accuracy $i
  ./evaluate.o flops_given_accuracy $i
done
