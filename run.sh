#!/bin/bash

GREEN='\033[0;32m'
NC='\033[0m'

make

#for i in $( find 'data' -name '*csv' ); do
for i in $( ls data | grep .csv ); do
  echo -e ${GREEN}Run evaluation for $i${NC}
  ./evaluate.o speed_vs_accuracy $i
  ./evaluate.o flops_given_accuracy $i
done
