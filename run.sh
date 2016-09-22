#!/bin/bash

GREEN='\033[0;32m'
NC='\033[0m'

make

# clean old experiment results
rm -f res/speed-vs-accuracy/*
rm -f res/flops-given-accuracy/*

#for i in $( find 'data' -name '*csv' ); do
for i in $( ls data | grep .csv ); do
  echo -e ${GREEN}Run evaluation for $i${NC}
  ./evaluate.o speed-vs-accuracy $i
  ./evaluate.o flops-given-accuracy $i
done
