#!/bin/bash

for i in {0..0};
do
  for j in {0..8};
  do
    for k in {0..3};
    do
        root.exe "CalculateYield.C($i,$k,$j)" -q -b
        # root.exe "CompareYieldSpectra.C($i,$j)" -q -b
    done
  done
done
