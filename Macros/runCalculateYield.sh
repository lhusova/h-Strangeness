#!/bin/bash

for i in {0..1};
do
  for j in {0..10};
  do
      root.exe "CalculateYield.C($i,$j)" -q -b

  done
done
