#!/bin/bash

array_AssocParticles=("K0Short" "Lambda" "AntiLambda")
# "XiMinus", "XiPlus", "OmegaMinus", "OmegaPlus")
array_InvMassRegions=("LeftBg" "Signal" "RightBg")
for i in ${!array_AssocParticles[@]};
do
  for j in ${!array_InvMassRegions[@]};
  do
    root.exe "MixingCorrection.C(\"${array_AssocParticles[$i]}\",\"${array_InvMassRegions[$j]}\")" -q -b
  done
done
