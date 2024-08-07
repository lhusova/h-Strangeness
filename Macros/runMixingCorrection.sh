#!/bin/bash

array_AssocParticles=("K0Short")
# "Lambda" "AntiLambda")
# array_AssocParticles=("Pion")
# array_AssocParticles=( "XiMinus" "XiPlus" "OmegaMinus" "OmegaPlus")
array_InvMassRegions=("LeftBg" "Signal" "RightBg")
# array_InvMassRegions=("p" )
for i in ${!array_AssocParticles[@]};
do
  for j in ${!array_InvMassRegions[@]};
  do
    for imult in {0..7}
    do
      root.exe "MixingCorrection.C(\"${array_AssocParticles[$i]}\",\"${array_InvMassRegions[$j]}\",$imult)" -q -b
    done
  done
done
