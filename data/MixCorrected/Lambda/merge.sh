#!/bin/bash

if [[ -f "./MixCorrected_Signal_Lambda_minBias.root" ]]
then
  rm MixCorrected_Signal_Lambda_minBias.root
fi
hadd MixCorrected_Signal_Lambda_minBias.root MixCorrected_Signal_Lambda_*
if [[ -f "./MixCorrected_RightBg_Lambda_minBias.root" ]]
then
  rm MixCorrected_RightBg_Lambda_minBias.root
fi
hadd MixCorrected_RightBg_Lambda_minBias.root MixCorrected_RightBg_Lambda_*
if [[ -f "./MixCorrected_LeftBg_Lambda_minBias.root" ]]
then
  rm MixCorrected_LeftBg_Lambda_minBias.root
fi
hadd MixCorrected_LeftBg_Lambda_minBias.root MixCorrected_LeftBg_Lambda_*
