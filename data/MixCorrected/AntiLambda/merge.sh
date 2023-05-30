#!/bin/bash

if [[ -f "./MixCorrected_Signal_AntiLambda_minBias.root" ]]
then
  rm MixCorrected_Signal_AntiLambda_minBias.root
fi
hadd MixCorrected_Signal_AntiLambda_minBias.root MixCorrected_Signal_AntiLambda_*
if [[ -f "./MixCorrected_RightBg_AntiLambda_minBias.root" ]]
then
  rm MixCorrected_RightBg_AntiLambda_minBias.root
fi
hadd MixCorrected_RightBg_AntiLambda_minBias.root MixCorrected_RightBg_AntiLambda_*
if [[ -f "./MixCorrected_LeftBg_AntiLambda_minBias.root" ]]
then
  rm MixCorrected_LeftBg_AntiLambda_minBias.root
fi
hadd MixCorrected_LeftBg_AntiLambda_minBias.root MixCorrected_LeftBg_AntiLambda_*
