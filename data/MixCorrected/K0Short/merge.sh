#!/bin/bash

if [[ -f "./MixCorrected_Signal_K0Short_minBias.root" ]]
then
  rm MixCorrected_Signal_K0Short_minBias.root
fi
hadd MixCorrected_Signal_K0Short_minBias.root MixCorrected_Signal_K0Short_*
if [[ -f "./MixCorrected_RightBg_K0Short_minBias.root" ]]
then
  rm MixCorrected_RightBg_K0Short_minBias.root
fi
hadd MixCorrected_RightBg_K0Short_minBias.root MixCorrected_RightBg_K0Short_*
if [[ -f "./MixCorrected_LeftBg_K0Short_minBias.root" ]]
then
  rm MixCorrected_LeftBg_K0Short_minBias.root
fi
hadd MixCorrected_LeftBg_K0Short_minBias.root MixCorrected_LeftBg_K0Short_*
