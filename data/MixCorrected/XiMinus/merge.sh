#!/bin/bash

if [[ -f "./MixCorrected_Signal_XiMinus_minBias.root" ]]
then
  rm MixCorrected_Signal_XiMinus_minBias.root
fi
hadd MixCorrected_Signal_XiMinus_minBias.root MixCorrected_Signal_XiMinus_*
if [[ -f "./MixCorrected_RightBg_XiMinus_minBias.root" ]]
then
  rm MixCorrected_RightBg_XiMinus_minBias.root
fi
hadd MixCorrected_RightBg_XiMinus_minBias.root MixCorrected_RightBg_XiMinus_*
if [[ -f "./MixCorrected_LeftBg_XiMinus_minBias.root" ]]
then
  rm MixCorrected_LeftBg_XiMinus_minBias.root
fi
hadd MixCorrected_LeftBg_XiMinus_minBias.root MixCorrected_LeftBg_XiMinus_*
