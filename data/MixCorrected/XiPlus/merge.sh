#!/bin/bash

if [[ -f "./MixCorrected_Signal_XiPlus_minBias.root" ]]
then
  rm MixCorrected_Signal_XiPlus_minBias.root
fi
hadd MixCorrected_Signal_XiPlus_minBias.root MixCorrected_Signal_XiPlus_*
if [[ -f "./MixCorrected_RightBg_XiPlus_minBias.root" ]]
then
  rm MixCorrected_RightBg_XiPlus_minBias.root
fi
hadd MixCorrected_RightBg_XiPlus_minBias.root MixCorrected_RightBg_XiPlus_*
if [[ -f "./MixCorrected_LeftBg_XiPlus_minBias.root" ]]
then
  rm MixCorrected_LeftBg_XiPlus_minBias.root
fi
hadd MixCorrected_LeftBg_XiPlus_minBias.root MixCorrected_LeftBg_XiPlus_*
