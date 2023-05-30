#!/bin/bash

if [[ -f "./MixCorrected_Signal_OmegaMinus_minBias.root" ]]
then
  rm MixCorrected_Signal_OmegaMinus_minBias.root
fi
hadd MixCorrected_Signal_OmegaMinus_minBias.root MixCorrected_Signal_OmegaMinus_*
if [[ -f "./MixCorrected_RightBg_OmegaMinus_minBias.root" ]]
then
  rm MixCorrected_RightBg_OmegaMinus_minBias.root
fi
hadd MixCorrected_RightBg_OmegaMinus_minBias.root MixCorrected_RightBg_OmegaMinus_*
if [[ -f "./MixCorrected_LeftBg_OmegaMinus_minBias.root" ]]
then
  rm MixCorrected_LeftBg_OmegaMinus_minBias.root
fi
hadd MixCorrected_LeftBg_OmegaMinus_minBias.root MixCorrected_LeftBg_OmegaMinus_*
