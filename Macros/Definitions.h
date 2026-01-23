#include <TString.h>
#include "TROOT.h"

const Int_t nPtTriggBins = 4;
const Int_t nMultBins = 9;
const Int_t nRegions = 3;
const Int_t nPtBins = 10;

TString nameSave[]={"K0s","Lam","Xi","Omega","Pion","Hadron"};
TString finalNames[]={"K_{S}^{0}","(#Lambda+#bar{#Lambda})","(#Xi^{+}+#Xi^{-})","(#Omega^{+}+#Omega^{-})","#pi^{+}+#pi^{-}","h^{#pm}"};
TString multiplicityNames[]={"MB","0_1Mult","1_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_70Mult","70_100Mult"};
// TString multiplicityNames[]={"MB","0_1Mult","1_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_70Mult","70_90Mult","70_100Mult"};
TString namesRegions[] = {"fHistNear", "fHistAway", "fHistUE"};
TString paveRegions[] = {"Near-side, |#Delta#varphi|<#pi/2", "Away-side, |#Delta#varphi-#pi|<#pi/2", "Underlying event"};
TString name[]={"K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus","Pion","Hadron"};
TString finalNamesSeparate[]={"K_{S}^{0}","#Lambda","#bar{#Lambda}","#Xi^{-}","#bar{#Xi^{+}}","#Omega^{-}","#bar{#Omega^{+}}","#pi^{+}+#pi^{-}","h"};
TString multiplicityPave[]={"MB","0#minus1%","1#minus10%","10#minus20%","20#minus30%","30#minus40%","40#minus50%","50#minus70%","70#minus100%"};
// TString multiplicityPave[]={"MB","0#minus1%","1#minus10%","10#minus20%","20#minus30%","30#minus40%","40#minus50%","50#minus70%","70#minus90%"};
TString nameBackground[] = {"flat","long_range","flow_modulation_v2","flow_modulation_v3","Signal"};
TString InvMassRanges[] = {"Signal", "LeftBg", "RightBg"};
TString namesRegionsShort[nRegions] = {"Near", "Away", "UE"};
TString multLabels[] = {"70#minus100%", "50#minus70%", "40#minus50%", "30#minus40%", "20#minus30%", "10#minus20%", "1#minus10%", "0#minus1%", "0#minus100%"};
// TString multLabels[] = {"70#minus90%", "50#minus70%", "40#minus50%", "30#minus40%", "20#minus30%", "10#minus20%", "1#minus10%", "0#minus1%", "0#minus90%"};
TString PhiRegions[nRegions] = {"|#Delta#it{#varphi}| < #pi/2", "#pi/2 < |#Delta#it{#varphi}| < 3/2#pi", "-#pi/2 < |#Delta#it{#varphi}| < 3/2#pi"};
TString triggHistName[]={"V0","Cascade","Pion","Hadron"};

TString nameEff[]={"Trigger","K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus","Pion"};
TString finalNamesEff[]={"charged hadrons","K_{S}^{0}","#Lambda","#bar{#Lambda}","#Xi^{-}","#Xi^{+}","#Omega^{-}","#Omega^{+}","#pi^{+}+#pi^{-}"};


Float_t particleMass[] = {0.497, 1.115, 1.321, 1.672, 0.1396, 0.1396};

Double_t ptTriggBins[]={2.,4.,6.,10.,50.};
Double_t ptBins[]={0.2,0.5,1.,2.,3.,4.,6.,8.,10.,12,15};
Double_t ptBinsHadron[]={0.2,0.6,0.8,1,1.5,2,2.5,3,3.5,4,5,6,7,8,9.5,11,13,15};
// Double_t ptBinsHadronRebin[]={0.2,0.6,0.8,1,1.5,2,2.5,3,3.5,4,5,6,8,11,15};
Int_t binRebin[]={1,2,3,4,5,6,7,8,9,10,11,12,14,18};
// Double_t ptBins_err[nPtBins]={0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};

Color_t colorRegion[nRegions]={kRed+1, kBlue+1, kGreen+2};
Color_t colRegions[nRegions][5] = {{kRed - 7, kRed - 4, kRed + 1, kRed + 2, kRed + 4}, {kAzure + 6, kAzure + 7, kBlue, kBlue + 1, kBlue + 3}, {kSpring + 7, kGreen + 1, kGreen + 2, kGreen + 3, kGreen + 5}};
Int_t markerTrigg[nPtTriggBins]={20,21,34,47};

// Int_t colorsMultiplicity[nMultBins];
// SetColorPalete(nMultBins,colorsMultiplicity);