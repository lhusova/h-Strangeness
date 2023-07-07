#include "Expander.C"

void ExpandTHns(){

  TFile *file = new TFile("../data/AnalysisResults_Hyperloop_24_06_Pion.root", "READ");

  // TString region[]={"Signal","RightBg","LeftBg"};
  TString region[]={"p"};
  // TString particle[]={"K0Short","Lambda","AntiLambda"};
  // TString particle[]={"XiMinus","XiPlus","OmegaMinus","OmegaPlus"};
  TString particle[]={"Pion"};

  Int_t nPart = sizeof(particle) / sizeof(TString);
  Int_t nReg = sizeof(region) / sizeof(TString);
  THnF *hNSameFull[3][3];
  TString histoName;

  TFile *fileNew = new TFile("../data/AnalysisResults_Hyperloop_24_06_Pion_Expanded.root", "RECREATE");

  for (size_t j = 0; j < nPart; j++) {
    for (size_t i = 0; i < nReg; i++) {
      histoName = "correlate-strangeness/sameEvent/";
      if(nReg>1){
        histoName+=region[i];
        histoName+="/";
      }
      histoName+=particle[j];
      hNSameFull[i][j] = GetTHnF(file,Form("Same_%s_%s_Expanded",particle[j].Data(),region[i].Data()),histoName.Data());
      hNSameFull[i][j]->Write();
    }
  }

  for (size_t j = 0; j < nPart; j++) {
    for (size_t i = 0; i < nReg; i++) {
      histoName = "correlate-strangeness/mixedEvent/";
      if(nReg>1){
        histoName+=region[i];
        histoName+="/";
      }
      histoName+=particle[j];
      hNSameFull[i][j] = GetTHnF(file,Form("Mixed_%s_%s_Expanded",particle[j].Data(),region[i].Data()),histoName.Data());
      hNSameFull[i][j]->Write();
    }
  }


}
