#include "Expander.C"

void ExpandTHns(){

  TFile *file = new TFile("../data/AnalysisResults_hhClosure_LHC25b4b_multEff.root", "READ");

  // TString region[]={"Signal","LeftBg","RightBg"};
  TString region[]={"Signal"};
  // TString particle[]={"K0Short"};//,"Lambda","AntiLambda"};
  // TString particle[]={"XiMinus","XiPlus"};//,"OmegaMinus","OmegaPlus"};
  // TString particle[]={"Pion"};
  TString particle[]={"Hadron"};

  Int_t nPart = sizeof(particle) / sizeof(TString);
  Int_t nReg = sizeof(region) / sizeof(TString);
  THnF *hNSameFull[3][3];
  TString histoName;

  TFile *fileNew = new TFile("../data/Expanded/AnalysisResults_hhrec_Mixing_LHC25b4b_multEff.root", "RECREATE");

  for (size_t j = 0; j < nPart; j++) {
    for (size_t i = 0; i < nReg; i++) {
      histoName = "h-strange-correlation/sameEvent/";
      // histoName = "h-strange-correlation__MCgen/ClosureTest/sameEvent/";
      histoName+=region[i];
      histoName+="/";
      histoName+=particle[j];
      cout << histoName << endl;
      hNSameFull[i][j] = GetTHnF(file,Form("Same_%s_%s_Expanded",particle[j].Data(),region[i].Data()),histoName.Data(),false);
      hNSameFull[i][j]->Write();
    }
  }

  for (size_t j = 0; j < nPart; j++) {
    for (size_t i = 0; i < nReg; i++) {
      histoName = "h-strange-correlation/mixedEvent/";
      histoName+=region[i];
      histoName+="/";
      histoName+=particle[j];
      hNSameFull[i][j] = GetTHnF(file,Form("Mixed_%s_%s_Expanded",particle[j].Data(),region[i].Data()),histoName.Data(),false);
      hNSameFull[i][j]->Write();
    }
  }

fileNew->Close();
}
