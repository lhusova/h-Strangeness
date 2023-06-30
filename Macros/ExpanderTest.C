#include "Expander.C"

void ExpanderTest(){ 
  // Macro to test expansion of axes
  
  // Get the full histogram in the default way
  TFile *file = new TFile("AnalysisResultsFull.root", "READ");
  TString histoName = "correlate-strangeness/sameEvent/Signal/K0Short";
  THnF *hNdK0ShortSignalFull = (THnF*) file->Get(histoName.Data());
  
  // Get via expander
  TFile *file2 = new TFile("AnalysisResultsCompact.root", "READ");
  THnF *hNdK0ShortSignalCompact = GetTHnF(file2, "K0ShortSignal", histoName.Data());
  
  TH1D *hProjFull = hNdK0ShortSignalFull->Projection(4);
  hProjFull->SetName("hProjFull");
  TH1D *hProjComp = hNdK0ShortSignalCompact->Projection(4);
  hProjComp->SetName("hProjComp");
  
  hProjFull->SetLineWidth(5);
  hProjFull->SetLineStyle(2);
  hProjFull->SetLineColor(kBlack);
  hProjFull->Draw();
  hProjComp->SetLineWidth(1);
  hProjComp->SetLineColor(kRed);
  hProjComp->Draw("same");
} 


