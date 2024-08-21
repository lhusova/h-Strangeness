
#include "Plotter.h"

void PlotSinglePartQA(){

  TFile *file = new TFile("../data/AnalysisResultsLongTrain.root", "READ");


  TH3F *fHist3d = (TH3F *) file->Get("correlate-strangeness/hK0ShortEtaVsPtVsPhi");
  fHist3d->Sumw2();
  TH3F *fHist3dBg = (TH3F *) file->Get("correlate-strangeness/hK0ShortEtaVsPtVsPhiBg");
  fHist3dBg->Sumw2();
  fHist3d->Add(fHist3dBg,-1);

  TCanvas * canPt = Plotter::CreateCanvas("canPt");
  TH1F * ptProj = (TH1F *)fHist3d->Project3D("x");
  Plotter::SetHistAxes(ptProj,"#font[52]{p}_{T} (GeV/#font[52]{c})","dN/dp_{T}");
  Plotter::SetHist(ptProj,"",20,kBlack,1.2);
  for (int i = 0; i < ptProj->GetXaxis()->GetNbins(); i++) {
    ptProj->SetBinContent(i+1,ptProj->GetBinContent(i+1)/ptProj->GetBinWidth(i+1));
    ptProj->SetBinError(i+1,ptProj->GetBinError(i+1)/ptProj->GetBinWidth(i+1));
  }
  // ptProj->GetXaxis()->SetRangeUser(2,50);
  gPad->SetLogy();
  ptProj->DrawCopy();

  TCanvas * canEta = Plotter::CreateCanvas("canEta");
  TH1F * etaProj = (TH1F *)fHist3d->Project3D("y");
  Plotter::SetHistAxes(etaProj,"#eta","#");
  Plotter::SetHist(etaProj,"",20,kBlack,1.2);
  // gPad->SetLogy();
  etaProj->DrawCopy();

  TCanvas * canPhi = Plotter::CreateCanvas("canPhi");
  TH1F * phiProj = (TH1F *)fHist3d->Project3D("z");
  Plotter::SetHistAxes(phiProj,"#varphi","#");
  Plotter::SetHist(phiProj,"",20,kBlack,1.2);
  // gPad->SetLogy();
  phiProj->DrawCopy();
}
