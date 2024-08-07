
#include "Plotter.h"
void spectrum(){

  TFile *file = new TFile("../data/AnalysisResults_ForMC_minBiasRuns.root", "READ");
  TH3F * genK0 = (TH3F*) file->Get("correlate-strangeness/ClosureTest/hK0Short");
  genK0->Sumw2();
  TH1F* genK0PojPt = (TH1F*) genK0->ProjectionX();
  Plotter::SetHist(genK0PojPt,"",20,kRed+1,1.);


  TH3F * recK0 = (TH3F*) file->Get("correlate-strangeness/hK0ShortEtaVsPtVsPhi");
  recK0->Sumw2();
  TH3F * recK0Bg = (TH3F*) file->Get("correlate-strangeness/hK0ShortEtaVsPtVsPhiBg");
  recK0Bg->Sumw2();
  recK0->Add(recK0Bg,-1);
  TH1F* recK0PojPt = (TH1F*) recK0->ProjectionX();
  Plotter::SetHist(recK0PojPt,"",24,kBlue+1,1.);
  Plotter::SetHistAxes(genK0PojPt,"p_{T}","dN/dp_{T}");
  for (size_t i = 0; i < recK0PojPt->GetXaxis()->GetNbins(); i++) {
    recK0PojPt->SetBinContent(i+1,recK0PojPt->GetBinContent(i+1)/recK0PojPt->GetBinWidth(i+1));
    recK0PojPt->SetBinError(i+1,recK0PojPt->GetBinError(i+1)/recK0PojPt->GetBinWidth(i+1));
    genK0PojPt->SetBinContent(i+1,genK0PojPt->GetBinContent(i+1)/genK0PojPt->GetBinWidth(i+1));
    genK0PojPt->SetBinError(i+1,genK0PojPt->GetBinError(i+1)/genK0PojPt->GetBinWidth(i+1));
  }

  TCanvas *c = Plotter::CreateCanvas("c",700,750,true);
  TPad * padRatio = new TPad(Form("padRatio"),"",0.001,0.001,0.999,0.3);
  padRatio->SetMargin(0.12,0.02,0.25,0.01);
  padRatio->Draw();

  c->GetPadSave()->cd();
  genK0PojPt->GetXaxis()->SetRangeUser(0,20);
  genK0PojPt->DrawCopy();
  recK0PojPt->DrawCopy("same");
  padRatio->cd();
  TH1F* recK0PojPtCopy = (TH1F*) recK0PojPt->Clone();
  recK0PojPtCopy->SetName("recK0PojPtCopy");
  recK0PojPtCopy->Divide(genK0PojPt);
  Plotter::SetHistAxesSmallPad(recK0PojPtCopy,"p_{T}","rec/gen");
  recK0PojPtCopy->GetXaxis()->SetRangeUser(0,20);
  recK0PojPtCopy->DrawCopy();
}
