#include "Plotter.h"

void PrepareEfficiencyPlots(Int_t iPart=0, TString suffix = ""){

  TString particleNames[]={"K0Short"};
  TString partSym[]={"K_{S}^{0}"};
  TFile * fFile = new TFile("../data/EffHist.root");

  TH2F* h2d = (TH2F*) fFile->Get(Form("%s%s",particleNames[iPart].Data(),suffix.Data()));
  TGaxis::SetMaxDigits(4);
  TCanvas * can2d = Plotter::CreateCanvas("can2d");
  gPad->SetTheta(45);
  gPad->SetPhi(40);
  gPad->GetFrame()->SetLineColor(0);
  Plotter::Set2DHistAxes(h2d,"#font[52]{p}_{T} (GeV/#font[52]{c})","#eta","#varepsilon","");
  h2d->DrawCopy("surf2 fb");
  h2d->DrawCopy("surf same fb");

  TPaveText * paveGen = new TPaveText();
  Plotter::SetPaveText(paveGen,42,0.05, 0, 0,12,0,0.02,0.45, 0.8,0.97);
  paveGen->AddText(Form("ALICE, Work in progress"));
  paveGen->AddText("pp, 13.6 TeV");
  paveGen->AddText(Form("%s",partSym[iPart].Data()));
  paveGen->Draw();

  TH1F * h1dPt = (TH1F *)fFile->Get(Form("pt_%s%s",particleNames[iPart].Data(),suffix.Data()));
  TCanvas * canPt = Plotter::CreateCanvas("canPt");
  Plotter::SetHist(h1dPt,"",24,kBlue+1,1.);
  Plotter::SetHistAxes(h1dPt,"#font[52]{p}_{T} (GeV/#font[52]{c})","#varepsilon");
  h1dPt->DrawCopy();
  paveGen->Draw();

  TH1F * h1dEta = (TH1F *)fFile->Get(Form("eta_%s%s",particleNames[iPart].Data(),suffix.Data()));
  TCanvas * canEta = Plotter::CreateCanvas("canEta");
  Plotter::SetHist(h1dEta,"",24,kBlue+1,1.);
  Plotter::SetHistAxes(h1dEta,"#eta","#varepsilon");
  h1dEta->DrawCopy();
  paveGen->Draw();



}
