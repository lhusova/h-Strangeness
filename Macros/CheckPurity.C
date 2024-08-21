#include "Plotter.h"

void CheckPurity(){

  TFile * fFile = new TFile(Form("../data/AnalysisResults_ForMCclosure_5sigma.root"));
  TH3F * recoAllTriggers  = (TH3F *) fFile->Get("correlate-strangeness/hTriggerAllSelectedEtaVsPt");
  recoAllTriggers->Sumw2();
  TH3F * recoPrimary  = (TH3F *) fFile->Get("correlate-strangeness/hTriggerPrimaryEtaVsPt");
  recoPrimary->Sumw2();

  recoPrimary->GetYaxis()->SetRangeUser(-0.8,0.8);
  recoAllTriggers->GetYaxis()->SetRangeUser(-0.8,0.8);
  TH1F * recoPrimaryMB_pt = (TH1F *)recoPrimary->Project3D("y");
  recoPrimaryMB_pt->SetName("recoPrimaryMB_pt");
  TH1F * recoAllMB_pt = (TH1F *)recoAllTriggers->Project3D("y");
  recoAllMB_pt->SetName("recoAllMB_pt");

  TCanvas *c = Plotter::CreateCanvas("c");
  recoPrimaryMB_pt->Divide(recoAllMB_pt);
  Plotter::SetHistAxes(recoPrimaryMB_pt,"#font[52]{p}_{T} (GeV/#font[52]{c})","purity");
  Plotter::SetHist(recoPrimaryMB_pt,"",20,kBlack,1.2);
  recoPrimaryMB_pt->GetXaxis()->SetRangeUser(1,50);
  recoPrimaryMB_pt->DrawCopy();

  TH1F * recoPrimary_pt[8];
  TH1F * recoAll_pt[8];
  Color_t col[] = {kRed+2,kOrange+7,kOrange-2,kYellow+1,kSpring+9,kTeal-5,kAzure+10, kBlue+2};
  TString multiplicityPave[]={"0-1%","1-10%","10-20%","20-30%","30-40%","40-50%","50-70%","70-100%"};
  TLegend *leg = Plotter::CreateLegend(0.15, 0.45, 0.15, 0.45,0.05);
  leg->AddEntry(recoPrimaryMB_pt,"MB","pl");

  for (Int_t i = 0; i < 8; i++) {
    recoPrimary->GetZaxis()->SetRange(i+1,i+1);
    recoPrimary_pt[i] = (TH1F *)recoPrimary->Project3D("y");
    recoPrimary_pt[i]->SetName(Form("recoPrimary%i_pt",i));
    recoAllTriggers->GetZaxis()->SetRange(i+1,i+1);
    recoAll_pt[i] = (TH1F *)recoAllTriggers->Project3D("y");
    recoAll_pt[i]->SetName(Form("recoAll%i_pt",i));
    Plotter::SetHist(recoPrimary_pt[i],"",20,col[i],1.2);
    recoPrimary_pt[i]->Divide(recoAll_pt[i]);
    recoPrimary_pt[i]->DrawCopy("same");
    leg->AddEntry(recoPrimary_pt[i],multiplicityPave[i].Data(),"pl");
  }
  recoPrimaryMB_pt->DrawCopy("same");
  leg->Draw();

  TPaveText * pave = new TPaveText();
  Plotter::SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.97, 0.7,0.95);
  pave->AddText("ALICE, Work in Progress");
  pave->AddText("pp, 13.6 TeV");
  // pave->AddText(" |#eta| < 0.8");
  pave->AddText(Form("2 < #font[52]{p}_{T} < 50 GeV/#font[52]{c}"));
  pave->Draw("same");

}
