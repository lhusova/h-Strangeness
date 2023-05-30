#include "Plotter.h"

void Plot1DProjections(){

  TFile * fFile = new TFile("../data/Yields_K0s.root");

  TH1D * fullProj = (TH1D *)fFile->Get("phiProj_withBckg_pT2");
  TH1D * UEProj = (TH1D *)fFile->Get("underlying_ev_pT2");
  TH1D * noBckgProj = (TH1D *)fFile->Get("phiProj_noBckg_pT2");

  TCanvas *can = Plotter::CreateCanvas("c");
  fullProj->SetLineColor(kBlack);
  fullProj->GetYaxis()->SetMaxDigits(2);
  fullProj->SetTitle("");
  fullProj->DrawCopy();
  TH1D * fullProjNoErr = (TH1D *) fullProj->Clone();
  fullProjNoErr->SetName("fullProjNoErr");
  for(Int_t i=0;i<72;i++){
    fullProjNoErr->SetBinError(i+1,0);
  }
  TH1D * fullProjNear = (TH1D *) fullProjNoErr->Clone();
  fullProjNear->SetName("fullProjNear");
  fullProjNear->GetXaxis()->SetRangeUser(-0.9,0.9);
  fullProjNear->SetFillColor(kRed+1);
  fullProjNear->SetLineColor(kBlack);
  fullProjNear->Draw("same");
  //
  TH1D * fullProjAway = (TH1D *) fullProjNoErr->Clone();
  fullProjAway->SetName("fullProjNear");
  fullProjAway->GetXaxis()->SetRangeUser(TMath::Pi()-1.4,TMath::Pi()+1.4);
  fullProjAway->SetFillColor(kBlue+1);
  fullProjAway->SetLineColor(kBlack);
  fullProjAway->DrawCopy("same");
  //
  TH1D * UEProjNoErr = (TH1D *) UEProj->Clone();
  UEProjNoErr->SetName("UEProjNoErr");
  for(Int_t i=0;i<72;i++){
    UEProjNoErr->SetBinError(i+1,0);
  }
  UEProjNoErr->SetFillColor(kGreen+1);
  UEProj->SetLineColor(kBlack);
  UEProj->DrawCopy("same");
  UEProjNoErr->SetLineColor(kBlack);
  UEProjNoErr->DrawCopy("same");

}
