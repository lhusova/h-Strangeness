#include "Plotter.h"

void CorrectSpectra(Int_t part = 0){

  TString name[]={"K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus","Pion"};
  TString finalNames[]={"K_{S}^{0}","#Lambda","#bar{#Lambda}","#Xi^{-}","#Xi^{+}","#Omega^{-}","#Omega^{+})","#pi^{+}+#pi^{-}"};

  TFile * fFile = new TFile(Form("../data/AnalysisResults_V0_26_09.root"));

  TH3F* h3dSpectrum = (TH3F *) fFile->Get(Form("correlate-strangeness/h3d%sSpectrum",name[part].Data()));
  h3dSpectrum->Sumw2();

  TH3F* h3dSpectrumY = (TH3F *) fFile->Get(Form("correlate-strangeness/h3d%sSpectrumY",name[part].Data()));
  h3dSpectrumY->Sumw2();

  TH1F* hPvz = (TH1F *) fFile->Get("correlate-strangeness/EventQA/hPvz");

  TFile * fFileEff = new TFile(Form("../data/Efficiency/Efficiency.root"));
  TH1F * pt_eff = (TH1F *) fFileEff->Get(Form("pt_%s_woBckg",name[part].Data()));
  TH1F * pt_eff_Y = (TH1F *) fFileEff->Get(Form("pt_%s_ycut_woBckg_wPV",name[part].Data()));

  TFile * fFilePub = new TFile(Form("../data/Spectra_MB_13TeV.root"));
  TH1F * pubK0s = (TH1F * )fFilePub->Get("Table\ 4/Hist1D_y1");
  TH1F * pubK0s_systErr = (TH1F * )fFilePub->Get("Table\ 4/Hist1D_y1");
  TH1F * pubK0s_stats = (TH1F * )fFilePub->Get("Table\ 4/Hist1D_y1_e1");
  TH1F * pubK0s_syst = (TH1F * )fFilePub->Get("Table\ 4/Hist1D_y1_e2");

  for (size_t i = 1; i < pubK0s->GetXaxis()->GetNbins()+1; i++) {
    pubK0s->SetBinError(i,pubK0s_stats->GetBinContent(i));
    pubK0s_systErr->SetBinError(i,pubK0s_syst->GetBinContent(i));
  }

  h3dSpectrum->GetZaxis()->SetRange(2,2);
  TH1F * signal = (TH1F *) h3dSpectrum->Project3D("x");
  signal->SetName("signal");
  h3dSpectrum->GetZaxis()->SetRange(1,1);
  TH1F * bckg = (TH1F *) h3dSpectrum->Project3D("x");
  bckg->SetName("bckg");
  h3dSpectrum->GetZaxis()->SetRange(3,3);
  bckg->Add ((TH1F *) h3dSpectrum->Project3D("x"));

  signal->Add(bckg,-1);
  signal->Scale(1./hPvz->GetEntries());

  h3dSpectrumY->GetZaxis()->SetRange(2,2);
  TH1F * signalY = (TH1F *) h3dSpectrumY->Project3D("x");
  signalY->SetName("signal");
  h3dSpectrumY->GetZaxis()->SetRange(1,1);
  TH1F * bckgY = (TH1F *) h3dSpectrumY->Project3D("x");
  bckgY->SetName("bckgY");
  h3dSpectrumY->GetZaxis()->SetRange(3,3);
  bckgY->Add ((TH1F *) h3dSpectrumY->Project3D("x"));

  signalY->Add(bckgY,-1);
  signalY->Scale(1./hPvz->GetEntries());

  cout << signal->GetXaxis()->GetNbins() << endl;
  for (size_t i = 1; i < signal->GetXaxis()->GetNbins()+1; i++) {
    if(pt_eff->Interpolate (signal->GetBinCenter(i))==0){
      signal->SetBinContent(i,0.0001);
    }
    else
    {
      signal->SetBinContent(i,signal->GetBinContent(i)/(signal->GetBinWidth(i)*pt_eff->Interpolate (signal->GetBinCenter(i))));
      signal->SetBinError(i,signal->GetBinError(i)/(signal->GetBinWidth(i)*pt_eff->Interpolate (signal->GetBinCenter(i))));
    }
    if(pt_eff_Y->Interpolate (signal->GetBinCenter(i))==0){
      signalY->SetBinContent(i,0.0001);
    }
    else
    {
      // bckgY->SetBinContent(i,bckgY->GetBinContent(i)/(bckgY->GetBinWidth(i)*pt_eff_Y->Interpolate (signalY->GetBinCenter(i))));
      // bckgY->SetBinError(i,bckgY->GetBinError(i)/(bckgY->GetBinWidth(i)*pt_eff_Y->Interpolate (signalY->GetBinCenter(i))));
      signalY->SetBinContent(i,signalY->GetBinContent(i)/(signalY->GetBinWidth(i)*pt_eff_Y->Interpolate (signalY->GetBinCenter(i))));
      signalY->SetBinError(i,signalY->GetBinError(i)/(signalY->GetBinWidth(i)*pt_eff_Y->Interpolate (signalY->GetBinCenter(i))));
    }

  }
  Double_t meanPt = 0.;
  Float_t norm=0.;
  for (size_t i = 1; i < signalY->GetXaxis()->GetNbins()+1; i++) {
    if(signalY->GetBinCenter(i)>10)continue;
    meanPt+=signalY->GetBinCenter(i)*signalY->GetBinContent(i);
    norm+=signalY->GetBinContent(i);
  }

  meanPt= meanPt/norm;

  TCanvas * canSpectrum = Plotter::CreateCanvas("canSpectrum");
  // canSpectrum->GetPadSave()->SetLogy();
  Plotter::SetHistAxes(signal,"#font[12]{p}_{T} (GeV/#font[12]{c})","1/N_{ev}dN/d#font[12]{p}_{T}");
  Plotter::SetHist(signal,"",20,kBlack,1.2);
  signal->GetXaxis()->SetRangeUser(0,10);
  signal->DrawCopy();

  TPaveText *pave = new TPaveText();
  Plotter::SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.95, 0.65,0.95);
  pave->AddText("ALICE, Work in Progress");
  pave->AddText("pp, 13.6 TeV");
  pave->AddText(Form("%s",finalNames[part].Data()));
  pave->AddText("|#eta| < 0.8");
  pave->Draw("same");

  TCanvas * canSpectrumY = Plotter::CreateCanvas("canSpectrumY");
  canSpectrumY->GetPadSave()->SetLogy();
  Plotter::SetHistAxes(signalY,"#font[12]{p}_{T} (GeV/#font[12]{c})","1/N_{ev}dN/d#font[12]{p}_{T}");
  Plotter::SetHist(signalY,"",20,kBlack,1.2);
  Plotter::SetHist(pubK0s_systErr,"",24,kBlack,1.2);
  Plotter::SetHist(pubK0s,"",24,kBlack,1.2);
  signalY->GetXaxis()->SetRangeUser(0,10);
  signalY->DrawCopy();
  pubK0s->DrawCopy("same");
  pubK0s_systErr->DrawCopy("same e2");

  TArrow* arrowMeanPt = new TArrow(meanPt, signalY->GetMinimum()-2, meanPt, signalY->GetMinimum()+0.001);
  arrowMeanPt->SetOption("<");
  arrowMeanPt->SetLineColor(kBlack);
  arrowMeanPt->SetAngle(40);
  arrowMeanPt->Draw("");

  TPaveText *paveY = new TPaveText();
  Plotter::SetPaveText(paveY,42,0.05, 0, 0,33,0,0.55,0.95, 0.65,0.95);
  paveY->AddText("ALICE, Work in Progress");
  paveY->AddText("pp, 13.6 TeV");
  paveY->AddText(Form("%s,",finalNames[part].Data()));
  paveY->AddText("|y| < 0.5");
  paveY->Draw("same");

  TLegend *leg= Plotter::CreateLegend(0.55, 0.92, 0.45, 0.65,0.05);
  leg->AddEntry(signalY,"This analysis","p");
  leg->AddEntry(pubK0s,"pp 13TeV, Eur. Phys. J. C 81 (2021) 256","p");
  leg->Draw();
}
