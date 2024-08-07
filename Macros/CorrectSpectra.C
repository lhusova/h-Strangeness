#include "Plotter.h"

void CorrectSpectra(Int_t part = 3){

  TString name[]={"K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus","Pion"};
  TString finalNames[]={"K_{S}^{0}","#Lambda","#bar{#Lambda}","#Xi^{-}","#Xi^{+}","#Omega^{-}","#Omega^{+})","#pi^{+}+#pi^{-}"};

  TFile * fFile = new TFile(Form("../data/AnalysisResults_Cascades_26_09.root"));

  TH3F* h3dSpectrum = (TH3F *) fFile->Get(Form("correlate-strangeness/h3d%sSpectrum",name[part].Data()));
  h3dSpectrum->Sumw2();

  TH3F* h3dSpectrumY = (TH3F *) fFile->Get(Form("correlate-strangeness/h3d%sSpectrumY",name[part].Data()));
  h3dSpectrumY->Sumw2();

  TH3F* h3dSpectrumYPlus;
  TH3F* h3dSpectrumPlus;
  TH1F * pt_effPlus;
  TH1F * pt_eff_YPlus;
  TH1F * signalPlus;
  TH1F * signalYPlus;

  if(part==3){
    h3dSpectrumPlus = (TH3F *) fFile->Get(Form("correlate-strangeness/h3d%sSpectrum",name[part+1].Data()));
    h3dSpectrumPlus->Sumw2();
    // h3dSpectrum->Add(h3dSpectrumPlus);

    h3dSpectrumYPlus = (TH3F *) fFile->Get(Form("correlate-strangeness/h3d%sSpectrumY",name[part+1].Data()));
    h3dSpectrumYPlus->Sumw2();
    // h3dSpectrumY->Add(h3dSpectrumYPlus);
  }

  TH1F* hPvz = (TH1F *) fFile->Get("correlate-strangeness/EventQA/hPvz");

  TFile * fFileEff = new TFile(Form("../data/Efficiency/Efficiency.root"));
  TH1F * pt_eff = (TH1F *) fFileEff->Get(Form("pt_%s_woBckg",name[part].Data()));
  TH1F * pt_eff_Y = (TH1F *) fFileEff->Get(Form("pt_%s_ycut_woBckg_wPV",name[part].Data()));

  if(part==3){
    pt_effPlus = (TH1F *) fFileEff->Get(Form("pt_%s_woBckg",name[part+1].Data()));
    pt_eff_YPlus = (TH1F *) fFileEff->Get(Form("pt_%s_ycut_woBckg_wPV",name[part+1].Data()));
  }

  TFile * fFileSpec = new TFile(Form("../data/Spectra/SpectrumFromInvMass_K0Short.root",name[part].Data()));
  TH1F * hFinalSpectrumInvMass = (TH1F* ) fFileSpec->Get("hFinalSpectrum");

  TFile * fFilePub = new TFile(Form("../data/Spectra/Spectra_MB_13TeV.root"));
  TH1F * pubK0s = (TH1F * )fFilePub->Get("Table\ 3/Hist1D_y1");
  TH1F * pubK0s_systErr = (TH1F * )fFilePub->Get("Table\ 3/Hist1D_y1");
  TH1F * pubK0s_stats = (TH1F * )fFilePub->Get("Table\ 3/Hist1D_y1_e1");
  TH1F * pubK0s_syst = (TH1F * )fFilePub->Get("Table\ 3/Hist1D_y1_e2");

  TFile * fFileChiaraXi = new TFile(Form("../data/Spectra/YieldEffCorrLHC22o_pass4_MinBias_Train108123_Xi_BkgParab_Mult0-100_Train109827_INELgt0_New.root"));
  TH1F * histXiChiara = (TH1F *) fFileChiaraXi->Get("histoYieldCorr");

  for (size_t i = 1; i < pubK0s->GetXaxis()->GetNbins()+1; i++) {
    pubK0s->SetBinError(i,pubK0s_stats->GetBinContent(i));
    pubK0s_systErr->SetBinError(i,pubK0s_syst->GetBinContent(i));
  }

  h3dSpectrum->GetYaxis()->SetRangeUser(0,100);
  h3dSpectrum->GetZaxis()->SetRange(2,2);
  TH1F * signal = (TH1F *) h3dSpectrum->Project3D("x");
  signal->SetName("signal");
  h3dSpectrum->GetZaxis()->SetRange(1,1);
  TH1F * bckg = (TH1F *) h3dSpectrum->Project3D("x");
  bckg->SetName("bckg");
  h3dSpectrum->GetZaxis()->SetRange(3,3);
  bckg->Add ((TH1F *) h3dSpectrum->Project3D("x"));
  signal->Add(bckg,-1);

  if(part==3){
    h3dSpectrumPlus->GetYaxis()->SetRangeUser(0,100);
    h3dSpectrumPlus->GetZaxis()->SetRange(2,2);
    signalPlus = (TH1F *) h3dSpectrumPlus->Project3D("x");
    signalPlus->SetName("signalPlus");
    h3dSpectrumPlus->GetZaxis()->SetRange(1,1);
    TH1F * bckgPlus = (TH1F *) h3dSpectrumPlus->Project3D("x");
    bckgPlus->SetName("bckgPlus");
    h3dSpectrumPlus->GetZaxis()->SetRange(3,3);
    bckgPlus->Add ((TH1F *) h3dSpectrumPlus->Project3D("x"));
    signalPlus->Add(bckgPlus,-1);
  }


  h3dSpectrumY->GetYaxis()->SetRangeUser(0,100);
  h3dSpectrumY->GetZaxis()->SetRange(2,2);
  TH1F * signalY = (TH1F *) h3dSpectrumY->Project3D("x");
  signalY->SetName("signal");
  h3dSpectrumY->GetZaxis()->SetRange(1,1);
  TH1F * bckgY = (TH1F *) h3dSpectrumY->Project3D("x");
  bckgY->SetName("bckgY");
  h3dSpectrumY->GetZaxis()->SetRange(3,3);
  bckgY->Add ((TH1F *) h3dSpectrumY->Project3D("x"));
  signalY->Add(bckgY,-1);

  if(part==3){
    h3dSpectrumYPlus->GetYaxis()->SetRangeUser(0,100);
    h3dSpectrumYPlus->GetZaxis()->SetRange(2,2);
    signalYPlus = (TH1F *) h3dSpectrumYPlus->Project3D("x");
    signalYPlus->SetName("signalYPlus");
    h3dSpectrumYPlus->GetZaxis()->SetRange(1,1);
    TH1F * bckgYPlus = (TH1F *) h3dSpectrumYPlus->Project3D("x");
    bckgYPlus->SetName("bckgYPlus");
    h3dSpectrumYPlus->GetZaxis()->SetRange(3,3);
    bckgYPlus->Add ((TH1F *) h3dSpectrumYPlus->Project3D("x"));
    signalYPlus->Add(bckgYPlus,-1);
  }


  TH1F * signalCopy = (TH1F *) hFinalSpectrumInvMass->Clone();
  signalCopy->SetName("signalCopy");
  TH1F * signalYCopy = (TH1F *) hFinalSpectrumInvMass->Clone();
  signalYCopy->SetName("signalYCopy");

  for (size_t i = 0; i < hFinalSpectrumInvMass->GetXaxis()->GetNbins(); i++) {
    signalCopy->SetBinContent(i+1,signal->GetBinContent(i+1)/signal->GetBinWidth(i+1));
    signalCopy->SetBinError(i+1,signal->GetBinError(i+1)/signal->GetBinWidth(i+1));
    signalYCopy->SetBinContent(i+1,signalY->GetBinContent(i+1)/signalY->GetBinWidth(i+1));
    signalYCopy->SetBinError(i+1,signalY->GetBinError(i+1)/signalY->GetBinWidth(i+1));
  }

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
  if(part==3) {
    for (size_t i = 1; i < signal->GetXaxis()->GetNbins()+1; i++) {
      if(pt_effPlus->Interpolate (signal->GetBinCenter(i))==0){
        signalPlus->SetBinContent(i,0.0001);
      }
      else
      {
        signalPlus->SetBinContent(i,signalPlus->GetBinContent(i)/(signalPlus->GetBinWidth(i)*pt_effPlus->Interpolate (signal->GetBinCenter(i))));
        signalPlus->SetBinError(i,signalPlus->GetBinError(i)/(signalPlus->GetBinWidth(i)*pt_effPlus->Interpolate (signal->GetBinCenter(i))));
      }
      if(pt_eff_YPlus->Interpolate (signalPlus->GetBinCenter(i))==0){
        signalYPlus->SetBinContent(i,0.0001);
      }
      else
      {
        signalYPlus->SetBinContent(i,signalYPlus->GetBinContent(i)/(signalYPlus->GetBinWidth(i)*pt_eff_YPlus->Interpolate (signalYPlus->GetBinCenter(i))));
        signalYPlus->SetBinError(i,signalYPlus->GetBinError(i)/(signalYPlus->GetBinWidth(i)*pt_eff_YPlus->Interpolate (signalYPlus->GetBinCenter(i))));
      }
    }
    signal->Add(signalPlus);
    signalY->Add(signalYPlus);
  }
  signal->Scale(1./hPvz->GetEntries());
  signalY->Scale(1./hPvz->GetEntries());

  Double_t meanPt = 0.;
  Float_t norm=0.;
  for (size_t i = 1; i < signalY->GetXaxis()->GetNbins()+1; i++) {
    if(signalY->GetBinCenter(i)>10)continue;
    meanPt+=signalY->GetBinCenter(i)*signalY->GetBinContent(i);
    norm+=signalY->GetBinContent(i);
  }

  meanPt= meanPt/norm;

  TCanvas * canSpectrum = Plotter::CreateCanvas("canSpectrum");
  canSpectrum->GetPadSave()->SetLogy();
  Plotter::SetHistAxes(signal,"#font[12]{p}_{T} (GeV/#font[12]{c})","1/N_{ev}dN/d#font[12]{p}_{T}");
  Plotter::SetHist(signal,"",20,kBlack,1.2);
  Plotter::SetHist(hFinalSpectrumInvMass,"",24,kBlack,1.2);
  Plotter::SetHist(signalYCopy,"",20,kRed,1.2);
  signal->GetXaxis()->SetRange(5,40);
  // signalY->GetXaxis()->SetRangeUser(0.5,10);
  signal->DrawCopy();
  // hFinalSpectrumInvMass->DrawCopy("same");
  // signalYCopy->DrawCopy("same");

  // TCanvas * canRatio = Plotter::CreateCanvas("canRatio");
  // signalCopy->Divide(hFinalSpectrumInvMass);
  // signalCopy->DrawCopy();
  // Plotter::SetHist(signalYCopy,"",20,kRed,1.2);
  // signalYCopy->Divide(hFinalSpectrumInvMass);
  // signalYCopy->DrawCopy("same");

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
  if(part==3) histXiChiara->DrawCopy("same");
  if(part==0){
    pubK0s->DrawCopy("same");
    pubK0s_systErr->DrawCopy("same e2");
  }

  // TArrow* arrowMeanPt = new TArrow(meanPt, signalY->GetMinimum()-2, meanPt, signalY->GetMinimum()+0.001);
  // arrowMeanPt->SetOption("<");
  // arrowMeanPt->SetLineColor(kBlack);
  // arrowMeanPt->SetAngle(40);
  // arrowMeanPt->Draw("");

  Float_t particleMass[]={0.497,1.115,1.115,1.321,1.321,1.672,0.1396};
  TF1 * levy = new TF1("levy","[3]*x/([0]*[1])*(([0]-1)*([0]-2))/([0]*[1]+[2]*([0]-2))*TMath::Power(1+(TMath::Sqrt([2]*[2]+x*x)-[2])/([0]*[1]),-[0])",0.6,6);
  levy->SetParameter(0,7);
  levy->SetParameter(1,0.8);
  levy->FixParameter(2,particleMass[part]);
  levy->SetParameter(3,0.03);

  signalY->Fit(levy,"R0");

  levy->SetLineColor(kBlue);
  levy->Draw("same");

  TPaveText *paveY = new TPaveText();
  Plotter::SetPaveText(paveY,42,0.05, 0, 0,33,0,0.55,0.95, 0.65,0.95);
  paveY->AddText("ALICE, Work in Progress");
  paveY->AddText("pp, 13.6 TeV");
  if(part==3)paveY->AddText(Form("%s+%s,",finalNames[part].Data(),finalNames[part+1].Data()));
  else paveY->AddText(Form("%s,",finalNames[part].Data()));
  paveY->AddText("|y| < 0.5");
  paveY->Draw("same");

  TLegend *leg= Plotter::CreateLegend(0.55, 0.92, 0.45, 0.65,0.05);
  leg->AddEntry(signalY,"This analysis","p");
  if(part==0) leg->AddEntry(pubK0s,"pp 13TeV, Eur. Phys. J. C 81 (2021) 256","p");
  if(part==3) leg->AddEntry(histXiChiara,"Chiara et al., 13.6 TeV","p");
  leg->AddEntry(levy,"Fit: This analysis","l");
  leg->Draw();

  TCanvas * canDataToFit = Plotter::CreateCanvas("canDataToFit");
  TH1F * signalYRatio = (TH1F *) signalY->Clone();
  signalYRatio->SetName("signalYRatio");
  signalYRatio->Divide(levy);
  Plotter::SetHistAxes(signalYRatio,"#font[12]{p}_{T} (GeV/#font[12]{c})","Data/Fit");
  signalYRatio->DrawCopy();
  TLegend *legRatio= Plotter::CreateLegend(0.55, 0.92, 0.45, 0.65,0.05);
  legRatio->AddEntry(signalYRatio,"This analysis","p");
  if(part==0){
    TH1F * pubK0sRatio = (TH1F *) pubK0s->Clone();
    pubK0sRatio->SetName("signalYRatio");
    pubK0sRatio->Divide(levy);
    pubK0sRatio->DrawCopy("same");
    legRatio->AddEntry(pubK0sRatio,"pp 13TeV, Eur. Phys. J. C 81 (2021) 256","p");
  }
  if(part==3){
    TH1F * histXiChiaraRatio = (TH1F *) histXiChiara->Clone();
    histXiChiaraRatio->SetName("histXiChiaraRatio");
    histXiChiaraRatio->Divide(levy);
    histXiChiaraRatio->DrawCopy("same");
    if(part==3) legRatio->AddEntry(histXiChiaraRatio,"Chiara et al., 13.6 TeV","p");
  }
  paveY->Draw("same");

  Plotter::DrawUnity(kBlue, 1.,0.,10);


  legRatio->Draw();

}
