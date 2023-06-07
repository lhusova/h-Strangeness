#include "Plotter.h"

void FitDeltaPhiProjection(Int_t iPart = 0,Int_t iPt =6, Bool_t sigma = kTRUE, TString function = "vMises"){

  if(iPt>6) {
    cout<< "Error!!! Pt range not defined !!!" << endl;
    return;
  }

  if(iPart>3) {
    cout<< "Error!!! Assoc particle not defined !!!" << endl;
    return;
  }

  TString particleName[]={"K0s","Lam","Xi","Omega"};
  TString finalNames[]={"K_{S}^{0}","(#Lambda+#bar{#Lambda})","(#Xi^{+}+#Xi^{-})","(#Omega^{+}+#Omega^{-})"};
  Float_t ptRanges[]={0.5,1.,1.5,2.,3.,4.,6.,10.};

  TString multiplicityNames[]={"minBias","0_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_60Mult","60_70Mult","70_80Mult","80_90Mult","90_100Mult"};
  TString multiplicityPave[]={"MB","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"};

  TF1 * funcMisesNear = new TF1("funcMisesNear", "[0]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x - 2*TMath::Pi() - [1]))",-1.2,1.2);
  funcMisesNear->SetParNames("norm","mu","kappa");
  funcMisesNear->SetParameter(1,0);

  TF1 * funcMisesAway = new TF1("funcMisesAway", "[0]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x - 2*TMath::Pi() - [1]))",TMath::Pi()/2,+3./2.*TMath::Pi());
  funcMisesAway->SetParNames("norm","mu","kappa");
  funcMisesAway->SetParameter(1,TMath::Pi());

  TF1 * funcGaussNear = new TF1("funcGaussNear", "gaus",-1.2,1.2);
  funcMisesNear->SetParameter(1,0);

  TF1 * funcGaussAway = new TF1("funcGaussAway", "gaus",TMath::Pi()/2,+3./2.*TMath::Pi());
  funcGaussAway->SetParameter(1,TMath::Pi());

  TFile * file[11];
  TH1F* fProj[11];
  TCanvas *can[11];
  TPaveText *pave[11];
  Double_t fwhm_near[11];
  Double_t fwhm_away[11];
  Double_t sigma2_near, sigma2_away, max_near, max_away;

  const char *labels[11] = {"90-100%","80-90%","70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","0-10%","MB"};
  TH1F *histWidt_Near = new TH1F("histWidt_Near","",11,0,11);
  TH1F *histWidt_Away = new TH1F("histWidt_Away","",11,0,11);

  for (Int_t iMult = 0; iMult < 11; iMult++) {
    file[iMult] = new TFile(Form("../data/Yields_%s_%s.root",particleName[iPart].Data(),multiplicityNames[iMult].Data()));
    fProj[iMult] = (TH1F *) file[iMult] -> Get(Form("phiProj_noBckg_pT%d",iPt));
    fProj[iMult]->RebinX(2);
    if(function=="vMises"){
      fProj[iMult]->Fit(funcMisesNear,"R");
      fProj[iMult]->Fit(funcMisesAway,"R");
    }else if(function=="gaus"){
      fProj[iMult]->Fit(funcGaussNear,"R");
      fProj[iMult]->Fit(funcGaussAway,"R");
    }

    can[iMult]= Plotter::CreateCanvas(Form("c%i",iMult));
    Plotter::SetHist(fProj[iMult],"",20,kBlack,0.7);
    fProj[iMult]->GetYaxis()-> SetMaxDigits(2);
    fProj[iMult]->GetYaxis()->SetTitleOffset(0.8);
    fProj[iMult]->DrawCopy();
    if(function=="vMises"){
      funcMisesNear->Draw("same");
      funcMisesAway->SetLineColor(kBlue);
      funcMisesAway->Draw("same");
    }else if(function=="gaus"){
      funcGaussNear->Draw("same");
      funcGaussAway->SetLineColor(kBlue);
      funcGaussAway->Draw("same");
    }

    pave[iMult] = new TPaveText();
    Plotter::SetPaveText(pave[iMult],42,0.05, 0, 0,33,0,0.55,0.97, 0.5,0.95);
    pave[iMult]->AddText("ALICE, Work in Progress");
    pave[iMult]->AddText("pp, 13.6 TeV");
    pave[iMult]->AddText(Form("%s",multiplicityPave[iMult].Data()));
    pave[iMult]->AddText(Form("h-%s",finalNames[iPart].Data()));
    pave[iMult]->AddText(Form("3 < #font[12]{p}^{trigg}_{T} < 20 GeV/#font[12]{c}"));
    pave[iMult]->AddText(Form("%g < #font[12]{p}^{assoc}_{T} < %g GeV/#font[12]{c}",ptRanges[iPt],ptRanges[iPt+1]));
    pave[iMult]->AddText("|#Delta#eta| < 1");
    pave[iMult]->Draw("same");

    can[iMult]->SaveAs(Form("../Plots/Fits/%s/DeltaPhiFit_%i_%s_%s.pdf",particleName[iPart].Data(),iPt,multiplicityNames[iMult].Data(),function.Data()));

    if(sigma){
      if(function=="vMises")sigma2_near = 1./funcMisesNear->GetParameter(2);
      else if(function=="gaus")sigma2_near = TMath::Power(funcGaussNear->GetParameter(2),2);
      fwhm_near[iMult] = 2*TMath::Sqrt(2*TMath::Log(2)*sigma2_near);

      if(function=="vMises")sigma2_away = 1./funcMisesAway->GetParameter(2);
      else if(function=="gaus")sigma2_away = TMath::Power(funcGaussAway->GetParameter(2),2);
      fwhm_away[iMult] = 2*TMath::Sqrt(2*TMath::Log(2)*sigma2_away);
    }else {
      if(function=="vMises"){
        max_near = funcMisesNear->GetMaximum();
        fwhm_near[iMult] = funcMisesNear->GetX(max_near/2,0.1,1) - funcMisesNear->GetX(max_near/2,-1.,-0.1);

        max_away = funcMisesAway->GetMaximum();
        fwhm_away[iMult] = funcMisesAway->GetX(max_away/2,3.15,5) - funcMisesAway->GetX(max_away/2,1,3.1);
      }else if(function=="gaus"){
        max_near = funcGaussNear->GetMaximum();
        fwhm_near[iMult] = funcGaussNear->GetX(max_near/2,0.1,1) - funcGaussNear->GetX(max_near/2,-1.,-0.1);

        max_away = funcGaussAway->GetMaximum();
        fwhm_away[iMult] = funcGaussAway->GetX(max_away/2,3.15,5) - funcGaussAway->GetX(max_away/2,1,3.1);
      }

    }


    histWidt_Near->GetXaxis()->SetBinLabel(iMult+1,labels[iMult]);
    if(iPart==3&&iMult>4)histWidt_Near->SetBinContent(11-iMult,0);
    else histWidt_Near->SetBinContent(11-iMult,fwhm_near[iMult]);
    histWidt_Near->SetBinError(11-iMult,0.0001);

    histWidt_Away->GetXaxis()->SetBinLabel(iMult+1,labels[iMult]);
    if(iPart==3&&iMult>3) histWidt_Away->SetBinContent(11-iMult,0);
    else histWidt_Away->SetBinContent(11-iMult,fwhm_away[iMult]);
    histWidt_Away->SetBinError(11-iMult,0.0001);
  }

  TCanvas *canNear= Plotter::CreateCanvas("canNear");
  histWidt_Near->GetYaxis()->SetRangeUser(0.45,2.1999);
  Plotter::SetHistAxes(histWidt_Near,"","FWHM");
  Plotter::SetHist(histWidt_Near,"",20,kRed+1,1.);
  histWidt_Near->DrawCopy();

  Plotter::SetHist(histWidt_Away,"",20,kBlue+1,1.);
  histWidt_Away->DrawCopy("same");

  TPaveText *paveWidth = new TPaveText();
  Plotter::SetPaveText(paveWidth,42,0.05, 0, 0,33,0,0.55,0.97, 0.63,0.95);
  paveWidth->AddText("ALICE, Work in Progress");
  paveWidth->AddText("pp, 13.6 TeV");
  paveWidth->AddText(Form("h-%s",finalNames[iPart].Data()));
  paveWidth->AddText(Form("3 < #font[12]{p}^{trigg}_{T} < 20 GeV/#font[12]{c}"));
  paveWidth->AddText(Form("%g < #font[12]{p}^{assoc}_{T} < %g GeV/#font[12]{c}",ptRanges[iPt],ptRanges[iPt+1]));
  paveWidth->AddText("|#Delta#eta| < 1");
  paveWidth->Draw("same");

  TLegend *leg = Plotter::CreateLegend(0.15, 0.45, 0.7, 0.95,0.05);
  leg->AddEntry(histWidt_Near,"Near-side","pl");
  leg->AddEntry(histWidt_Away,"Away-side","pl");
  leg->Draw("same");

  if(sigma)canNear->SaveAs(Form("../Plots/Fits/%s/Width_%i_sigma_%s.pdf",particleName[iPart].Data(),iPt,function.Data()));
  else canNear->SaveAs(Form("../Plots/Fits/%s/Width_%i_max_%s.pdf",particleName[iPart].Data(),iPt,function.Data()));

}
