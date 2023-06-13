#include "Plotter.h"

void FitDeltaPhiProjection(Int_t iPart = 3,Int_t iPt =3, Bool_t sigma = kFALSE, TString function = "doublegaus"){

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

  TF1 * funcMisesNear = new TF1("funcMisesNear", "[0]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x - 2*TMath::Pi() - [1]))",-1.4,1.4);
  funcMisesNear->SetParNames("norm","mu","kappa");
  // funcMisesNear->SetParameter(0,1);
  funcMisesNear->SetParameter(1,0);
  // funcMisesNear->SetParLimits(2,1,20);

  TF1 * funcMisesAway = new TF1("funcMisesAway", "[0]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x - 2*TMath::Pi() - [1]))",TMath::Pi()/2-0.3,+3./2.*TMath::Pi());
  funcMisesAway->SetParNames("norm","mu","kappa");
  funcMisesAway->SetParameter(1,TMath::Pi());

  TF1 * funcGaussNear = new TF1("funcGaussNear", "gaus",-0.4,0.4);
  funcGaussNear->SetParameter(1,0);

  TF1 * funcGaussAway = new TF1("funcGaussAway", "gaus",TMath::Pi()/2-0.3,+3./2.*TMath::Pi());
  funcGaussAway->SetParameter(1,TMath::Pi());

  TF1 * funcGenGausNear = new TF1("funcGenGausNear","[2]/([0]*TMath::Sqrt(TMath::Pi()))*TMath::Exp(-TMath::Power((x-[1])/[0],2))",-1.2,1.2);
  funcGenGausNear->SetParameter(0,0.5);
  funcGenGausNear->SetParameter(1,0);
  funcGenGausNear->SetParameter(2,0.05);

  TF1 * funcDoubleGausNear = new TF1("funcDoubleGausNear","[0]*exp(-0.5*((x-[1])/[2]*(x-[1])/[2])) + [3]*exp(-0.5*((x-[1])/[4]*(x-[1])/[4]))",-1.2,1.2);
  // funcDoubleGausNear->SetParameter(0,0.02);
  funcDoubleGausNear->SetParameter(1,0);
  funcDoubleGausNear->SetParameter(2,0.2);
  funcDoubleGausNear->SetParLimits(2,0,1);
  funcDoubleGausNear->SetParameter(4,0.3);
  funcDoubleGausNear->SetParLimits(4,0,1);

  TFile * file[11];
  TH1F* fProj[11];
  TCanvas *can[11];
  TPaveText *pave[11];
  Double_t fwhm_near[11];
  Double_t fwhm_away[11];
  Double_t sigma2_near, sigma2_away, max_near, max_away, chi2_near, chi2_away;

  const char *labels[11] = {"90-100%","80-90%","70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","0-10%","MB"};
  TH1F *histWidth_Near = new TH1F("histWidth_Near","",11,0,11);
  TH1F *histWidth_Away = new TH1F("histWidth_Away","",11,0,11);

  TH1F *histChi2_Near = new TH1F("histChi2_Near","",11,0,11);
  TH1F *histChi2_Away = new TH1F("histChi2_Away","",11,0,11);

  for (Int_t iMult = 0; iMult < 11; iMult++) {
    file[iMult] = new TFile(Form("../data/Yields_%s_%s.root",particleName[iPart].Data(),multiplicityNames[iMult].Data()));
    fProj[iMult] = (TH1F *) file[iMult] -> Get(Form("phiProj_noBckg_pT%d",iPt));
    fProj[iMult]->RebinX(2);
    if(function=="vMises"){
      fProj[iMult]->Fit(funcMisesNear,"MER");
      fProj[iMult]->Fit(funcMisesAway,"MER");
    }else if(function=="gaus"){
      fProj[iMult]->Fit(funcGaussNear,"MER");
      fProj[iMult]->Fit(funcGaussAway,"MER");
    }else if(function=="Gengaus"){
      fProj[iMult]->Fit(funcGenGausNear,"MER");
      fProj[iMult]->Fit(funcGaussAway,"MER");
    }else if(function=="doublegaus"){
      if(iMult==7) {
        // funcDoubleGausNear->SetParameter(0,3.91938e-03);
        // funcDoubleGausNear->SetParameter(3,3.31162e-03);
        // funcDoubleGausNear->SetParameter(2,2.55384e-01);
      }
      fProj[iMult]->Fit(funcDoubleGausNear,"MER");
      fProj[iMult]->Fit(funcGaussAway,"MER");
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
      chi2_near = funcMisesNear->GetChisquare()/funcMisesNear->GetNDF();
      chi2_away = funcMisesAway->GetChisquare()/funcMisesAway->GetNDF();
    }else if(function=="gaus"){
      funcGaussNear->Draw("same");
      funcGaussAway->SetLineColor(kBlue);
      funcGaussAway->Draw("same");
      chi2_near = funcGaussNear->GetChisquare()/funcGaussNear->GetNDF();
      chi2_away = funcGaussAway->GetChisquare()/funcGaussAway->GetNDF();
    }else if(function=="Gengaus"){
      funcGenGausNear->Draw("same");
      funcGaussAway->SetLineColor(kBlue);
      funcGaussAway->Draw("same");
      chi2_near = funcGenGausNear->GetChisquare()/funcGenGausNear->GetNDF();
      chi2_away = funcGaussAway->GetChisquare()/funcGaussAway->GetNDF();
    }else if(function=="doublegaus"){
      funcDoubleGausNear->Draw("same");
      funcGaussAway->SetLineColor(kBlue);
      funcGaussAway->Draw("same");
      chi2_near = funcDoubleGausNear->GetChisquare()/funcDoubleGausNear->GetNDF();
      chi2_away = funcGaussAway->GetChisquare()/funcGaussAway->GetNDF();
    }
    cout << "____________________________________________________________________" << chi2_near << endl;
    histChi2_Near->SetBinContent(11-iMult,chi2_near);
    histChi2_Away->SetBinContent(11-iMult,chi2_away);
    histChi2_Near->SetBinError(11-iMult,0.001);
    histChi2_Away->SetBinError(11-iMult,0.001);
    histChi2_Near->GetXaxis()->SetBinLabel(iMult+1,labels[iMult]);

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
      }else if(function=="Gengaus"){
        max_near = funcGenGausNear->GetMaximum();
        fwhm_near[iMult] = funcGenGausNear->GetX(max_near/2,0.1,1) - funcGenGausNear->GetX(max_near/2,-1.,-0.1);

        max_away = funcGaussAway->GetMaximum();
        fwhm_away[iMult] = funcGaussAway->GetX(max_away/2,3.15,5) - funcGaussAway->GetX(max_away/2,1,3.1);
      }else if(function=="doublegaus"){
        max_near = funcDoubleGausNear->GetMaximum();
        fwhm_near[iMult] = funcDoubleGausNear->GetX(max_near/2,0.1,1) - funcDoubleGausNear->GetX(max_near/2,-1.,-0.1);

        max_away = funcGaussAway->GetMaximum();
        fwhm_away[iMult] = funcGaussAway->GetX(max_away/2,3.15,5) - funcGaussAway->GetX(max_away/2,1,3.1);
      }

    }


    histWidth_Near->GetXaxis()->SetBinLabel(iMult+1,labels[iMult]);
    if(iPart==3&&iMult>4)histWidth_Near->SetBinContent(11-iMult,0);
    else histWidth_Near->SetBinContent(11-iMult,fwhm_near[iMult]);
    histWidth_Near->SetBinError(11-iMult,0.0001);

    histWidth_Away->GetXaxis()->SetBinLabel(iMult+1,labels[iMult]);
    if(iPart==3&&iMult>3) histWidth_Away->SetBinContent(11-iMult,0);
    else histWidth_Away->SetBinContent(11-iMult,fwhm_away[iMult]);
    histWidth_Away->SetBinError(11-iMult,0.0001);
  }

  TCanvas *canNear= Plotter::CreateCanvas("canNear");
  histWidth_Near->GetYaxis()->SetRangeUser(0.45,2.1999);
  Plotter::SetHistAxes(histWidth_Near,"","FWHM");
  Plotter::SetHist(histWidth_Near,"",20,kRed+1,1.);
  histWidth_Near->DrawCopy();

  Plotter::SetHist(histWidth_Away,"",20,kBlue+1,1.);
  histWidth_Away->DrawCopy("same");

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
  leg->AddEntry(histWidth_Near,"Near-side","pl");
  leg->AddEntry(histWidth_Away,"Away-side","pl");
  leg->Draw("same");

  if(sigma)canNear->SaveAs(Form("../Plots/Fits/%s/Width_%i_sigma_%s.pdf",particleName[iPart].Data(),iPt,function.Data()));
  else canNear->SaveAs(Form("../Plots/Fits/%s/Width_%i_max_%s.pdf",particleName[iPart].Data(),iPt,function.Data()));

  TCanvas *canChi2= Plotter::CreateCanvas("canChi2");
  Plotter::SetHistAxes(histChi2_Near,"","#chi^{2}/NDF");
  Plotter::SetHist(histChi2_Near,"",20,kRed+1,1.);
  histChi2_Near->DrawCopy();

  Plotter::SetHist(histChi2_Away,"",20,kBlue+1,1.);
  histChi2_Away->DrawCopy("same");
  canChi2->SaveAs(Form("../Plots/Fits/%s/Chi2_%i_sigma_%s.pdf",particleName[iPart].Data(),iPt,function.Data()));
}
