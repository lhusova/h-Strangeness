#include "Plotter.h"

void CalculateRawSpectraFromInvMass(Int_t part = 0){

  TString name[]={"K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus","Pion"};
  TString finalNames[]={"K_{S}^{0}","#Lambda","#bar{#Lambda}","#Xi^{-}","#Xi^{+}","#Omega^{-}","#Omega^{+})","#pi^{+}+#pi^{-}"};

  TFile * fFile = new TFile(Form("../data/AnalysisResults_V0_26_09.root"));

  TH3F * hMass = (TH3F *) fFile->Get(Form("hstrangecorrelationfilter/h3dMass%s",name[part].Data()));
  hMass->Sumw2();

  TH1F * hEventCounter = (TH1F *) fFile->Get("event-selection-task/hColCounterAcc");

  TH2F * h2dMass = (TH2F *)hMass->Project3D("yx");
  h2dMass->DrawCopy("colz");

  Float_t pt[] = {0.,0.1f, 0.2f,0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f};
  const Int_t nPtBins = sizeof(pt) / sizeof(Float_t)-1;

  TH1F *h1DMass[nPtBins];
  TCanvas *canMass[nPtBins];

  TF1 * fitFunc = new TF1("fitFunc","[3]+[4]*x+[2]*TMath::Gaus(x,[0],[1])",0.4,0.6);
  fitFunc->SetParameter(0,0.49);
  fitFunc->SetParLimits(0,0.47,0.52);
  fitFunc->SetParameter(1,0.02);
  // fitFunc->SetParLimits(1,0.47,0.52);
  fitFunc->SetParameter(2,20000);

  Double_t signal,sideleft,sideright,error,errorleft,errorright;

  TH1F* hFinalSpectrum = new TH1F("hFinalSpectrum","",nPtBins,pt);
  for (size_t iPt = 0; iPt < nPtBins; iPt++) {
    if(iPt==0){
      hFinalSpectrum->SetBinContent(iPt+1,-9999);
      hFinalSpectrum->SetBinError(iPt+1,1);
      continue;
    }
    h2dMass->GetXaxis()->SetRangeUser(pt[iPt],pt[iPt+1]);
    canMass[iPt] = Plotter::CreateCanvas(Form("canMass%i",iPt));
    h1DMass[iPt] = (TH1F *)h2dMass->ProjectionY();
    h1DMass[iPt]->Scale(1./h1DMass[iPt]->GetXaxis()->GetBinWidth(10));
    Plotter::SetHistAxes(h1DMass[iPt],"m (GeV/#font[12]{c}^{2})","Counts / (0.001 GeV/#font[12]{c}^{2})");
    Plotter::SetHist(h1DMass[iPt],"",20,kBlack,1);
    h1DMass[iPt]->DrawCopy();
    h1DMass[iPt]->Fit(fitFunc);

    signal = h1DMass[iPt]->IntegralAndError(h1DMass[iPt]->FindBin(fitFunc->GetParameter(0)-3*fitFunc->GetParameter(1)),h1DMass[iPt]->FindBin(fitFunc->GetParameter(0)+3*fitFunc->GetParameter(1)),error,"width");
    sideleft = h1DMass[iPt]->IntegralAndError(h1DMass[iPt]->FindBin(fitFunc->GetParameter(0)-8*fitFunc->GetParameter(1)),h1DMass[iPt]->FindBin(fitFunc->GetParameter(0)-5*fitFunc->GetParameter(1)),errorleft,"width");
    sideright = h1DMass[iPt]->IntegralAndError(h1DMass[iPt]->FindBin(fitFunc->GetParameter(0)+5*fitFunc->GetParameter(1)),h1DMass[iPt]->FindBin(fitFunc->GetParameter(0)+8*fitFunc->GetParameter(1)),errorright,"width");

    hFinalSpectrum->SetBinContent(iPt+1,(signal-sideleft-sideright)/hFinalSpectrum->GetBinWidth(iPt+1));
    hFinalSpectrum->SetBinError(iPt+1,TMath::Sqrt(error*error+errorleft*errorleft+errorright*errorright)/hFinalSpectrum->GetBinWidth(iPt+1));
  }

  hFinalSpectrum->Scale(1./hEventCounter->GetBinContent(1));

  TCanvas * canSpec = Plotter::CreateCanvas("canSpec");
  Plotter::SetHistAxes(hFinalSpectrum,"#font[12]{p}_{T} (GeV/#font[12]{c})","1/N_{ev} dN/d#font[12]{p}_{T}");
  Plotter::SetHist(hFinalSpectrum,"",20,kBlack,1);
  hFinalSpectrum->DrawCopy();

  TFile *fileNew = new TFile(Form("../data/Spectra/SpectrumFromInvMass_%s.root",name[part].Data()), "RECREATE");
  hFinalSpectrum->Write();
  fileNew->Close();
}
