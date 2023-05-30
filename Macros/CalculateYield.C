#include "Plotter.h"

TH1F * GetBachgroundHist(TH1F* hist);

void CalculateYield(Int_t part=0){

  Int_t particleType =part;
  if(part==2)particleType=3;
  if(part==3)particleType=5;
  TString name[]={"K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus"};
  TString finalNames[]={"K_{S}^{0}","(#Lambda+#bar{#Lambda})","(#Xi^{+}+#Xi^{-})","(#Omega^{+}+#Omega^{-})"};
  TString nameSave[]={"K0s","Lam","Xi","Omega"};

  TFile * fFile[3];
  TFile * fFile2[3];
  TString InvMassRanges[] = {"Signal", "LeftBg", "RightBg"};
  for (Int_t i = 0; i < 3; i++) {
    fFile[i] = new TFile(Form("../data/MixCorrected/%s/MixCorrected_%s_%s_minBias.root",name[particleType].Data(),InvMassRanges[i].Data(),name[particleType].Data()));
    if(part>0) fFile2[i] = new TFile(Form("../data/MixCorrected/%s/MixCorrected_%s_%s_minBias.root",name[particleType].Data(),InvMassRanges[i].Data(),name[particleType+1].Data()));
  }
  TFile * fFileTrigger = new TFile("../data/AnalysisResults_Hyperloop_15_05.root");
  TH1F* histTriggers;
  if(part<3) histTriggers = (TH1F*) fFileTrigger->Get("correlate-strangeness/sameEvent/TriggerParticlesV0");
  else histTriggers = (TH1F*) fFileTrigger->Get("correlate-strangeness/sameEvent/TriggerParticlesCascades");

  Float_t nTrigg= histTriggers->Integral();

  const Int_t nPtBins = 7;
  // Double_t ptBins[nPtBins]={0.75,1.25,1.75,2.5,3.5,5.,8.};
  Double_t ptBins[nPtBins+1]={0.5,1.,1.5,2.,3.,4.,6.,10.};
  Double_t ptBins_err[nPtBins]={0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  TH2F * fHistCorrected[nPtBins][3];
  TH1F * fHistCorrectedProjection[nPtBins];
  TH1F * fHistBckg[nPtBins];
  Double_t yield_UE[nPtBins];
  Double_t yield_UE_error[nPtBins];
  Double_t yields[2][nPtBins];
  Double_t yields_err[2][nPtBins];

  TH1D * fHistNear = new TH1D("fHistNear","",nPtBins,ptBins);
  TH1D * fHistAway = new TH1D("fHistAway","",nPtBins,ptBins);
  TH1D * fHistUE = new TH1D("fHistUE","",nPtBins,ptBins);

  TFile * fFileNew = TFile::Open (Form("../data/Yields_%s.root",nameSave[part].Data()),"RECREATE");

  for (Int_t iPt = 0; iPt < nPtBins; iPt++) {
    for (Int_t iFile = 0; iFile < 3; iFile++) { // side-band subtraction
      fHistCorrected[iPt][iFile] = (TH2F*)fFile[iFile]->Get(Form("fHistCorrected_%s_pt%d",name[particleType].Data(),iPt));
      if(part>0)fHistCorrected[iPt][iFile]->Add ((TH2F*)fFile2[iFile]->Get(Form("fHistCorrected_%s_pt%d",name[particleType+1].Data(),iPt)));
      if(iFile>0) fHistCorrected[iPt][0]->Add(fHistCorrected[iPt][iFile],-1);
    }

    fHistCorrected[iPt][0]->GetXaxis()->SetRangeUser(-1,1);

    fHistCorrectedProjection[iPt] = (TH1F *) fHistCorrected[iPt][0]->ProjectionY();
    fHistCorrectedProjection[iPt]->Scale(fHistCorrected[iPt][0]->GetXaxis()->GetBinWidth(2));
    fHistCorrectedProjection[iPt]->Scale(1./nTrigg);

    Plotter::SetHistAxes(fHistCorrectedProjection[iPt],"#Delta#varphi","#frac{1}{N_{trigg}} #frac{dN}{d#Delta#varphi}");

    fHistCorrectedProjection[iPt]->SetName(Form("phiProj_withBckg_pT%d",iPt));
    fHistCorrectedProjection[iPt]->Write();

    fHistBckg[iPt] = GetBachgroundHist(fHistCorrectedProjection[iPt]);
    fHistBckg[iPt]->SetName(Form("underlying_ev_pT%d",iPt));

    fHistBckg[iPt]->Write();
    yield_UE[iPt]=fHistBckg[iPt]->IntegralAndError(fHistBckg[iPt]->GetXaxis()->GetFirst(),fHistBckg[iPt]->GetXaxis()->GetLast(),yield_UE_error[iPt],"width");
    yield_UE[iPt]=yield_UE[iPt]/50;
    yield_UE_error[iPt]=yield_UE_error[iPt]/50;
    fHistUE->SetBinContent(iPt+1,yield_UE[iPt]);
    fHistUE->SetBinError(iPt+1,yield_UE_error[iPt]);

    fHistCorrectedProjection[iPt]->Add(fHistBckg[iPt],-1);
    fHistCorrectedProjection[iPt]->SetName(Form("phiProj_noBckg_pT%d",iPt));
    fHistCorrectedProjection[iPt]->Write();

    yields[0][iPt]=fHistCorrectedProjection[iPt]->IntegralAndError(fHistCorrectedProjection[iPt]->FindBin(-0.9),fHistCorrectedProjection[iPt]->FindBin(0.9),yields_err[0][iPt],"width");
    yields[1][iPt]=fHistCorrectedProjection[iPt]->IntegralAndError(fHistCorrectedProjection[iPt]->FindBin(TMath::Pi()-1.4),fHistCorrectedProjection[iPt]->FindBin(TMath::Pi()+1.4),yields_err[1][iPt],"width");
    fHistNear->SetBinContent(iPt+1,yields[0][iPt]);
    fHistNear->SetBinError(iPt+1,yields_err[0][iPt]);
    fHistAway->SetBinContent(iPt+1,yields[1][iPt]);
    fHistAway->SetBinError(iPt+1,yields_err[1][iPt]);
  }

  Plotter::SetHistAxes(fHistNear,"#font[12]{p}^{assoc}_{T} (GeV/#font[12]{c})","Y");
  Plotter::SetHist(fHistNear,"",20,kRed+1,1.);

  Plotter::SetHistAxes(fHistAway,"#font[12]{p}^{assoc}_{T} (GeV/#font[12]{c})","Y");
  Plotter::SetHist(fHistAway,"",20,kBlue+1,1.);

  Plotter::SetHistAxes(fHistUE,"#font[12]{p}^{assoc}_{T} (GeV/#font[12]{c})","Y");
  Plotter::SetHist(fHistUE,"",20,kGreen+1,1.);

  TCanvas *can = Plotter::CreateCanvas("c");
  gStyle->SetErrorX(0.01);
  can->GetPadSave()->SetLogy();
  fHistNear->GetYaxis()->SetRangeUser(0.5*fHistUE->GetMinimum(),1.5*fHistNear->GetMaximum());
  fHistNear->DrawCopy("");
  fHistAway->DrawCopy("p same");
  fHistUE->DrawCopy("p same");

  fHistNear->Write();
  fHistAway->Write();
  fHistUE->Write();

  TPaveText *pave = new TPaveText();
  Plotter::SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.95, 0.65,0.95);
  pave->AddText("ALICE, Work in Progress");
  pave->AddText("pp, 13.6 TeV");
  pave->AddText(Form("h-%s",finalNames[part].Data()));
  pave->AddText(Form("3 < #font[12]{p}^{trigg}_{T} < 20 GeV/#font[12]{c}"));
  pave->AddText("|#Delta#eta| < 1");
  pave->Draw("same");

  TLegend *leg = Plotter::CreateLegend(0.15, 0.45, 0.15, 0.45,0.05);
  leg->AddEntry(fHistNear,"Near-side, |#Delta#varphi|<0.9","pl");
  leg->AddEntry(fHistAway,"Away-side, |#Delta#varphi-#pi|<1.4","pl");
  leg->AddEntry(fHistUE,"Underlying event #times 1/50","pl");
  leg->Draw("same");

  can->SaveAs(Form("../Plots/Yield_%s.pdf",nameSave[part].Data()));
  fFileNew->Close();
}
//____________________________________________________________________
TH1F * GetBachgroundHist(TH1F* hist){

  Int_t bins[8]={1,2,3,4,33,34,35,36};
  Double_t value =0;
  Double_t err=0;

  for (size_t i = 0; i < 8; i++) {
    value+=hist->GetBinContent(bins[i]);
    err+=TMath::Power(hist->GetBinError(bins[i]),2);
  }
  err=TMath::Sqrt(err);
  value=value/8;
  err=err/8;

  TH1F* bckg = (TH1F*) hist->Clone();
  for (size_t iPhi = 1; iPhi < hist->GetXaxis()->GetNbins()+1; iPhi++) {
    bckg->SetBinContent(iPhi,value);
    bckg->SetBinError(iPhi,err);
  }
  return bckg;

}
