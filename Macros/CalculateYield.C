#include "Plotter.h"

TH1F * GetBachgroundHist(TH1F* hist);

void CalculateYield(Int_t part=2){

  Int_t particleType =part;
  if(part==2)particleType=3;
  if(part==3)particleType=5;
  TString name[]={"K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus"};
  TString finalNames[]={"K_{S}^{0}","(#Lambda+#bar{#Lambda})","(#Xi^{+}+#Xi^{-})","(#Omega^{+}+#Omega^{-})"};
  TString nameSave[]={"K0s","Lam","Xi","Omega"};

  TFile * fFile = new TFile(Form("../data/MixCorrected_%s.root",name[particleType].Data()));
  TFile * fFile2;
  if(part>0) fFile2 = new TFile(Form("../data/MixCorrected_%s.root",name[particleType+1].Data()));

  const Int_t nPtBins = 10;
  Double_t ptBins[nPtBins]={0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,4.25,4.75};
  Double_t ptBins_err[nPtBins]={0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  TH2F * fHistCorrected[nPtBins];
  TH1F * fHistCorrectedProjection[nPtBins];
  TH1F * fHistBckg[nPtBins];
  Double_t yield_UE[nPtBins];
  Double_t yield_UE_error[nPtBins];
  Double_t yields[2][nPtBins];
  Double_t yields_err[2][nPtBins];

  TFile * fFileNew = TFile::Open (Form("../data/Yields_%s.root",nameSave[part].Data()),"RECREATE");

  for (Int_t iPt = 0; iPt < nPtBins; iPt++) {
    fHistCorrected[iPt] = (TH2F*)fFile->Get(Form("fHistCorrected_%s_pt%d",name[particleType].Data(),iPt));
    if(part>0)fHistCorrected[iPt]->Add ((TH2F*)fFile2->Get(Form("fHistCorrected_%s_pt%d",name[particleType+1].Data(),iPt)));

    fHistCorrected[iPt]->GetXaxis()->SetRangeUser(-1,1);

    fHistCorrectedProjection[iPt] = (TH1F *) fHistCorrected[iPt]->ProjectionY();
    fHistCorrectedProjection[iPt]->Scale(fHistCorrected[iPt]->GetXaxis()->GetBinWidth(2));

    Plotter::SetHistAxes(fHistCorrectedProjection[iPt],"#Delta#varphi","#frac{dN}{d#Delta#varphi}");

    fHistCorrectedProjection[iPt]->SetName(Form("phiProj_withBckg_pT%d",iPt));
    fHistCorrectedProjection[iPt]->Write();

    fHistBckg[iPt] = GetBachgroundHist(fHistCorrectedProjection[iPt]);
    fHistBckg[iPt]->SetName(Form("underlying_ev_pT%d",iPt));

    fHistBckg[iPt]->Write();
    yield_UE[iPt]=fHistBckg[iPt]->IntegralAndError(fHistBckg[iPt]->GetXaxis()->GetFirst(),fHistBckg[iPt]->GetXaxis()->GetLast(),yield_UE_error[iPt],"width");
    yield_UE[iPt]=yield_UE[iPt]/50;
    yield_UE_error[iPt]=yield_UE_error[iPt]/50;

    fHistCorrectedProjection[iPt]->Add(fHistBckg[iPt],-1);
    fHistCorrectedProjection[iPt]->SetName(Form("phiProj_noBckg_pT%d",iPt));
    fHistCorrectedProjection[iPt]->Write();

    yields[0][iPt]=fHistCorrectedProjection[iPt]->IntegralAndError(fHistCorrectedProjection[iPt]->FindBin(-0.9),fHistCorrectedProjection[iPt]->FindBin(0.9),yields_err[0][iPt],"width");
    yields[1][iPt]=fHistCorrectedProjection[iPt]->IntegralAndError(fHistCorrectedProjection[iPt]->FindBin(TMath::Pi()-1.4),fHistCorrectedProjection[iPt]->FindBin(TMath::Pi()+1.4),yields_err[1][iPt],"width");

  }
  fFileNew->Close();

  TGraphErrors * fGraphNear = new TGraphErrors(nPtBins,ptBins,yields[0],ptBins_err,yields_err[0]);
  fGraphNear->SetName("fGraphNear");
  Plotter::SetGraphAxes(fGraphNear,"#font[12]{p}^{assoc}_{T} (GeV/#font[12]{c})","");
  Plotter::SetGraph(fGraphNear,"",20,kRed+1,1.);

  TGraphErrors * fGraphAway = new TGraphErrors(nPtBins,ptBins,yields[1],ptBins_err,yields_err[1]);
  fGraphAway->SetName("fGraphAway");
  Plotter::SetGraphAxes(fGraphAway,"#font[12]{p}^{assoc}_{T} (GeV/#font[12]{c})","");
  Plotter::SetGraph(fGraphAway,"",20,kBlue+1,1.);

  TGraphErrors * fGraphUE = new TGraphErrors(nPtBins,ptBins,yield_UE,ptBins_err,yield_UE_error);
  fGraphUE->SetName("fGraphUE");
  Plotter::SetGraphAxes(fGraphUE,"#font[12]{p}^{assoc}_{T} (GeV/#font[12]{c})","");
  Plotter::SetGraph(fGraphUE,"",20,kMagenta+1,1.);

  TCanvas *can = Plotter::CreateCanvas("c");
  can->GetPadSave()->SetLogy();
  fGraphNear->GetYaxis()->SetRangeUser(0.5*fGraphNear->GetPointY(0),1.5*fGraphNear->GetPointY(4));
  fGraphNear->Draw("ap");
  fGraphAway->Draw("p same");
  fGraphUE->Draw("p same");

  TPaveText *pave = new TPaveText();
  Plotter::SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.95, 0.65,0.95);
  pave->AddText("ALICE, Work in Progress");
  pave->AddText("pp, 13.6 TeV");
  pave->AddText(Form("h-%s",finalNames[part].Data()));
  pave->AddText("1 < #font[12]{p}^{trigg}_{T} < 3 GeV/#font[12]{c}"); //TODO check the values
  pave->AddText("|#Delta#eta| < 1");
  pave->Draw("same");

  TLegend *leg = Plotter::CreateLegend(0.25, 0.45, 0.25, 0.55,0.05);
  leg->AddEntry(fGraphNear,"Near-side, |#Delta#varphi|<0.9","pl");
  leg->AddEntry(fGraphAway,"Away-side, |#Delta#varphi-#pi|<1.4","pl");
  leg->AddEntry(fGraphUE,"Underlying event #times 1/50","pl");
  leg->Draw("same");

  can->SaveAs(Form("../Plots/Yield_%s.pdf",nameSave[part].Data()));

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
