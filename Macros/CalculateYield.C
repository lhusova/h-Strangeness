#include "Plotter.h"

TH1F * GetBackgroundHist(TH1F* hist);
TH1F * GetBackgroundLongRange(TH1F* hist,TH2F * hist2d);
TH1F * GetBackgroundFlow(TH1F* hist,Int_t part, Double_t meanPtTrigg, Double_t ptAssoc, TFile *  fileFlow, TFile *  fileFlow_v3,Bool_t include_v3);

void CalculateYield(Int_t part=4,Int_t multClass =1, Int_t bckg = 0){

  Int_t particleType =part;
  if(part==2)particleType=3;
  if(part==3)particleType=5;
  if(part==4)particleType=7;
  TString name[]={"K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus","Pion"};
  TString finalNames[]={"K_{S}^{0}","(#Lambda+#bar{#Lambda})","(#Xi^{+}+#Xi^{-})","(#Omega^{+}+#Omega^{-})","#pi^{+}+#pi^{-}"};
  TString nameSave[]={"K0s","Lam","Xi","Omega","Pion"};
  // TString multiplicityNames[]={"minBias","0_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_60Mult","60_70Mult","70_80Mult","80_90Mult","90_100Mult"};
  // TString multiplicityPave[]={"MB","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"};
  TString multiplicityNames[]={"minBias","0_1Mult","1_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_70Mult","70_100Mult"};
  TString multiplicityPave[]={"MB","0-1%","1-10%","10-20%","20-30%","30-40%","40-50%","50-70%","70-100%"};
  TString nameBackground[] = {"flat","long_range","flow_modulation_v2","flow_modulation_v3"};

  TFile *fileFlow = new TFile("../data/Flow/prapared_flow_v2.root");
  TFile *fileFlow_v3 = new TFile("../data/Flow/prapared_flow_v3.root");

  TFile * fFile[3];
  TFile * fFile2[3];
  // TString InvMassRanges[] = {"Signal", "LeftBg", "RightBg"};
  TString InvMassRanges[] = {"p"};
  Int_t nFile = sizeof(InvMassRanges) / sizeof(TString);

  for (Int_t i = 0; i < nFile; i++) {
    fFile[i] = new TFile(Form("../data/MixCorrected/%s/MixCorrected_%s_%s_%s.root",name[particleType].Data(),InvMassRanges[i].Data(),name[particleType].Data(),multiplicityNames[multClass].Data()));
    if(part>0&&part<4) fFile2[i] = new TFile(Form("../data/MixCorrected/%s/MixCorrected_%s_%s_%s.root",name[particleType+1].Data(),InvMassRanges[i].Data(),name[particleType+1].Data(),multiplicityNames[multClass].Data()));
  }
  TFile * fFileTrigger = new TFile("../data/AnalysisResults_Hyperloop_24_06_Pion.root");
  TH2F* histTriggers;
  if(part<2) histTriggers = (TH2F*) fFileTrigger->Get("correlate-strangeness/sameEvent/TriggerParticlesV0");
  else if(part<4)histTriggers = (TH2F*) fFileTrigger->Get("correlate-strangeness/sameEvent/TriggerParticlesCascade");
  else histTriggers = (TH2F*) fFileTrigger->Get("correlate-strangeness/sameEvent/TriggerParticlesPion");

  if(multClass>0)histTriggers->GetYaxis()->SetRange(multClass,multClass);
  TH1F * hist1DTriggers = (TH1F *) histTriggers->ProjectionX();
  Float_t nTrigg= hist1DTriggers->Integral();

  const Int_t nPtBins = 15;
  // Double_t ptBins[nPtBins]={0.75,1.25,1.75,2.5,3.5,5.,8.};
  // Double_t ptBins[nPtBins+1]={0.5,1.,1.5,2.,3.,4.,6.,10.};
  Double_t ptBins[]={0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,5.,6.,7.,8.,10.,12,15};
  Double_t ptBins_err[nPtBins]={0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
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

  TFile * fFileNew = TFile::Open (Form("../data/Yields/%s/Yields_%s_%s_fullrangePeak_11_%s.root",nameSave[part].Data(),nameSave[part].Data(),multiplicityNames[multClass].Data(),nameBackground[bckg].Data()),"RECREATE");

  for (Int_t iPt = 0; iPt < nPtBins; iPt++) {
    for (Int_t iFile = 0; iFile < nFile; iFile++) { // side-band subtraction
      fHistCorrected[iPt][iFile] = (TH2F*)fFile[iFile]->Get(Form("fHistCorrected_%s_pt%d",name[particleType].Data(),iPt));
      if(part>0&&part<4)fHistCorrected[iPt][iFile]->Add ((TH2F*)fFile2[iFile]->Get(Form("fHistCorrected_%s_pt%d",name[particleType+1].Data(),iPt)));
      if(iFile>0) fHistCorrected[iPt][0]->Add(fHistCorrected[iPt][iFile],-1);
    }
    if(bckg==1&&iPt==0){
      fHistUE->SetBinContent(iPt+1,0.001/(ptBins[iPt+1]-ptBins[iPt]));
      fHistNear->SetBinContent(iPt+1,0.001/(ptBins[iPt+1]-ptBins[iPt]));
      fHistAway->SetBinContent(iPt+1,0.001/(ptBins[iPt+1]-ptBins[iPt]));
      continue;
    }
    fHistCorrected[iPt][0]->RebinY(2);
    fHistCorrected[iPt][0]->Scale(1./2);
    fHistCorrected[iPt][0]->Write();
    fHistCorrected[iPt][0]->GetXaxis()->SetRangeUser(-1.1,1.1);

    fHistCorrectedProjection[iPt] = (TH1F *) fHistCorrected[iPt][0]->ProjectionY();
    fHistCorrectedProjection[iPt]->Scale(fHistCorrected[iPt][0]->GetXaxis()->GetBinWidth(2));
    fHistCorrectedProjection[iPt]->Scale(1./nTrigg);

    Plotter::SetHistAxes(fHistCorrectedProjection[iPt],"#Delta#varphi","#frac{1}{N_{trigg}} #frac{dN}{d#Delta#varphi}");

    fHistCorrectedProjection[iPt]->SetName(Form("phiProj_withBckg_pT%d",iPt));
    fHistCorrectedProjection[iPt]->Write();

    if(bckg==0)fHistBckg[iPt] = GetBackgroundHist(fHistCorrectedProjection[iPt]);
    else if(bckg==1)fHistBckg[iPt] = GetBackgroundLongRange(fHistCorrectedProjection[iPt],fHistCorrected[iPt][0]);
    else if (bckg==2)fHistBckg[iPt] = GetBackgroundFlow(fHistCorrectedProjection[iPt],part,hist1DTriggers->GetMean(),(ptBins[iPt]+ptBins[iPt+1])/2,fileFlow,fileFlow_v3,kFALSE);
    else if (bckg==3)fHistBckg[iPt] = GetBackgroundFlow(fHistCorrectedProjection[iPt],part,hist1DTriggers->GetMean(),(ptBins[iPt]+ptBins[iPt+1])/2,fileFlow,fileFlow_v3,kTRUE);
    else {
      cout << "ERROR: Underlying Event not defined!!!!" << endl;
      return;
    }
    fHistBckg[iPt]->SetName(Form("underlying_ev_pT%d",iPt));

    fHistBckg[iPt]->Write();
    yield_UE[iPt]=fHistBckg[iPt]->IntegralAndError(fHistBckg[iPt]->GetXaxis()->GetFirst(),fHistBckg[iPt]->GetXaxis()->GetLast(),yield_UE_error[iPt],"width");
    yield_UE[iPt]=yield_UE[iPt]/50;
    yield_UE_error[iPt]=yield_UE_error[iPt]/50;
    cout << (ptBins[iPt+1]-ptBins[iPt]) << endl;
    fHistUE->SetBinContent(iPt+1,yield_UE[iPt]/(ptBins[iPt+1]-ptBins[iPt]));
    fHistUE->SetBinError(iPt+1,yield_UE_error[iPt]/(ptBins[iPt+1]-ptBins[iPt]));

    fHistCorrectedProjection[iPt]->Add(fHistBckg[iPt],-1);
    fHistCorrectedProjection[iPt]->SetName(Form("phiProj_noBckg_pT%d",iPt));
    fHistCorrectedProjection[iPt]->Write();

    yields[0][iPt]=fHistCorrectedProjection[iPt]->IntegralAndError(fHistCorrectedProjection[iPt]->FindBin(-TMath::Pi()/2),fHistCorrectedProjection[iPt]->FindBin(TMath::Pi()/2),yields_err[0][iPt],"width");
    yields[1][iPt]=fHistCorrectedProjection[iPt]->IntegralAndError(fHistCorrectedProjection[iPt]->FindBin(TMath::Pi()-TMath::Pi()/2),fHistCorrectedProjection[iPt]->FindBin(TMath::Pi()+TMath::Pi()/2),yields_err[1][iPt],"width");
    fHistNear->SetBinContent(iPt+1,yields[0][iPt]/(ptBins[iPt+1]-ptBins[iPt]));
    fHistNear->SetBinError(iPt+1,yields_err[0][iPt]/(ptBins[iPt+1]-ptBins[iPt]));
    fHistAway->SetBinContent(iPt+1,yields[1][iPt]/(ptBins[iPt+1]-ptBins[iPt]));
    fHistAway->SetBinError(iPt+1,yields_err[1][iPt]/(ptBins[iPt+1]-ptBins[iPt]));
  }

  Plotter::SetHistAxes(fHistNear,"#font[12]{p}^{assoc}_{T} (GeV/#font[12]{c})","1/#font[12]{N}_{Trigg}d#font[12]{N}/d#font[12]{p}_{T}");
  Plotter::SetHist(fHistNear,"",20,kRed+1,1.);

  Plotter::SetHistAxes(fHistAway,"#font[12]{p}^{assoc}_{T} (GeV/#font[12]{c})","1/#font[12]{N}_{Trigg}d#font[12]{N}/d#font[12]{p}_{T}");
  Plotter::SetHist(fHistAway,"",20,kBlue+1,1.);

  Plotter::SetHistAxes(fHistUE,"#font[12]{p}^{assoc}_{T} (GeV/#font[12]{c})","1/#font[12]{N}_{Trigg}d#font[12]{N}/d#font[12]{p}_{T}");
  Plotter::SetHist(fHistUE,"",20,kGreen+1,1.);

  TCanvas *can = Plotter::CreateCanvas("c");
  gStyle->SetErrorX(0.01);
  // can->GetPadSave()->SetLogy();
  fHistNear->GetYaxis()->SetRangeUser(0.5*fHistUE->GetMinimum(),1.5*fHistNear->GetMaximum());
  fHistNear->DrawCopy("");
  if(bckg!=1)fHistAway->DrawCopy("p same");
  fHistUE->DrawCopy("p same");

  hist1DTriggers->Write();
  fHistNear->Write();
  fHistAway->Write();
  fHistUE->Write();

  TPaveText *pave = new TPaveText();
  Plotter::SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.95, 0.65,0.95);
  pave->AddText("ALICE, Work in Progress");
  pave->AddText("pp, 13.6 TeV");
  pave->AddText(Form("%s",multiplicityPave[multClass].Data()));
  pave->AddText(Form("h-%s",finalNames[part].Data()));
  pave->AddText(Form("3 < #font[12]{p}^{trigg}_{T} < 20 GeV/#font[12]{c}"));
  pave->AddText("|#Delta#eta| < 1.1");
  pave->Draw("same");

  TLegend *leg = Plotter::CreateLegend(0.15, 0.45, 0.15, 0.45,0.05);
  leg->AddEntry(fHistNear,"Near-side, |#Delta#varphi|<#pi/2","pl");
  leg->AddEntry(fHistAway,"Away-side, |#Delta#varphi-#pi|<#pi/2","pl");
  leg->AddEntry(fHistUE,"Underlying event #times 1/50","pl");
  leg->Draw("same");

  can->SaveAs(Form("../Plots/Yields/%s/Yield_%s_%s_fullrangePeak_11_%s.pdf",nameSave[part].Data(),nameSave[part].Data(),multiplicityNames[multClass].Data(),nameBackground[bckg].Data()));
  fFileNew->Close();
}
//____________________________________________________________________
TH1F * GetBackgroundHist(TH1F* hist){

  // Int_t bins[8]={1,2,3,4,33,34,35,36};
  Int_t bins[6]={1,2,3,16,17,18};
  Double_t value =0;
  Double_t err=0;

  for (size_t i = 0; i < 6; i++) {
    value+=hist->GetBinContent(bins[i]);
    err+=TMath::Power(hist->GetBinError(bins[i]),2);
  }
  err=TMath::Sqrt(err);
  value=value/6;
  err=err/6;

  TH1F* bckg = (TH1F*) hist->Clone();
  for (size_t iPhi = 1; iPhi < hist->GetXaxis()->GetNbins()+1; iPhi++) {
    bckg->SetBinContent(iPhi,value);
    bckg->SetBinError(iPhi,err);
  }
  return bckg;

}
//____________________________________________________________________
TH1F * GetBackgroundLongRange(TH1F* hist,TH2F * hist2d){

  hist2d->GetXaxis()->SetRangeUser(-1.5,-1.1);
  TH1F * projLeft = (TH1F*)hist2d->ProjectionY();

  hist2d->GetXaxis()->SetRangeUser(1.1,1.5);
  TH1F * projRight = (TH1F*)hist2d->ProjectionY();

  TH1F * histBckg = (TH1F*) projLeft->Clone();
  histBckg->Add(projRight);

  // histBckg->Scale(hist->GetBinContent(36)/histBckg->GetBinContent(36));
   histBckg->Scale(hist->GetBinContent(17)/histBckg->GetBinContent(17));

  return histBckg;
}
//____________________________________________________________________
TH1F * GetBackgroundFlow(TH1F* hist, Int_t part, Double_t meanPtTrigg, Double_t ptAssoc, TFile *  fileFlow, TFile *  fileFlow_v3, Bool_t include_v3){

  // Int_t bins[8]={1,2,3,4,33,34,35,36};
  Int_t bins[6]={1,2,3,16,17,18};
  Double_t value =0;

  for (size_t i = 0; i < 6; i++) {
    value+=hist->GetBinContent(bins[i]);
  }
  value=value/6;

  TGraphErrors * trigg = (TGraphErrors * ) fileFlow->Get("v2_charged");
  TGraphErrors * assoc;
  if(part==0) assoc = (TGraphErrors * ) fileFlow->Get("v2_K0_Combined");
  else if(part==1) assoc = (TGraphErrors * ) fileFlow->Get("v2_Lambda_Combined");
  else if(part==4) assoc = (TGraphErrors * ) fileFlow->Get("v2_Pion_Combined");
  else {
    cout << "ERROR!! Cascade v2 not available!!! ";
    return 0x0;
  }

  TF1 * bcgkFunction;
  if(include_v3){
    TGraphErrors * trigg_v3 = (TGraphErrors * ) fileFlow_v3->Get("v3_charged");
    TGraphErrors * assoc_3;
    if(part==0) assoc_3 = (TGraphErrors * ) fileFlow_v3->Get("v3_K0_Combined");
    else if(part==1) assoc_3 = (TGraphErrors * ) fileFlow_v3->Get("v3_Lambda_Combined");
    else if(part==4) assoc_3 = (TGraphErrors * ) fileFlow_v3->Get("v3_Pion_Combined");
    else {
      cout << "ERROR!! Cascade v2 not available!!! ";
      return 0x0;
    }

    bcgkFunction = new TF1("bcgkFunction","[0]*(1+2*[1]*([2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))",-TMath::Pi()/2,3/2*TMath::Pi());
    bcgkFunction->SetParameter(0,value);

    bcgkFunction->SetParameter(1,trigg->Eval(meanPtTrigg));
    cout << assoc->Eval(ptAssoc) << endl;
    bcgkFunction->SetParameter(2,assoc->Eval(ptAssoc));
    bcgkFunction->SetParameter(3,assoc_3->Eval(ptAssoc));

  }else{
    bcgkFunction = new TF1("bcgkFunction","[0]*(1+2*[1]*[2]*TMath::Cos(2*x))",-TMath::Pi()/2,3/2*TMath::Pi());
    cout << value << endl;
    cout << trigg->Eval(meanPtTrigg) << endl;
    bcgkFunction->SetParameter(0,value);

    bcgkFunction->SetParameter(1,trigg->Eval(meanPtTrigg));
    cout << assoc->Eval(ptAssoc) << endl;
    bcgkFunction->SetParameter(2,assoc->Eval(ptAssoc));
  }


  TH1F * bckg = (TH1F *)hist->Clone();

  for (Int_t i = 0; i < bckg->GetXaxis()->GetNbins(); i++) {
    bckg->SetBinContent(i+1,bcgkFunction->Eval(bckg->GetBinCenter(i+1)));
  }
  return bckg;

}
