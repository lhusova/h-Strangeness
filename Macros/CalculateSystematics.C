#include "Plotter.h"

void CalculateSystematics(TString variation = "WithPV", TString region = "UE", Int_t part=0,Int_t ptTriggBin = 2, Int_t multClass =1){

  TString nameSave[]={"K0s","Lam","Xi","Omega","Pion"};
  TString multiplicityNames[]={"minBias","0_1Mult","1_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_70Mult","70_100Mult"};
  TString finalNames[]={"K_{S}^{0}","(#Lambda+#bar{#Lambda})","(#Xi^{+}+#Xi^{-})","(#Omega^{+}+#Omega^{-})","#pi^{+}+#pi^{-}"};
  TString multiplicityPave[]={"MB","0-1%","1-10%","10-20%","20-30%","30-40%","40-50%","50-70%","70-100%"};
  Double_t ptTriggBins[]={2.,4.,6.,10.,100.};

  TFile * fFileDefault = new TFile(Form("../data/Yields/%s/K0_Yields/Yields_longTrain_%s_%s_fullrangePeak_11_flat_ptTrigg%d.root",nameSave[part].Data(),nameSave[part].Data(),multiplicityNames[multClass].Data(),ptTriggBin));
  TFile * fFileVar = new TFile(Form("../data/Yields_Systematics/%s_%s/Yields_longTrain_%s_%s_fullrangePeak_11_flat_ptTrigg%d.root",nameSave[part].Data(),variation.Data(),nameSave[part].Data(),multiplicityNames[multClass].Data(),ptTriggBin));

  TH1F * fHistDef = (TH1F*) fFileDefault->Get(Form("fHist%s",region.Data()));
  TH1F * fHistVar = (TH1F*) fFileVar->Get(Form("fHist%s",region.Data()));

  fHistDef->SetLineColor(kBlack);
  fHistDef->SetMarkerColor(kBlack);

  Plotter::SetHist(fHistVar,"",20,kRed,1.);

  TCanvas *canCompare = Plotter::CreateCanvas("padCompare");
  gPad->SetLogy();
  gStyle->SetErrorX(0);
  fHistVar->DrawCopy();
  fHistDef->DrawCopy("same");

  TPaveText *pave = new TPaveText();
  Plotter::SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.95, 0.65,0.95);
  pave->AddText("ALICE, Work in progress");
  pave->AddText("pp, 13.6 TeV");
  pave->AddText(Form("%s",multiplicityPave[multClass].Data()));
  pave->AddText(Form("h-%s",finalNames[part].Data()));
  pave->AddText(Form("%g < #font[52]{p}^{trigg}_{T} < %g GeV/#font[52]{c}",ptTriggBins[ptTriggBin],ptTriggBins[ptTriggBin+1]));
  if(region=="Near")pave->AddText("Near-Side");
  if(region=="Away")pave->AddText("Away-Side");
  if(region=="UE")pave->AddText("Underlying event");
  pave->Draw("same");
  //
  TLegend *lg = Plotter::CreateLegend(0.3,0.67,0.25,0.45,0.04);
  lg->AddEntry(fHistDef,"Default","pl");
  if(variation=="WithPV")lg->AddEntry(fHistVar,"Mixed correction in PV bins","pl");
  lg->Draw();
  //
  TCanvas *canRatio = Plotter::CreateCanvas("padRatio");
  TH1D * ratio = (TH1D*) fHistVar->Clone();
  ratio->SetName("ratio");
  ratio->GetYaxis()->SetRangeUser(0.8001,1.3999);
  ratio->GetYaxis()->SetTitle("Y^{var} / Y^{def}");
  ratio->Divide(fHistDef);
  ratio->DrawCopy();
  pave->Draw("same");
  Plotter::DrawUnity(kBlack,1,0,15);
  //
  TCanvas *canRelativeUncer = Plotter::CreateCanvas("padRelativeUncer");
  TH1D * relatUncer = (TH1D*) fHistVar->Clone();
  relatUncer->SetName("relatUncer");
  relatUncer->GetYaxis()->SetTitle("Relative Uncertainty");
  relatUncer->Add(fHistDef,-1);
  relatUncer->Divide(fHistDef);
  for (size_t i = 1; i < relatUncer->GetXaxis()->GetNbins()+1; i++) {
    relatUncer->SetBinContent(i,TMath::Abs(relatUncer->GetBinContent(i)));
    relatUncer->SetBinError(i,0.0001);
  }
  relatUncer->GetYaxis()->SetRangeUser(-0.01,0.09);
  relatUncer->DrawCopy();
  pave->Draw("same");

  TCanvas *canBarlow = Plotter::CreateCanvas("padBarlow");
  TH1D * histbarlow = (TH1D*) fHistVar->Clone();
  histbarlow->SetName("histbarlow");
  histbarlow->GetYaxis()->SetTitle("#sigma_{Barlow}");
  histbarlow->GetYaxis()->SetRangeUser(0.1,9.99);
  Double_t num, denom, uncerDef2, uncerVar2;
  for (size_t i = 1; i < relatUncer->GetXaxis()->GetNbins()+1; i++){
    num=TMath::Abs(fHistVar->GetBinContent(i)-fHistDef->GetBinContent(i));
    uncerDef2=TMath::Power(fHistDef->GetBinError(i),2);
    uncerVar2=TMath::Power(fHistVar->GetBinError(i),2);
    denom = TMath::Sqrt(TMath::Abs(uncerDef2-uncerVar2));
    histbarlow->SetBinContent(i,num/denom);
  }
  histbarlow->DrawCopy();
  pave->Draw("same");
  Plotter::DrawUnity(kBlack,1,0,15);

  // TCanvas *canRelativeUncerSmoothed = Plotter::CreateCanvas("padRelativeUncerSmooth"); // if smoothing necessary
  // TH1D * relatUncerSmooth = (TH1D*) relatUncer->Clone();
  // relatUncerSmooth->SetName("relatUncerSmooth");
  // relatUncerSmooth->GetYaxis()->SetRangeUser(-0.01,0.19);
  // relatUncerSmooth->GetYaxis()->SetTitle("Relative Uncertainty, Smoothed");
  //
  // Double_t avUncer  = 0.;
  // for (size_t i = 1; i < relatUncerSmooth->GetXaxis()->GetNbins()+1; i++) {
  //   if((species=="Proton"||(species=="Kaon"&&system=="pPb"&&(variation=="FMDcut"||variation=="Pvz"))||(species=="Pion"&&system=="pPb"&&(variation=="Bayes"))||(species=="Pion"&&system=="pp"&&(variation=="FMDcut"||variation=="Bayes"))||(species=="Kaon"&&system=="pp"&&(variation=="Pvz"))||(species=="Lambda"&&system=="pPb"&&(variation=="2Sigma"||variation=="Fit"))||(species=="Lambda"&&system=="pp"&&(variation=="LooseCuts"||variation=="FMDcut")))&&i==1) continue;
  //
  //   if(i<relatUncerSmooth->GetXaxis()->GetNbins()-1&&(species=="Proton"||(species=="Pion"&&system=="pp"&&(variation=="FMDcut"||variation=="Bayes"))||(species=="Kaon"&&system=="pp"&&(variation=="Pvz"||variation=="Base"||variation=="Bayes"))||((species=="K0s"||species=="Kaon")&&system=="pPb"&&variation=="Base")||(species=="Kaon"&&system=="pp"&&(variation=="FMDcut"||variation=="FB"))))avUncer+=relatUncer->GetBinContent(i);
  //   else if(i<relatUncerSmooth->GetXaxis()->GetNbins())avUncer+=relatUncer->GetBinContent(i);
  //
  //   if(i==1) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+1)) / 2);
  //   else if(i==2&&(species=="Proton"||(species=="Pion"&&system=="pPb"&&(variation=="Bayes"))||(species=="Kaon"&&system=="pPb"&&(variation=="FMDcut"||variation=="Pvz"))||(species=="Pion"&&system=="pp"&&(variation=="FMDcut"||variation=="Bayes"))||(species=="Kaon"&&system=="pp"&&(variation=="Pvz"))||(species=="Lambda"&&system=="pPb"&&(variation=="2Sigma"||variation=="Fit"))||(species=="Lambda"&&system=="pp"&&(variation=="LooseCuts"||variation=="FMDcut")))) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+1)) / 2);
  //   else if(i==relatUncerSmooth->GetXaxis()->GetNbins()-2&&(species=="Proton"||(species=="Kaon"&&system=="pPb"&&(variation=="FMDcut"))||(species=="Pion"&&(variation=="FMDcut"||variation=="Bayes"))||(species=="Kaon"&&system=="pp"&&(variation=="Pvz"||variation=="Base"||variation=="Bayes"))||((species=="K0s"||species=="Kaon")&&system=="pPb"&&variation=="Base")||(species=="Kaon"&&system=="pp"&&(variation=="FMDcut"||variation=="FB"))))relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i-1)) / 2);
  //   else if(i==relatUncerSmooth->GetXaxis()->GetNbins()-1&&(species=="Proton"||(species=="Kaon"&&system=="pPb"&&(variation=="FMDcut"))||(species=="Pion"&&(variation=="FMDcut"||variation=="Bayes"))||(species=="Kaon"&&system=="pp"&&(variation=="Pvz"||variation=="Base"||variation=="Bayes"))||((species=="K0s"||species=="Kaon")&&system=="pPb"&&variation=="Base")||(species=="Kaon"&&system=="pp"&&(variation=="FMDcut"||variation=="FB"))))relatUncerSmooth ->SetBinContent(i,avUncer/(relatUncerSmooth->GetXaxis()->GetNbins()-2));
  //   else if(i<relatUncerSmooth->GetXaxis()->GetNbins()-1)relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i-1) + relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+1)) / 3);
  //   else if(i==relatUncerSmooth->GetXaxis()->GetNbins()-1)relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i-1)) / 2);
  //   else if(i==relatUncerSmooth->GetXaxis()->GetNbins()&&(species=="Proton"||(species=="Kaon"&&system=="pPb"&&(variation=="FMDcut"))||(species=="Pion"&&(variation=="FMDcut"||variation=="Bayes"))||(species=="Kaon"&&system=="pp"&&(variation=="Pvz"||variation=="Base"||variation=="Bayes"))||((species=="K0s"||species=="Kaon")&&system=="pPb"&&variation=="Base")||(species=="Kaon"&&system=="pp"&&(variation=="FMDcut"||variation=="FB"))))relatUncerSmooth ->SetBinContent(i,avUncer/(relatUncerSmooth->GetXaxis()->GetNbins()-2));
  //   else relatUncerSmooth ->SetBinContent(i,avUncer/(relatUncerSmooth->GetXaxis()->GetNbins()-1));
  //
  //   relatUncerSmooth->SetBinError(i,0.0001);
  // }
  // if((species=="Lambda"&&(variation=="2Sigma"||variation=="Fit")&&system=="pPb")||(species=="Pion"&&system=="pPb"&&(variation=="Bayes"))||(species=="Pion"&&system=="pp"&&(variation=="FMDcut"||variation=="Bayes"))||(species=="Kaon"&&system=="pp"&&(variation=="Pvz"))||(species=="Kaon"&&system=="pPb"&&(variation=="FMDcut"||variation=="Pvz"))||(species=="Lambda"&&system=="pp"&&(variation=="LooseCuts"||variation=="FMDcut")))relatUncerSmooth ->SetBinContent(1,relatUncerSmooth->GetBinContent(relatUncerSmooth->GetXaxis()->GetNbins()));
  // if((species=="Lambda"&&variation=="2Sigma"&&system=="pp")||(species=="Proton"&&variation=="FMDcut"&&system=="pPb")){
  //   avUncer=0.;
  //   for (size_t i = 1; i < relatUncerSmooth->GetXaxis()->GetNbins()+1; i++) {
  //     if(species=="Proton"&&i==1) continue;
  //     if(i!=3&&i<relatUncerSmooth->GetXaxis()->GetNbins())avUncer+=relatUncer->GetBinContent(i);
  //     if(i==1) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+1)) / 2);
  //     if(i==2&&species=="Proton") relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+2)) / 2);
  //     if(i==2&&species=="Lambda") relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+2)+relatUncer->GetBinContent(i-1)) / 3);
  //     if(i==4) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i-2)+relatUncer->GetBinContent(i+1)) / 3);
  //     if(i>4&&i<relatUncerSmooth->GetXaxis()->GetNbins()-1) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i-1)+relatUncer->GetBinContent(i+1)) / 3);
  //     if(i==relatUncerSmooth->GetXaxis()->GetNbins()-1) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i-1)) / 2);
  //     if(i==relatUncerSmooth->GetXaxis()->GetNbins())relatUncerSmooth ->SetBinContent(i,avUncer/(relatUncerSmooth->GetXaxis()->GetNbins()-2));
  //   }
  //   relatUncerSmooth ->SetBinContent(3,avUncer/(relatUncerSmooth->GetXaxis()->GetNbins()-2));
  // }
  // if((species=="Kaon"&&variation=="FB"&&system=="pp")){
  //   avUncer=0.;
  //   for (size_t i = 1; i < relatUncerSmooth->GetXaxis()->GetNbins()+1; i++) {
  //     if(i==2)continue;
  //     if(i<relatUncerSmooth->GetXaxis()->GetNbins()-1)avUncer+=relatUncer->GetBinContent(i);
  //     if(i==1) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+2)) / 2);
  //     else if(i==3) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+1)+relatUncer->GetBinContent(i-2)) / 3);
  //     else if(i>3&&i<relatUncerSmooth->GetXaxis()->GetNbins()-2) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i-1)+relatUncer->GetBinContent(i+1)) / 3);
  //     else if(i==relatUncerSmooth->GetXaxis()->GetNbins()-2) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i-1)) / 2);
  //     else relatUncerSmooth ->SetBinContent(i,avUncer/(relatUncerSmooth->GetXaxis()->GetNbins()-3));
  //   }
  //   relatUncerSmooth ->SetBinContent(2,avUncer/(relatUncerSmooth->GetXaxis()->GetNbins()-3));
  // }
  // if(species=="Proton"&&variation=="Pvz"&&system=="pPb"){
  //   avUncer=0.;
  //   for (size_t i = 2; i < relatUncerSmooth->GetXaxis()->GetNbins()+1; i++) {
  //     if(i==4) continue;
  //     if(i<relatUncerSmooth->GetXaxis()->GetNbins())avUncer+=relatUncer->GetBinContent(i);
  //     if(i==2) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+1)) / 2);
  //     if(i==3) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+2)+relatUncer->GetBinContent(i-1)) / 3);
  //     if(i==5) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i-2)+relatUncer->GetBinContent(i+1)) / 3);
  //     if(i>5&&i<relatUncerSmooth->GetXaxis()->GetNbins()-1) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i-1)+relatUncer->GetBinContent(i+1)) / 3);
  //     if(i==relatUncerSmooth->GetXaxis()->GetNbins()-1) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i-1)) / 2);
  //     if(i==relatUncerSmooth->GetXaxis()->GetNbins())relatUncerSmooth ->SetBinContent(i,avUncer/(relatUncerSmooth->GetXaxis()->GetNbins()-2));
  //   }
  //   relatUncerSmooth ->SetBinContent(4,avUncer/(relatUncerSmooth->GetXaxis()->GetNbins()-2));
  // }
  // if(species=="K0s"&&(variation=="Base")&&system=="pp"){
  //   relatUncerSmooth ->SetBinContent(6,relatUncerSmooth->GetBinContent(relatUncerSmooth->GetXaxis()->GetNbins()));
  //   relatUncerSmooth ->SetBinContent(5,relatUncerSmooth->GetBinContent(relatUncerSmooth->GetXaxis()->GetNbins()));
  // }
  //
  // if(species=="Proton"&&(variation=="Bayes"||variation=="Pvz")&&system=="pp"){
  //   avUncer=0.;
  //   for (size_t i = 2; i < relatUncerSmooth->GetXaxis()->GetNbins()+1; i++) {
  //     if(i<=relatUncerSmooth->GetXaxis()->GetNbins()-3)avUncer+=relatUncer->GetBinContent(i);
  //     if(i==2) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+1)) / 2);
  //     else if(i<relatUncerSmooth->GetXaxis()->GetNbins()-3)relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+1)+relatUncer->GetBinContent(i-1)) / 3);
  //     else if(i==relatUncerSmooth->GetXaxis()->GetNbins()-3)relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i-1)) / 2);
  //     else relatUncerSmooth ->SetBinContent(i,avUncer/(relatUncer->GetXaxis()->GetNbins()-4));
  //   }
  // }
  // if(species=="Pion"&&(variation=="Base")&&system=="pPb"){
  //   avUncer=0.;
  //   for (size_t i = 1; i < relatUncerSmooth->GetXaxis()->GetNbins()+1; i++) {
  //     if(i<=relatUncerSmooth->GetXaxis()->GetNbins()-3)avUncer+=relatUncer->GetBinContent(i);
  //     if(i==1) relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+1)) / 2);
  //     else if(i<relatUncerSmooth->GetXaxis()->GetNbins()-3)relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i+1)+relatUncer->GetBinContent(i-1)) / 3);
  //     else if(i==relatUncerSmooth->GetXaxis()->GetNbins()-3)relatUncerSmooth ->SetBinContent(i, (relatUncer->GetBinContent(i)+relatUncer->GetBinContent(i-1)) / 2);
  //     else relatUncerSmooth ->SetBinContent(i,avUncer/(relatUncer->GetXaxis()->GetNbins()-3));
  //   }
  // }
  // if(species=="Proton"&&(variation=="Base")&&system=="pp")relatUncerSmooth->SetBinContent(9,relatUncerSmooth->GetBinContent(8));
  // relatUncerSmooth->DrawCopy();
  // pave->Draw("same");
  //
  // TH1D * relatUncerSmooth1;
  // if(variation=="LooseCuts"){
  //   relatUncerSmooth1 = (TH1D*) relatUncer1->Clone();
  //   relatUncerSmooth1->SetName("relatUncerSmooth1");
  //   relatUncerSmooth1->GetYaxis()->SetRangeUser(-0.01,0.19);
  //   relatUncerSmooth1->GetYaxis()->SetTitle("Relative Uncertainty, Smoothed");
  //   Double_t avUncer1  = 0.;
  //   for (size_t i = 1; i < relatUncerSmooth1->GetXaxis()->GetNbins()+1; i++) {
  //     if(i<relatUncerSmooth1->GetXaxis()->GetNbins())avUncer1+=relatUncer1->GetBinContent(i);
  //     if(i==1) relatUncerSmooth1 ->SetBinContent(i, (relatUncer1->GetBinContent(i)+relatUncer1->GetBinContent(i+1)) / 2);
  //     else if(species=="K0s"&&system=="pp"&&i==relatUncerSmooth1->GetXaxis()->GetNbins()-2)relatUncerSmooth1 ->SetBinContent(i, (relatUncer1->GetBinContent(i)+relatUncer1->GetBinContent(i-1)) / 2);
  //     else if(species=="K0s"&&system=="pp"&&i==relatUncerSmooth1->GetXaxis()->GetNbins()-1) relatUncerSmooth1 ->SetBinContent(i,avUncer1/(relatUncerSmooth1->GetXaxis()->GetNbins()-1));
  //     else if(i<relatUncerSmooth1->GetXaxis()->GetNbins()-1)relatUncerSmooth1 ->SetBinContent(i, (relatUncer1->GetBinContent(i-1) + relatUncer1->GetBinContent(i)+relatUncer1->GetBinContent(i+1)) / 3);
  //     else if(i==relatUncerSmooth1->GetXaxis()->GetNbins()-1)relatUncerSmooth1 ->SetBinContent(i, (relatUncer1->GetBinContent(i)+relatUncer1->GetBinContent(i-1)) / 2);
  //     else relatUncerSmooth1 ->SetBinContent(i,avUncer1/(relatUncerSmooth1->GetXaxis()->GetNbins()-1));
  //   }
  //   relatUncerSmooth1->DrawCopy("same");
  // }


  //
  canBarlow->SaveAs(Form("../Plots/Systematics/%s/BarlowCheck_%s_%s_%s_%i.pdf",variation.Data(),region.Data(),nameSave[part].Data(),multiplicityNames[multClass].Data(),ptTriggBin));
  canCompare->SaveAs(Form("../Plots/Systematics/%s/SystCompare_%s_%s_%s_%i.pdf",variation.Data(),region.Data(),nameSave[part].Data(),multiplicityNames[multClass].Data(),ptTriggBin));
  canRatio->SaveAs(Form("../Plots/Systematics/%s/SystRatio_%s_%s_%s_%i.pdf",variation.Data(),region.Data(),nameSave[part].Data(),multiplicityNames[multClass].Data(),ptTriggBin));
  canRelativeUncer->SaveAs(Form("../Plots/Systematics/%s/SystUncer_%s_%s_%s_%i.pdf",variation.Data(),region.Data(),nameSave[part].Data(),multiplicityNames[multClass].Data(),ptTriggBin));
  // canRelativeUncerSmoothed->SaveAs(Form("/Users/ltarasovic/alice/FlowAnalysisData/%s_Systematics/%s/TF_%s%s/SystUncerSmoothed_%s_%s.pdf",system.Data(),variation.Data(),system.Data(),pid.Data(),species.Data(),mult.Data()));
  //
  // TFile* fout = TFile::Open(Form("/Users/ltarasovic/alice/FlowAnalysisData/%s_Systematics/%s/TF_%s%s/Systematics_%s_Final_%s_withEff.root",system.Data(),variation.Data(),system.Data(),pid.Data(),species.Data(),mult.Data()),"RECREATE");
  // relatUncer->Write();
  // relatUncerSmooth->Write();
  // if(variation=="LooseCuts") {
  //   relatUncer1->Write();
  //   relatUncerSmooth1->Write();
  // }
  // fout->Close();
}
