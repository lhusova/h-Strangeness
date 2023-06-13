#include "Plotter.h"

void IntegratedYieldVsMult(Int_t iPart = 1){

  if(iPart>3) {
    cout<< "Error!!! Assoc particle not defined !!!" << endl;
    return;
  }

  TString particleName[]={"K0s","Lam","Xi","Omega"};
  Float_t particleMass[]={0.497,1.115,1.321,1.672};
  TString finalNames[]={"K_{S}^{0}","(#Lambda+#bar{#Lambda})","(#Xi^{+}+#Xi^{-})","(#Omega^{+}+#Omega^{-})"};
  TString multiplicityNames[]={"0_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_60Mult","60_70Mult","70_80Mult","80_90Mult","90_100Mult"};
  TString multiplicityPave[]={"0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"};
  TString namesRegions[3] = {"fHistNear","fHistAway","fHistUE"};
  TString paveRegions[3] = {"Near-side","Away-side","Underlying event"};
  Color_t colRegions[3] = {kRed+1,kBlue+1,kGreen+2};

  TFile * file[10];
  TH1F* fProj[10][3];
  TCanvas *can[10][3];
  TPaveText *pave[10][3];
  TLegend *leg = Plotter::CreateLegend(0.65, 0.95, 0.25, 0.5,0.05);
  Double_t yield, yielderr, err;

  const char *labels[11] = {"90-100%","80-90%","70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","0-10%","MB"};
  TH1F *histYield[3];
  histYield[0] = new TH1F("histYieldNear","",10,0,10);
  histYield[1] = new TH1F("histYieldAway","",10,0,10);
  histYield[2] = new TH1F("histYieldUE","",10,0,10);

  TF1 * boltzmann = new TF1("boltzmann","x*[0]*TMath::Sqrt([2]*[2]+x*x)*TMath::Sqrt(-TMath::Sqrt([2]*[2]+x*x)/[1])",0,10);
  boltzmann->SetParameter(0,0.1);
  boltzmann->SetParameter(1,70);
  boltzmann->FixParameter(2,particleMass[iPart]);

  TF1 * fermiDir = new TF1("fermiDir","x*[0]/(TMath::Exp(TMath::Sqrt([2]*[2]+x*x)/[1])+1)",0,10);
  fermiDir->SetParameter(0,1);
  fermiDir->SetParameter(1,160);
  fermiDir->FixParameter(2,particleMass[iPart]);

  TF1 * levy = new TF1("levy","[3]*x/([0]*[1])*(([0]-1)*([0]-2))/([0]*[1]+[2]*([0]-2))*TMath::Power(1+(TMath::Sqrt([2]*[2]+x*x)-[2])/([0]*[1]),-[0])",0,10);
  levy->SetParameter(0,7);
  levy->SetParameter(1,0.8);
  levy->FixParameter(2,particleMass[iPart]);
  levy->SetParameter(3,0.05);

  TF1 * mT = new TF1("mT","x*[0]*TMath::Exp(-TMath::Sqrt([2]*[2]+x*x)/[1])",0,10);
  mT->SetParameter(0,0.01);
  mT->SetParameter(1,160);
  mT->FixParameter(2,particleMass[iPart]);
  for (Int_t iMult = 0; iMult < 10; iMult++) {
    file[iMult] = new TFile(Form("../data/Yields_%s_%s.root",particleName[iPart].Data(),multiplicityNames[iMult].Data()));
    for (Int_t iReg = 0; iReg < 3; iReg++) {
      fProj[iMult][iReg] = (TH1F *) file[iMult]->Get(namesRegions[iReg].Data());
      can[iMult][iReg] = Plotter::CreateCanvas(Form("can%i%i",iMult,iReg));
      if(iReg==2)fProj[iMult][iReg]->Scale(50);
      fProj[iMult][iReg]->Fit(boltzmann,"me");
      fProj[iMult][iReg]->Fit(levy,"me");
      fProj[iMult][iReg]->Fit(fermiDir,"me");
      fProj[iMult][iReg]->Fit(mT,"me");
      levy->SetLineColor(kBlue);
      boltzmann->SetLineColor(kMagenta);
      mT->SetLineColor(kGreen+2);
      fProj[iMult][iReg]->GetYaxis()->SetMaxDigits(2);
      fProj[iMult][iReg]->GetXaxis()->SetRangeUser(0,12);
      fProj[iMult][iReg]->DrawCopy();
      fermiDir->DrawCopy("same");
      levy->DrawCopy("same");
      boltzmann->DrawCopy("same");
      mT->DrawCopy("same");
      if(iMult==0&&iReg==0){
        leg->AddEntry(fermiDir,"Fermi-Dirac","l");
        leg->AddEntry(levy,"Levy-Tsallis","l");
        leg->AddEntry(boltzmann,"Boltzmann","l");
        leg->AddEntry(mT,"m_{T} scaling","l");
      }
      leg->Draw();

      pave[iMult][iReg] = new TPaveText();
      Plotter::SetPaveText(pave[iMult][iReg],42,0.05, 0, 0,33,0,0.55,0.97, 0.6,0.95);
      pave[iMult][iReg]->AddText("ALICE, Work in Progress");
      pave[iMult][iReg]->AddText("pp, 13.6 TeV");
      pave[iMult][iReg]->AddText(Form("%s",multiplicityPave[iMult].Data()));
      pave[iMult][iReg]->AddText(Form("h-%s",finalNames[iPart].Data()));
      pave[iMult][iReg]->AddText(Form("3 < #font[12]{p}^{trigg}_{T} < 20 GeV/#font[12]{c}"));
      pave[iMult][iReg]->AddText("|#Delta#eta| < 1");
      pave[iMult][iReg]->AddText(paveRegions[iReg].Data());
      pave[iMult][iReg]->Draw("same");

      yield=fProj[iMult][iReg]->Integral("width");
      yield+=fermiDir->Integral(0,0.5);
      yield+=fermiDir->Integral(10,100);

      histYield[iReg]->SetBinContent(10-iMult,yield);
      histYield[iReg]->SetBinError(10-iMult,0.001);
      histYield[iReg]->GetXaxis()->SetBinLabel(iMult+1,labels[iMult]);
      //TODO: error propagation!!!

      can[iMult][iReg]->SaveAs(Form("../Plots/FitPt/%s/Ptspectrum_%s_%s.pdf",particleName[iPart].Data(),namesRegions[iReg].Data(),multiplicityNames[iMult].Data()));
    }
  }
  TLegend *legYield = Plotter::CreateLegend(0.15, 0.35, 0.65, 0.8,0.05);
  TCanvas * canYield = Plotter::CreateCanvas(Form("canY"));
  Plotter::SetHistAxes(histYield[0],"","Y");
  histYield[0]->GetYaxis()->SetRangeUser(0.002,0.299);
  for (size_t i = 0; i < 3; i++) {
    Plotter::SetHist(histYield[i],"",20,colRegions[i],1.);
    if(i==0)histYield[0]->DrawCopy();
    else histYield[i]->DrawCopy("same");
    legYield->AddEntry(histYield[i],paveRegions[i].Data(),"pl");
  }

  TPaveText * paveYields = new TPaveText();
  Plotter::SetPaveText(paveYields,42,0.05, 0, 0,33,0,0.55,0.97, 0.6,0.95);
  paveYields->AddText("ALICE, Work in Progress");
  paveYields->AddText("pp, 13.6 TeV");
  paveYields->AddText(Form("h-%s",finalNames[iPart].Data()));
  paveYields->AddText(Form("3 < #font[12]{p}^{trigg}_{T} < 20 GeV/#font[12]{c}"));
  paveYields->AddText("|#Delta#eta| < 1");
  paveYields->Draw("same");
  legYield->Draw();



}
