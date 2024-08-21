#include "Plotter.h"

void IntegratedYieldVsMult(Int_t iPart = 0){

  if(iPart>4) {
    cout<< "Error!!! Assoc particle not defined !!!" << endl;
    return;
  }

  TString particleName[]={"K0s","Lam","Xi","Omega","Pion"};
  Float_t particleMass[]={0.497,1.115,1.321,1.672,0.1396};
  TString finalNames[]={"K_{S}^{0}","(#Lambda+#bar{#Lambda})","(#Xi^{+}+#Xi^{-})","(#Omega^{+}+#Omega^{-})","#pi^{+}+#pi^{-}"};
  // TString multiplicityNames[]={"0_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_60Mult","60_70Mult","70_80Mult","80_90Mult","90_100Mult"};
  // TString multiplicityPave[]={"0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"};
  TString multiplicityNames[]={"minBias","MB","1_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_70Mult","70_100Mult"};
  TString multiplicityPave[]={"MB","0-1%","1-10%","10-20%","20-30%","30-40%","40-50%","50-70%","70-100%"};
  TString namesRegions[3] = {"fHistNear","fHistAway","fHistUE"};
  TString paveRegions[3] = {"Near-side","Away-side","Underlying event"};
  Color_t colRegions[3] = {kRed+1,kBlue+1,kGreen+2};
  Int_t markers[4] = {20,21,29,33};
  Double_t ptTriggBins[]={2.,4.,6.,10.,100};
  // Double_t ptTriggBins[]={0.,1.,2.,3.,100};

  TFile * file[10][4];
  TH1F* fProj[10][4][3];
  TCanvas *can[10][4][3];
  TPaveText *pave[10][4][3];
  TLegend *leg = Plotter::CreateLegend(0.65, 0.95, 0.25, 0.5,0.05);
  Double_t yield, yielderr, err;

  // const char *labels[11] = {"90-100%","80-90%","70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","0-10%","MB"};
  const char *labels[9] = {"70-100%","50-70%","40-50%","30-40%","20-30%","10-20%","1-10%","0-1%","MB"};
  TH1F *histYield[4][3];
  for (Int_t iPtTrigg = 0; iPtTrigg < 4; iPtTrigg++) {
    histYield[iPtTrigg][0] = new TH1F(Form("histYieldNear%d",iPtTrigg),"",9,0,9);
    histYield[iPtTrigg][1] = new TH1F(Form("histYieldAway%d",iPtTrigg),"",9,0,9);
    histYield[iPtTrigg][2] = new TH1F(Form("histYieldUE%d",iPtTrigg),"",9,0,9);
  }


  TF1 * boltzmann = new TF1("boltzmann","x*[0]*TMath::Sqrt([2]*[2]+x*x)*TMath::Exp(-TMath::Sqrt([2]*[2]+x*x)/[1])",0,6);
  if(iPart==0){
    boltzmann->SetParameter(0,0.04);
    boltzmann->SetParameter(1,0.1);
  }
  if(iPart==4){
    boltzmann->SetParameter(0,0.1);
    boltzmann->SetParameter(1,1.5);
  }
  boltzmann->FixParameter(2,particleMass[iPart]);

  TF1 * fermiDir = new TF1("fermiDir","x*[0]/(TMath::Exp(TMath::Sqrt([2]*[2]+x*x)/[1])+1)",0,6);
  fermiDir->SetParameter(0,1);
  fermiDir->SetParameter(1,160);
  fermiDir->FixParameter(2,particleMass[iPart]);

  TF1 * levy = new TF1("levy","[3]*x/([0]*[1])*(([0]-1)*([0]-2))/([0]*[1]+[2]*([0]-2))*TMath::Power(1+(TMath::Sqrt([2]*[2]+x*x)-[2])/([0]*[1]),-[0])",0,6);
  levy->SetParameter(0,7);
  levy->SetParameter(1,0.8);
  levy->FixParameter(2,particleMass[iPart]);
  levy->SetParameter(3,0.03);

  TF1 * mT = new TF1("mT","x*[0]*TMath::Exp(-TMath::Sqrt([2]*[2]+x*x)/[1])",0,6);
  mT->SetParameter(0,0.01);
  mT->SetParameter(1,0.7);
  mT->SetParLimits(1,0.1,10);
  mT->FixParameter(2,particleMass[iPart]);

  for (Int_t iMult = 0; iMult < 9; iMult++) {
    for (Int_t iPtTrigg = 2; iPtTrigg < 3; iPtTrigg++) {
      if(iMult<2) file[iMult][iPtTrigg] = new TFile(Form("../data/Yields/%s/Yields_longTrain_%s_%s_fullrangePeak_11_flat_ptTrigg%d.root",particleName[iPart].Data(),particleName[iPart].Data(),multiplicityNames[iMult].Data(),iPtTrigg));
      else file[iMult][iPtTrigg] = new TFile(Form("../data/Yields/%s/K0_Yields/Yields_longTrain_%s_%s_fullrangePeak_11_flat_ptTrigg%d.root",particleName[iPart].Data(),particleName[iPart].Data(),multiplicityNames[iMult].Data(),iPtTrigg));
      for (Int_t iReg = 2; iReg < 3; iReg++) {
        fProj[iMult][iPtTrigg][iReg] = (TH1F *) file[iMult][iPtTrigg]->Get(namesRegions[iReg].Data());
        can[iMult][iPtTrigg][iReg] = Plotter::CreateCanvas(Form("can%i%i%i",iMult,iPtTrigg,iReg));
        if(iReg==2)fProj[iMult][iPtTrigg][iReg]->Scale(50);
        // fProj[iMult][iPtTrigg][iReg]->Fit(boltzmann,"mer");
        // fProj[iMult][iPtTrigg][iReg]->Fit(levy,"mer");
        // fProj[iMult][iPtTrigg][iReg]->Fit(fermiDir,"mer");
        // fProj[iMult][iPtTrigg][iReg]->Fit(mT,"mer");
        levy->SetLineColor(kBlue);
        boltzmann->SetLineColor(kMagenta);
        mT->SetLineColor(kGreen+2);
        fProj[iMult][iPtTrigg][iReg]->GetYaxis()->SetMaxDigits(2);
        fProj[iMult][iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0,fProj[iMult][iPtTrigg][iReg]->GetMaximum()*1.2);
        fProj[iMult][iPtTrigg][iReg]->GetXaxis()->SetRangeUser(0,15);
        fProj[iMult][iPtTrigg][iReg]->DrawCopy();
        // fermiDir->DrawCopy("same");
        // levy->DrawCopy("same");
        // boltzmann->DrawCopy("same");
        // mT->DrawCopy("same");
        if(iMult==0&&iReg==0&&iPtTrigg==0){
          leg->AddEntry(fermiDir,"Fermi-Dirac","l");
          leg->AddEntry(levy,"Levy-Tsallis","l");
          leg->AddEntry(boltzmann,"Boltzmann","l");
          leg->AddEntry(mT,"m_{T} scaling","l");
        }
        // leg->Draw();

        pave[iMult][iPtTrigg][iReg] = new TPaveText();
        Plotter::SetPaveText(pave[iMult][iPtTrigg][iReg],42,0.05, 0, 0,33,0,0.55,0.97, 0.6,0.95);
        pave[iMult][iPtTrigg][iReg]->AddText("ALICE, Work in Progress");
        pave[iMult][iPtTrigg][iReg]->AddText("pp, 13.6 TeV");
        pave[iMult][iPtTrigg][iReg]->AddText(Form("%s",multiplicityPave[iMult].Data()));
        pave[iMult][iPtTrigg][iReg]->AddText(Form("h-%s",finalNames[iPart].Data()));
        pave[iMult][iPtTrigg][iReg]->AddText(Form("%g < #font[52]{p}^{trigg}_{T} < %g GeV/#font[52]{c}",ptTriggBins[iPtTrigg],ptTriggBins[iPtTrigg+1]));
        pave[iMult][iPtTrigg][iReg]->AddText("|#Delta#eta| < 1");
        pave[iMult][iPtTrigg][iReg]->AddText(paveRegions[iReg].Data());
        pave[iMult][iPtTrigg][iReg]->Draw("same");

        // yield=fProj[iMult][iPtTrigg][iReg]->Integral(0,8,"width");
        yield=0.;
        for (size_t i = 1; i < 7; i++) {
          yield+=fProj[iMult][iPtTrigg][iReg]->GetBinContent(i)*fProj[iMult][iPtTrigg][iReg]->GetBinWidth(i);
        }
        cout << iMult << "___" << iPtTrigg << "_____" << yield << endl;
        // yield+=fermiDir->Integral(0,0.5);
        // yield+=fermiDir->Integral(6,100);

        histYield[iPtTrigg][iReg]->SetBinContent(9-iMult,yield);
        histYield[iPtTrigg][iReg]->SetBinError(9-iMult,0.001);
        histYield[iPtTrigg][iReg]->GetXaxis()->SetBinLabel(iMult+1,labels[iMult]);
        //TODO: error propagation!!!

        can[iMult][iPtTrigg][iReg]->SaveAs(Form("../Plots/FitPt/%s/Uncorrected/Ptspectrum_%s_%s_ptTrigg%d.pdf",particleName[iPart].Data(),namesRegions[iReg].Data(),multiplicityNames[iMult].Data(),iPtTrigg));
      }
    }
  }
  TLegend *legYield = Plotter::CreateLegend(0.15, 0.35, 0.65, 0.8,0.05);
  TLegend *legPtTrigg = Plotter::CreateLegend(0.35, 0.65, 0.65, 0.8,0.05);
  TCanvas * canYield = Plotter::CreateCanvas(Form("canY"));

  // histYield[2][0]->GetYaxis()->SetRangeUser(-0.002,0.199);
  for (size_t i = 2; i < 3; i++) {
    for (Int_t iPtTrigg = 2; iPtTrigg < 3; iPtTrigg++)
    {
      Plotter::SetHist(histYield[iPtTrigg][i],"",markers[iPtTrigg],colRegions[i],1.);
      if(i==2&&iPtTrigg==2){
        Plotter::SetHistAxes(histYield[iPtTrigg][i],"","Y");
        histYield[iPtTrigg][i]->DrawCopy();
      }
      else histYield[iPtTrigg][i]->DrawCopy("same");
      if(i==2) legPtTrigg->AddEntry(histYield[iPtTrigg][i],Form("%g < #font[52]{p}^{trigg}_{T} < %g GeV/#font[52]{c}",ptTriggBins[iPtTrigg],ptTriggBins[iPtTrigg+1]),"p");
    }
    legYield->AddEntry(histYield[0][i],paveRegions[i].Data(),"pl");
  }

  TPaveText * paveYields = new TPaveText();
  Plotter::SetPaveText(paveYields,42,0.05, 0, 0,33,0,0.55,0.97, 0.6,0.95);
  paveYields->AddText("ALICE, Work in Progress");
  paveYields->AddText("pp, 13.6 TeV");
  paveYields->AddText(Form("h-%s",finalNames[iPart].Data()));
  paveYields->AddText("|#Delta#eta| < 1.1");
  paveYields->AddText("0 < #font[52]{p}^{assoc}_{T} < 8 GeV/#font[52]{c}");
  paveYields->AddText("Underlying Event");
  // paveYields->AddText("Away-side");
  paveYields->Draw("same");
  // legYield->Draw();
  legPtTrigg->Draw();



}
