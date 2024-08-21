#include "Plotter.h"

void PlotDeltaEtaProjections(Int_t part=0,Int_t ptTriggBin = 0, Int_t multClass =0){

  Int_t particleType =part;
  if(part==2)particleType=3;
  if(part==3)particleType=5;
  if(part==4)particleType=7;

  TString name[]={"K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus","Pion"};
  TString finalNames[]={"K_{S}^{0}","(#Lambda+#bar{#Lambda})","(#Xi^{+}+#Xi^{-})","(#Omega^{+}+#Omega^{-})","#pi^{+}+#pi^{-}"};
  TString nameSave[]={"K0s","Lam","Xi","Omega","Pion"};

  TString multiplicityNames[]={"minBias","0_1Mult","1_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_70Mult","70_100Mult"};
  TString multiplicityPave[]={"MB","0-1%","1-10%","10-20%","20-30%","30-40%","40-50%","50-70%","70-100%"};

  TFile * fFile[3];
  TFile * fFile2[3];
  TString InvMassRanges[] = {"Signal", "LeftBg", "RightBg"};
  // TString InvMassRanges[] = {"p"};
  Int_t nFile = sizeof(InvMassRanges) / sizeof(TString);
  Double_t ptTriggBins[]={2.,4.,6.,10.,100.};
  const Int_t nPtBins = 9;
  Double_t ptBins[]={0.,1.,2.,3.,4.,6.,8.,10.,12,15};

  TH2F * fHistCorrected[nPtBins][3];
  TH1F * fHistCorrectedProjection[nPtBins];
  TCanvas * can[nPtBins];
  TPaveText *pave[nPtBins];
  TLine * leftLine[nPtBins];
  TLine * rightLine[nPtBins];



  for (Int_t i = 0; i < nFile; i++) {
    fFile[i] = new TFile(Form("../data/MixCorrected/%s/MixCorrected_data_longTrain_%s_%s_%s.root",name[particleType].Data(),InvMassRanges[i].Data(),name[particleType].Data(),multiplicityNames[multClass].Data()));
    // fFile[i] = new TFile(Form("../FromKai/RootFiles/Shoving/MixCorrected/MixCorrectedMC_3.root"));
    if(part>0&&part<4) fFile2[i] = new TFile(Form("../data/AutoCorrelations/MixCorrected_noAutoCorr/%s/MixCorrected_%s_%s_%s.root",name[particleType+1].Data(),InvMassRanges[i].Data(),name[particleType+1].Data(),multiplicityNames[multClass].Data()));
  }

  for (Int_t iPt = 0; iPt < nPtBins; iPt++) {
    for (Int_t iFile = 0; iFile < nFile; iFile++) { // side-band subtraction
      fHistCorrected[iPt][iFile] = (TH2F*)fFile[iFile]->Get(Form("fHistCorrected_%s_pt%d_ptTrigg%d",name[particleType].Data(),iPt,ptTriggBin));
      if(part>0&&part<4)fHistCorrected[iPt][iFile]->Add ((TH2F*)fFile2[iFile]->Get(Form("fHistCorrected_%s_pt%d_ptTrigg%d",name[particleType+1].Data(),iPt,ptTriggBin)));
      if(iFile>0) fHistCorrected[iPt][0]->Add(fHistCorrected[iPt][iFile],-1);
    }

    fHistCorrectedProjection[iPt] = (TH1F *) fHistCorrected[iPt][0]->ProjectionX();
    fHistCorrectedProjection[iPt]->SetName(Form("proj%i",iPt));
    fHistCorrectedProjection[iPt]->Scale(1./fHistCorrected[iPt][0]->GetXaxis()->GetBinWidth(2));

    Plotter::SetHistAxes(fHistCorrectedProjection[iPt],"#Delta#eta","d#font[52]{N}/d#Delta#eta");
    Plotter::SetHist(fHistCorrectedProjection[iPt],"",20,kBlack,1.);

    can[iPt] = Plotter::CreateCanvas(Form("can%i",iPt));
    fHistCorrectedProjection[iPt]->GetYaxis()->SetRangeUser(0.96*fHistCorrectedProjection[iPt]->GetMinimum(),1.1*fHistCorrectedProjection[iPt]->GetMaximum());
    fHistCorrectedProjection[iPt]->Draw();

    leftLine[iPt] = new TLine(-1.1,fHistCorrectedProjection[iPt]->GetMinimum(),-1.1,fHistCorrectedProjection[iPt]->GetMaximum());
    // leftLine[iPt]->SetName(Form("left%i",iPt));
    leftLine[iPt] ->SetLineColor(kMagenta);
    leftLine[iPt] ->Draw("same");

    rightLine[iPt] = new TLine(1.1,fHistCorrectedProjection[iPt]->GetMinimum(),1.1,fHistCorrectedProjection[iPt]->GetMaximum());
    // rightLine[iPt]->SetName(Form("right%i",iPt));
    rightLine[iPt] ->SetLineColor(kMagenta);
    rightLine[iPt] ->Draw("same");

    pave[iPt] = new TPaveText();
    Plotter::SetPaveText(pave[iPt],42,0.05, 0, 0,33,0,0.55,0.95, 0.6,0.95);
    pave[iPt]->AddText("ALICE, Work in Progress");
    pave[iPt]->AddText("pp, 13.6 TeV");
    pave[iPt]->AddText(Form("%s",multiplicityPave[multClass].Data()));
    pave[iPt]->AddText(Form("h-%s",finalNames[part].Data()));
    pave[iPt]->AddText(Form("%g < #font[52]{p}^{trigg}_{T} < %g GeV/#font[52]{c}",ptTriggBins[ptTriggBin],ptTriggBins[ptTriggBin+1]));
    pave[iPt]->AddText(Form("%g < #font[52]{p}^{assoc}_{T} < %g GeV/#font[52]{c}",ptBins[iPt],ptBins[iPt+1]));
    pave[iPt]->Draw("same");

    can[iPt]->SaveAs(Form("../Plots/DeltaEtaProj/DeltaEta_%s_%s_%i_%i.pdf",nameSave[part].Data(),multiplicityNames[multClass].Data(),iPt,ptTriggBin));


  }
}
