#include "Plotter.h"

void PerformClosureTest(Int_t part=0,Int_t ptTriggBin = 2){

  Int_t particleType =part;
  if(part==3)particleType=5;
  if(part==4)particleType=7;
  TString name[]={"K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus","Pion"};
  TString finalNames[]={"K_{S}^{0}","#Lambda","#bar{#Lambda}","(#Xi^{+}+#Xi^{-})","(#Omega^{+}+#Omega^{-})","#pi^{+}+#pi^{-}"};
  TString nameSave[]={"K0s","Lam","ALam","Xi","Omega","Pion"};
  // TString multiplicityNames[]={"minBias","0_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_60Mult","60_70Mult","70_80Mult","80_90Mult","90_100Mult"};
  // TString multiplicityPave[]={"MB","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"};
  // TString multiplicityNames[]={"minBias","0_1Mult","1_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_70Mult","70_100Mult"};
  // TString multiplicityPave[]={"MB","0-1%","1-10%","10-20%","20-30%","30-40%","40-50%","50-70%","70-100%"};
  TString nameBackground[] = {"flat","long_range","flow_modulation_v2","flow_modulation_v3"};

  TFile * fFileMCrec[3];
  // TFile * fFileMCGen[3];//
  TFile * fFileMCGen=new TFile(Form("../data/MixCorrected/%s/MixCorrected_MCgen__%s_MB.root",name[part].Data(),name[part].Data()));
  TString InvMassRanges[] = {"Signal", "LeftBg", "RightBg"};
  // TString InvMassRanges[] = {"p"};
  Int_t nFile = sizeof(InvMassRanges) / sizeof(TString);

  for (Int_t i = 0; i < nFile; i++) {
    fFileMCrec[i] = new TFile(Form("../data/MixCorrected/%s/MixCorrected_MCrec_TriggerAndAssocPrim_%s_%s_MB.root",name[part].Data(),InvMassRanges[i].Data(),name[part].Data()));
    // fFileMCGen[i] = new TFile(Form("../data/MCrec_MixCorrected/%s/MixCorrected_%s_%s_MB.root",name[part].Data(),InvMassRanges[i].Data(),name[part].Data()));

    // if(part>0&&part<4) fFile2[i] = new TFile(Form("../data/AutoCorrelations/MixCorrected_noAutoCorr/%s/MixCorrected_%s_%s_%s.root",name[particleType+1].Data(),InvMassRanges[i].Data(),name[particleType+1].Data(),multiplicityNames[multClass].Data()));
  }
  TFile * fFileTrigger = new TFile("../data/AnalysisResults_ForMCclosure_14_08.root");
  TH2F* histTriggersRec;
  if(part<3) histTriggersRec = (TH2F*) fFileTrigger->Get("correlate-strangeness_TriggerAndAssocPrim_id14337/sameEvent/TriggerParticlesV0");
  else if(part<4)histTriggersRec = (TH2F*) fFileTrigger->Get("correlate-strangeness/sameEvent/TriggerParticlesCascade");
  else histTriggersRec = (TH2F*) fFileTrigger->Get("correlate-strangeness/sameEvent/TriggerParticlesPion");

  Double_t ptTriggBins[]={2.,4.,6.,10.,100};
  // Double_t ptTriggBins[]={0.,1.,2.,3.,100};
  TH1F * hist1DTriggersRec = (TH1F *) histTriggersRec->ProjectionX();
  TCanvas *ccdddc = Plotter::CreateCanvas("dddsssd");
  hist1DTriggersRec->DrawCopy();
  Float_t nTrigg= hist1DTriggersRec->Integral(hist1DTriggersRec->FindBin(ptTriggBins[ptTriggBin]),hist1DTriggersRec->FindBin(ptTriggBins[ptTriggBin+1]-0.1));

  TH3F * histTriggersGen = (TH3F*) fFileTrigger->Get("correlate-strangeness_id14337/ClosureTest/hTrigger");
  TH1F * hist1DTriggerGen = (TH1F *) histTriggersGen->ProjectionX();
  // TCanvas *ccc = Plotter::CreateCanvas("dddd");
  hist1DTriggerGen->SetLineColor(kRed);
  hist1DTriggerGen->DrawCopy("same");
  Float_t nTriggGen= hist1DTriggerGen->Integral(hist1DTriggerGen->FindBin(ptTriggBins[ptTriggBin]),hist1DTriggerGen->FindBin(ptTriggBins[ptTriggBin+1]-0.1));


  const Int_t nPtBins = 9;
  Double_t ptBins[]={0.,1.,2.,3.,4.,6.,8.,10.,12,15};
  Double_t ptBins_err[nPtBins]={0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
  TH2F * fHistCorrected[nPtBins][3];
  TH2F * fHistGen[nPtBins];//[3];
  TH1F * fHistCorrectedProjection[nPtBins];
  TH1F * fHistGenProjection[nPtBins];
  TH1F * fHistBckg[nPtBins];
  TCanvas *canProjection[nPtBins];
  TPad * padRatio[nPtBins];
  TH1F* fHistRatio[nPtBins];
  TPaveText *pave [nPtBins];

  for (Int_t iPt = 0; iPt < nPtBins; iPt++) {
    cout << iPt << endl;
    for (Int_t iFile = 0; iFile < nFile; iFile++) { // side-band subtraction
      fHistCorrected[iPt][iFile] = (TH2F*)fFileMCrec[iFile]->Get(Form("fHistCorrected_%s_pt%d_ptTrigg%d",name[particleType].Data(),iPt,ptTriggBin));
      fHistCorrected[iPt][iFile]->SetName(Form("red2d%i%i",iPt,iFile));
      // fHistGen[iPt][iFile] = (TH2F*)fFileMCGen[iFile]->Get(Form("fHistCorrected_%s_pt%d_ptTrigg%d",name[particleType].Data(),iPt,ptTriggBin));
      // fHistGen[iPt][iFile]->SetName(Form("gen2d%i%i",iPt,iFile));
      // if(part>0&&part<4)fHistCorrected[iPt][iFile]->Add ((TH2F*)fFile2[iFile]->Get(Form("fHistCorrected_%s_pt%d_ptTrigg%d",name[particleType+1].Data(),iPt,ptTriggBin)));
      if(iFile>0) {
        fHistCorrected[iPt][0]->Add(fHistCorrected[iPt][iFile],-1);
        // fHistGen[iPt][0]->Add(fHistGen[iPt][iFile],-1);
      }
    }

    fHistCorrected[iPt][0]->GetXaxis()->SetRangeUser(-1.1,1.1);

    fHistCorrectedProjection[iPt] = (TH1F *) fHistCorrected[iPt][0]->ProjectionY();
    // fHistCorrectedProjection[iPt]->RebinX(2);
    fHistCorrectedProjection[iPt]->SetName(Form("rec%i",iPt));
    fHistCorrectedProjection[iPt]->Scale(1./fHistCorrected[iPt][0]->GetYaxis()->GetBinWidth(2));
    // fHistCorrectedProjection[iPt]->Scale(fHistCorrected[iPt][0]->GetXaxis()->GetBinWidth(2));
    fHistCorrectedProjection[iPt]->Scale(1./nTrigg);
    Plotter::SetHist(fHistCorrectedProjection[iPt],"",24,kBlue+1,1.);

    fHistCorrectedProjection[iPt]->SetName(Form("phiProj_withBckg_pT%d",iPt));

    fHistGen[iPt] = (TH2F*)fFileMCGen->Get(Form("fHistCorrected_%s_pt%d_ptTrigg%d",name[particleType].Data(),iPt,ptTriggBin));
    fHistGen[iPt]->GetXaxis()->SetRangeUser(-1.1,1.1);
    fHistGenProjection[iPt] = (TH1F *) fHistGen[iPt]->ProjectionY();
    // fHistGenProjection[iPt]->RebinX(2);
    fHistGenProjection[iPt]->SetName(Form("gen%i",iPt));
    fHistGenProjection[iPt]->Scale(1./fHistGen[iPt]->GetYaxis()->GetBinWidth(2));
    // fHistGenProjection[iPt]->Scale(fHistGen[iPt]->GetXaxis()->GetBinWidth(2));
    fHistGenProjection[iPt]->Scale(1./nTriggGen);
    Plotter::SetHist(fHistGenProjection[iPt],"",20,kRed+1,1.);
    Plotter::SetHistAxes(fHistGenProjection[iPt],"#Delta#varphi","#frac{1}{N_{trigg}} #frac{dN}{d#Delta#varphi}");

    canProjection[iPt]=Plotter::CreateCanvas(Form("c%i",iPt),700,750,true);
    padRatio[iPt] = new TPad(Form("padRatio%i",iPt),"",0.001,0.001,0.999,0.3);
    padRatio[iPt]->SetMargin(0.12,0.02,0.25,0.01);
    padRatio[iPt]->Draw();

    canProjection[iPt]->GetPadSave()->cd();
    fHistGenProjection[iPt]->DrawCopy();
    fHistCorrectedProjection[iPt]->DrawCopy("same");

    pave[iPt] = new TPaveText();
    Plotter::SetPaveText(pave[iPt],42,0.05, 0, 0,33,0,0.55,0.95, 0.72,0.95);
    pave[iPt]->AddText("ALICE, Work in Progress");
    pave[iPt]->AddText("pp, 13.6 TeV");
    pave[iPt]->AddText(Form("%g < p_{T}^{trig} < %g GeV/c",ptTriggBins[ptTriggBin],ptTriggBins[ptTriggBin+1]));
    pave[iPt]->AddText(Form("%g < p_{T}^{assoc} < %g GeV/c",ptBins[iPt],ptBins[iPt+1]));
    pave[iPt]->Draw("same");

    padRatio[iPt]->cd();
    fHistRatio[iPt] = (TH1F*)fHistCorrectedProjection[iPt]->Clone();
    // fHistRatio[iPt] = (TH1F*)fHistGenProjection[iPt]->Clone();
    fHistRatio[iPt]->SetName(Form("ratio_%i",iPt));
    fHistRatio[iPt]->Divide(fHistGenProjection[iPt]);
    // fHistRatio[iPt]->Divide(fHistCorrectedProjection[iPt]);
    Plotter::SetHistAxesSmallPad(fHistRatio[iPt],"#Delta#varphi","rec/gen");
    // Plotter::SetHistAxesSmallPad(fHistRatio[iPt],"#Delta#varphi","gen/rec");
    fHistRatio[iPt]->DrawCopy();


  }


}
