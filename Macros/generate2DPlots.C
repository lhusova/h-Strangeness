#include "Plotter.h"

void generate2DPlots( Int_t iPart = 1){

  TString partName[] = {"Pion","K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus"};
  TString partSym[] = {"#pi^{+} + #pi^{-}","K^{0}_{S}","#Lambda","#bar{#Lambda}","#Xi^{+}","#Xi^{-}","#Omega^{+}","#Omega^{-}"};
  TFile * fFile = new TFile(Form("../data/DataThin_Corrected_MixCorrected/%s/MixCorrected_Signal_%s_0_1Mult.root",partName[iPart].Data(),partName[iPart].Data()));
  TH2D * fHist2D = (TH2D *) fFile->Get(Form("fHistMixed_%s_pt2_ptTrigg1",partName[iPart].Data()));
  // if(iPart==4||iPart==2||iPart==6){
    // TFile * fFileConj = new TFile(Form("../data/DataThin_Corrected_MixCorrected/%s/MixCorrected_RightBg_%s_0_1Mult.root",partName[iPart].Data(),partName[iPart].Data()));
    // fHist2D->Add((TH2D *) fFileConj->Get(Form("fHistCorrected_%s_pt2_ptTrigg1",partName[iPart+1].Data())));
    // fHist2D->Divide((TH2D *) fFileConj->Get(Form("fHistCorrected_%s_pt2_ptTrigg1",partName[iPart].Data())));
    // TFile * fFileConj1 = new TFile(Form("../data/DataThin_Corrected_MixCorrected/%s/MixCorrected_LeftBg_%s_0_1Mult.root",partName[iPart].Data(),partName[iPart].Data()));
    // fHist2D->Add((TH2D *) fFileConj1->Get(Form("fHistCorrected_%s_pt2_ptTrigg1",partName[iPart].Data())),-1);
  // }

  TGaxis::SetMaxDigits(4);
  TCanvas * can = Plotter::CreateCanvas("can");
  gPad->SetTheta(45);
  gPad->SetPhi(40);
  gPad->GetFrame()->SetLineColor(0);
  // gPad->SetLeftMargin(0.18);
  Plotter::Set2DHistAxes(fHist2D,"#Delta#eta","#Delta#varphi","","");
  // Plotter::Set2DHistAxes(fHist2D,"#Delta#eta","#Delta#varphi","#frac{d#font[52]{N}_{pair}}{d#Delta#varphi d#Delta#eta}","");

  fHist2D->GetXaxis()->SetRangeUser(-1.3,1.3);
  // fHist2D->GetZaxis()->SetRangeUser(400,2000);
  // fHist2D->GetZaxis()->SetRangeUser(25000000,35000000);
  fHist2D->DrawCopy("surf2 fb");
  fHist2D->DrawCopy("surf same fb");

  TPaveText * paveGen = new TPaveText();
  Plotter::SetPaveText(paveGen,42,0.05, 0, 0,12,0,0.02,0.45, 0.8,0.97);
  paveGen->AddText(Form("ALICE, Work in progress"));
  paveGen->AddText("pp, 13.6 TeV");
  if(iPart==4||iPart==2||iPart==6)paveGen->AddText(Form("h-(%s+%s)",partSym[iPart].Data(),partSym[iPart+1].Data()));
  else paveGen->AddText(Form("h-%s",partSym[iPart].Data()));
  paveGen->Draw();

  TPaveText * pavePart = new TPaveText();
  Plotter::SetPaveText(pavePart,42,0.05, 0, 0,32,0,0.7,0.9999, 0.8,0.99);
  pavePart->AddText("4 < #font[52]{p}_{T}^{trigg} < 6 GeV/#font[52]{c}");
  pavePart->AddText("2 < #font[52]{p}_{T}^{assoc} < 3 GeV/#font[52]{c}");
  pavePart->AddText("FT0C: 0-1%");
  pavePart->Draw();
}
