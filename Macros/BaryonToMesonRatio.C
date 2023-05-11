#include "Plotter.h"

void BaryonToMesonRatio(){

  TFile * fFileK0s = new TFile("../data/Yields_K0s.root");
  TFile * fFileLam = new TFile("../data/Yields_Lam.root");

  TH1D *histK0s[3];
  TH1D *histLam[3];
  TH1D *histRatio[3];

  TString names[3] = {"fHistNear","fHistAway","fHistUE"};
  Color_t col[3] = {kRed+1,kBlue+1,kGreen+1};

  TCanvas *can = Plotter::CreateCanvas("c");
  gStyle->SetErrorX(0.01);

  for (Int_t iReg = 0; iReg < 3; iReg++) {
    histK0s[iReg]=(TH1D *) fFileK0s->Get(names[iReg].Data());
    histK0s[iReg]->SetName(Form("%s_K0s",names[iReg].Data()));
    histLam[iReg]=(TH1D *) fFileLam->Get(names[iReg].Data());
    histLam[iReg]->SetName(Form("%s_Lam",names[iReg].Data()));

    histK0s[iReg]->Scale(2.);
    histRatio[iReg] = (TH1D *) histLam[iReg]->Clone();
    histRatio[iReg]->SetName(Form("%s_ratio",names[iReg].Data()));
    histRatio[iReg]->Divide(histK0s[iReg]);

    Plotter::SetHistAxes(histRatio[iReg],"#font[12]{p}^{assoc}_{T} (GeV/#font[12]{c})","Y^{h-(#Lambda+#bar{#Lambda})}_{#Delta#varphi}/2Y^{h-K_{S}^{0}}_{#Delta#varphi}");
    Plotter::SetHist(histRatio[iReg],"",20,col[iReg],1.);

    if(iReg==0){
      histRatio[iReg]->GetYaxis()->SetRangeUser(0,0.4);
      histRatio[iReg]->DrawCopy();
    }
    else histRatio[iReg]->DrawCopy("same");
  }

  TPaveText *pave = new TPaveText();
  Plotter::SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.95, 0.65,0.95);
  pave->AddText("ALICE, Work in Progress");
  pave->AddText("pp, 13.6 TeV");
  pave->AddText(Form("3 < #font[12]{p}^{trigg}_{T} < 20 GeV/#font[12]{c}"));
  pave->AddText("|#Delta#eta| < 1");
  pave->Draw("same");

  TLegend *leg = Plotter::CreateLegend(0.15, 0.45, 0.15, 0.45,0.05);
  leg->AddEntry(histRatio[0],"Near-side, |#Delta#varphi|<0.9","pl");
  leg->AddEntry(histRatio[1],"Away-side, |#Delta#varphi-#pi|<1.4","pl");
  leg->AddEntry(histRatio[2],"Underlying event","pl");
  leg->Draw("same");

  can->SaveAs("../Plots/LambdaOverK0sRatio.pdf");
}
