#include "Plotter.h"

void Plot2DCorrelationFunction(TString particle="K0Short",TString region = "RightBg",Int_t ptBin=2, TString option = "corrected"){

  TFile * fFile = new TFile (Form("../data/MixCorrected_%s_%s.root",region.Data(),particle.Data()));

  TH2F * hist2D;
  TCanvas * can = Plotter::CreateCanvas("c");
  gPad->SetMargin(0.15,0.03,0.12,0.03);
  if(option=="corrected"){
    hist2D = (TH2F *) fFile->Get(Form("fHistCorrected_%s_pt%d",particle.Data(),ptBin));
    hist2D->GetXaxis()->SetRangeUser(-1,1);
    hist2D->SetTitle("");
  }
  else{
    if(option=="mix")hist2D = (TH2F *) fFile->Get(Form("fHistMixed_%s_pt%d",particle.Data(),ptBin));
    if(option=="same")  {
      hist2D = (TH2F *) fFile->Get(Form("fHistSame_%s_pt%d",particle.Data(),ptBin));
      hist2D->Scale(1./hist2D->GetXaxis()->GetBinWidth(2));
      hist2D->Scale(1./hist2D->GetYaxis()->GetBinWidth(2));
    }
    Plotter::Set2DHistAxes(hist2D,"#Delta#eta","#Delta#varphi","#frac{d^{2}N}{d#Delta#eta d#Delta#varphi}","");
  }

  hist2D->GetZaxis()->SetMaxDigits(2);
  hist2D->DrawCopy("surf2");
  hist2D->DrawCopy("surf same");

  can->SaveAs(Form("../Plots/2Dfunction_%s_%s_%s.pdf",particle.Data(),region.Data(),option.Data()));
}
