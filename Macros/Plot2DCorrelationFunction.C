#include "Plotter.h"

void Plot2DCorrelationFunction(Int_t pTTrigBin = 1,Int_t ptBin=1,TString particle="K0Short",TString region = "signal", TString option = "corrected"){

  TFile * fFile = new TFile (Form("../data/MixCorrected/%s/MixCorrected_%s_%s_1_10Mult.root",particle.Data(),region.Data(),particle.Data()));
  // Float_t ptRanges[]={0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,5.,6.,7.,8.,10.,12,15};
  Double_t ptRanges[]={0.,1.,2.,3.,4.,6.,8.,10.,12,15};
  Double_t ptRangesTrigg[]={0.,1.,2.,3.,100};

  TH2F * hist2D;
  TCanvas * can = Plotter::CreateCanvas("c");
  gPad->SetMargin(0.15,0.03,0.12,0.03);
  if(option=="corrected"){
    hist2D = (TH2F *) fFile->Get(Form("fHistCorrected_%s_pt%d_ptTrigg%d",particle.Data(),ptBin,pTTrigBin));
    // hist2D->GetXaxis()->SetRangeUser(-1.45,1.45);
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
  hist2D->SetTitle(Form("%g < p_{T}^{trigg} < %g GeV/#font[12]{c}, %g < p_{T}^{assoc} < %g GeV/#font[12]{c}",ptRangesTrigg[pTTrigBin],ptRangesTrigg[pTTrigBin+1],ptRanges[ptBin],ptRanges[ptBin+1]));
  hist2D->DrawCopy("surf2");
  hist2D->DrawCopy("surf same");

  can->SaveAs(Form("../Plots/2Dfunction_%s_%s_%s_%i_%i_1_10Mult.pdf",particle.Data(),region.Data(),option.Data(),ptBin,pTTrigBin));
}
