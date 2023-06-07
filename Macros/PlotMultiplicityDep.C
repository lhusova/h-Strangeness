#include "Plotter.h"

void PlotMultiplicityDep(Int_t region = 2, Int_t part = 0){

  TString multiplicityNames[]={"minBias","0_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_60Mult","60_70Mult","70_80Mult","80_90Mult","90_100Mult"};
  TString multiplicityLeg[]={"MB","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"};
  TString nameSave[]={"K0s","Lam","Xi","Omega"};
  TString histName[]={"fHistNear","fHistAway","fHistUE"};
  TString finalNames[]={"K_{S}^{0}","(#Lambda+#bar{#Lambda})","(#Xi^{+}+#Xi^{-})","(#Omega^{+}+#Omega^{-})"};
  TString regionName[]={"Near-side, |#Delta#varphi|<0.9","Away-side, |#Delta#varphi-#pi|<1.4","Underlying event"};
  TString regionSaveName[]={"Near","Away","UE"};
  TFile *fFile[11];
  TH1F * yieldHist[11];

  TCanvas *can = Plotter::CreateCanvas("c");
  gStyle->SetErrorX(0.01);
  can->GetPadSave()->SetLogy();

  gStyle->SetPalette(55);
  Double_t nColorptassMult = gStyle->GetNumberOfColors();
  Int_t iColptassMult =0.;
  Int_t colors[10];
  for (Int_t i=0; i<10; i++) {
      iColptassMult = nColorptassMult/10*i;
      colors[i]=gStyle->GetColorPalette(iColptassMult);
  }

  TLegend *leg = Plotter::CreateLegend(0.15, 0.45, 0.15, 0.6,0.05);

  for (size_t i = 0; i < 11; i++) {
    fFile[i] = new TFile(Form("../data/Yields_%s_%s.root",nameSave[part].Data(),multiplicityNames[i].Data()));
    yieldHist[i] = (TH1F *)fFile[i]->Get(Form("%s",histName[region].Data()));
    if(region==2)yieldHist[i]->Scale(50.);
    if(i==0){
      Plotter::SetHist(yieldHist[i],"",20,kBlack,1.);
      yieldHist[i]->DrawCopy();
    }
    else {
      Plotter::SetHist(yieldHist[i],"",20,colors[10-i],1.);
      yieldHist[i]->DrawCopy("same");
    }
    leg->AddEntry(yieldHist[i],Form("%s",multiplicityLeg[i].Data()),"pl");
  }
  yieldHist[0]->DrawCopy("same");
  leg->Draw();

  TPaveText *pave = new TPaveText();
  Plotter::SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.95, 0.65,0.95);
  pave->AddText("ALICE, Work in Progress");
  pave->AddText("pp, 13.6 TeV");
  pave->AddText(Form("h-%s",finalNames[part].Data()));
  pave->AddText(Form("3 < #font[12]{p}^{trigg}_{T} < 20 GeV/#font[12]{c}"));
  pave->AddText("|#Delta#eta| < 1");
  pave->AddText(Form("%s",regionName[region].Data()));
  pave->Draw("same");

  can->SaveAs(Form("../Plots/Yield_%s_mult_%s.pdf",nameSave[part].Data(),regionSaveName[region].Data()));

  TCanvas *canRatio = Plotter::CreateCanvas("c1");
  gStyle->SetErrorX(0.01);

  TLegend *legRatio = Plotter::CreateLegend(0.15, 0.45, 0.53, 0.93,0.05);

  for (size_t i = 1; i < 11; i++) {
    yieldHist[i]->Divide(yieldHist[0]);
    if(i==1){
      yieldHist[i]->GetYaxis()->SetRangeUser(-0.59,3.99);
      yieldHist[i]->GetYaxis()->SetTitle("Y/Y_{MB}");
      yieldHist[i]->DrawCopy();
    }
    else {
      yieldHist[i]->DrawCopy("same");
    }
    legRatio->AddEntry(yieldHist[i],Form("%s",multiplicityLeg[i].Data()),"pl");
  }
  legRatio->Draw();
  pave->Draw("same");
  Plotter::DrawUnity(kBlack);

  canRatio->SaveAs(Form("../Plots/Yield_%s_multRatio_%s.pdf",nameSave[part].Data(),regionSaveName[region].Data()));
}
