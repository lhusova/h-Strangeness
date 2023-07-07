#include "Plotter.h"

void CompareBackgroundDPhi(Int_t part=4,Int_t multClass =2,Int_t pt=7){

  TString finalNames[]={"K_{S}^{0}","(#Lambda+#bar{#Lambda})","(#Xi^{+}+#Xi^{-})","(#Omega^{+}+#Omega^{-})","(#pi^{+}+#pi^{-})"};
  TString nameSave[]={"K0s","Lam","Xi","Omega","Pion"};
  TString multiplicityNames[]={"minBias","0_1Mult","1_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_70Mult","70_100Mult"};
  TString multiplicityPave[]={"MB","0-1%","1-10%","10-20%","20-30%","30-40%","40-50%","50-70%","70-100%"};
  Double_t ptBins[]={0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,5.,6.,7.,8.,10.,12,15};

  const Int_t n = 4;
  TFile *fFile[n];
  TString file_Suffix[n]={"_fullrangePeak_11_flat","_fullrangePeak_11_long_range","_fullrangePeak_11_flow_modulation_v2","_fullrangePeak_11_flow_modulation_v3"};
  TString legend[n] = {"ZYAM","long range","flow: v2", "flow: v2+v3"};
  Color_t colors[n] = {kBlack, kBlue+1,kMagenta,kGreen+2};

  TH1F * phiProj_bck[n];

  for (Int_t i = 0; i < n; i++) {
    fFile[i] = new TFile(Form("../data/Yields/%s/Yields_%s_%s%s.root",nameSave[part].Data(),nameSave[part].Data(),multiplicityNames[multClass].Data(),file_Suffix[i].Data()));
    phiProj_bck[i] = (TH1F *) fFile[i]->Get(Form("underlying_ev_pT%d",pt));
    phiProj_bck[i]->SetName(Form("underlying_ev_pT%d%d",pt,i));
  }

  TH1F * phiProj = (TH1F *) fFile[0]->Get(Form("phiProj_withBckg_pT%d",pt));

  TCanvas *can = Plotter::CreateCanvas("c");
  gStyle->SetErrorX(0.01);
  phiProj->SetName("phiProj");
  Plotter::SetHistAxes(phiProj,"#Delta#varphi (rad)","1/#font[12]{N}_{Trigg}d#font[12]{N}/d#Delta#varphi (rad^{-1})");
  Plotter::SetHist(phiProj,"",20,kBlack,0.8);
  phiProj->DrawCopy();

  TLegend *leg = Plotter::CreateLegend(0.45, 0.65, 0.35, 0.65,0.05);

  for (Int_t i = 0; i < n; i++) {
    Plotter::SetHist(phiProj_bck[i],"",20,colors[i],0.8);
    phiProj_bck[i]->DrawCopy("same");
    leg->AddEntry(phiProj_bck[i],legend[i].Data(),"pl");
  }
  leg->Draw();

  TPaveText *pave = new TPaveText();
  Plotter::SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.97, 0.55,0.95);
  pave->AddText("ALICE, Work in Progress");
  pave->AddText("pp, 13.6 TeV");
  pave->AddText(Form("%s",multiplicityPave[multClass].Data()));
  pave->AddText(Form("h-%s",finalNames[part].Data()));
  pave->AddText(Form("3 < #font[12]{p}^{trigg}_{T} < 20 GeV/#font[12]{c}"));
  pave->AddText(Form("%g < #font[12]{p}^{assoc}_{T} < %g GeV/#font[12]{c}",ptBins[pt],ptBins[pt+1]));
  pave->AddText("|#Delta#eta| < 1.1");
  pave->Draw("same");

}
