#include "Plotter.h"

void CompareYieldSpectra(Int_t part=0,Int_t multClass =2){

  TString finalNames[]={"K_{S}^{0}","(#Lambda+#bar{#Lambda})","(#Xi^{+}+#Xi^{-})","(#Omega^{+}+#Omega^{-})","#pi^{+}+#pi^{-}"};
  TString nameSave[]={"K0s","Lam","Xi","Omega","Pion"};
  // TString multiplicityNames[]={"minBias","0_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_60Mult","60_70Mult","70_80Mult","80_90Mult","90_100Mult"};
  // TString multiplicityPave[]={"MB","0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"};
  TString multiplicityNames[]={"minBias","0_1Mult","1_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_70Mult","70_100Mult"};
  TString multiplicityPave[]={"MB","0-1%","1-10%","10-20%","20-30%","30-40%","40-50%","50-70%","70-100%"};

  const Int_t n = 3;
  TFile *fFile[n];
  TString file_Suffix[n]={"_fullrangePeak_11_flat","_fullrangePeak_11_long_range","_fullrangePeak_11_flow_modulation_v2"};//,"_fullrangePeak_11_flow_modulation_v3"};//{"_fullrange","_fullrangePeak_12","_fullrangePeak_11","_fullrangePeak_1","_fullrangePeak_09"};
  TString legend[n] = {"ZYAM","long range","flow: v2"};//, "flow: v2+v3"};//{"|#Delta#eta| < 1.4","|#Delta#eta| < 1.2","|#Delta#eta| < 1.1","|#Delta#eta| < 1","|#Delta#eta| < 0.9"};
  TH1F * spectrum[n];
  TH1F * ratio[n-1];
  Color_t col[n]={kGreen+2,kBlue,kMagenta};//,kViolet-5};//,kPink-5};
  Int_t marker[n]={20,20,20};//,30};//,28};

  TCanvas *can = Plotter::CreateCanvas(Form("c"),700,750,true);

  TPad * pad2 = new TPad("pad2","pad2",0.001,0.001,0.999,0.3);
  pad2->SetMargin(0.12,0.02,0.25,0.01);
  pad2->Draw();
  pad2->SetTicky();
  pad2->SetTickx();

  TLegend *leg = Plotter::CreateLegend(0.7, 0.95, 0.15, 0.45,0.05);

  for (Int_t i = 0; i < n; i++) {
    fFile[i] = new TFile(Form("../data/DataThin_Uncorrected_MixCorrected/Yields/%s/Yields_%s_%s%s_ptTrigg1.root",nameSave[part].Data(),nameSave[part].Data(),multiplicityNames[multClass].Data(),file_Suffix[i].Data()));
    spectrum[i] = (TH1F *) fFile[i]->Get("fHistUE");

    can->GetPadSave()->cd();
    spectrum[i]->GetYaxis()->SetMaxDigits(2);
    Plotter::SetHist(spectrum[i],"",marker[i],col[i],1.);
    if(i==0) {
      spectrum[i]->DrawCopy();
      TPaveText *pave = new TPaveText();
      Plotter::SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.97, 0.55,0.95);
      pave->AddText("ALICE, Work in Progress");
      pave->AddText("pp, 13.6 TeV");
      pave->AddText(Form("%s",multiplicityPave[multClass].Data()));
      pave->AddText(Form("h-%s",finalNames[part].Data()));
      pave->AddText(Form("1 < #font[52]{p}^{trigg}_{T} < 2 GeV/#font[52]{c}"));
      pave->AddText("|#Delta#eta| < 1.1");
      // pave->AddText("Near-side, |#Delta#varphi|<#pi/2");
      pave->AddText("Underlying event");
      pave->Draw("same");
    }
    else spectrum[i]->DrawCopy("same");
    leg->AddEntry(spectrum[i],legend[i].Data(),"pl");
    if(i==n-1)leg->Draw("same");

    if(i>0){
      pad2->cd();
      ratio[i-1] = (TH1F *) spectrum[i]->Clone();
      ratio[i-1] ->SetName(Form("ratio%i",i));
      ratio[i-1]->Divide(spectrum[0]);
      if(i==1){
        Plotter::SetHistAxesSmallPad(ratio[i-1],spectrum[0]->GetXaxis()->GetTitle(),"ratio to full range");
        ratio[i-1]->GetYaxis()->SetRangeUser(0.950001,1.04999);
        ratio[i-1]->DrawCopy();
      }
      else ratio[i-1]->DrawCopy("same");
      Plotter::DrawUnity(kBlack,1,0,15);
    }
  }
// can->SaveAs(Form("../Plots/Checks/background/Yields_Near_%s_%s.pdf",nameSave[part].Data(),multiplicityPave[multClass].Data()));

}
