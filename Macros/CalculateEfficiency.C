#include "Plotter.h"

void CalculateEfficiency(Int_t part =4){

  TString name[]={"K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus","Pion"};
  TString finalNames[]={"K_{S}^{0}","#Lambda","#bar{#Lambda}","#Xi^{-}","#Xi^{+}","#Omega^{-}","#Omega^{+})","#pi^{+}+#pi^{-}"};

  TFile * fFile = new TFile(Form("../data/AnalysisResults_MC_11_09_Cascades.root"));

  TH2F * genHist  = (TH2F *) fFile->Get(Form("correlate-strangeness/Generated/h%s",name[part].Data()));
  genHist->Sumw2();
  TH2F * genHist_wPV  = (TH2F *) fFile->Get(Form("correlate-strangeness/GeneratedWithPV/h%s",name[part].Data()));
  genHist_wPV->Sumw2();
  TH3F * recoHist  = (TH3F *) fFile->Get(Form("correlate-strangeness/h%sEtaVsPtVsPhi",name[part].Data()));
  recoHist->Sumw2();

  TH2F * recoEtaPtProj = (TH2F*) recoHist->Project3D("yx");

  recoEtaPtProj->RebinY(2);
  genHist->RebinY(2);
  genHist_wPV->RebinY(2);
  recoEtaPtProj->Divide(genHist);
  TCanvas * can = Plotter::CreateCanvas("c");
  Plotter::Set2DHistAxes(recoEtaPtProj,"p_{T}","#eta","eff","");

  recoEtaPtProj->DrawCopy("surf2");

  //1d Projections
  recoHist->GetYaxis()->SetRangeUser(-0.8,0.8);
  TH1F * eff = (TH1F*) recoHist->Project3D("x");
  eff->SetName("eff");
  TH1F * eff_wPv = (TH1F*) recoHist->Project3D("x");
  eff_wPv->SetName("eff_wPv");
  genHist->GetYaxis()->SetRangeUser(-0.8,0.8);
  TH1F * genPtProj = (TH1F*) genHist->ProjectionX();
  genPtProj->SetName("genPtProj");
  genHist_wPV->GetYaxis()->SetRangeUser(-0.8,0.8);
  TH1F * genPtProj_wPv = (TH1F*) genHist_wPV->ProjectionX();
  genPtProj_wPv->SetName("genPtProj_wPv");
  eff->Divide(genPtProj);
  eff_wPv->Divide(genPtProj_wPv);

  TCanvas * canPtEff = Plotter::CreateCanvas("canPtEff");

  Plotter::SetHist(eff_wPv,"",24,kBlack,1.2);
  // eff->DrawCopy("");
  // eff_wPv->DrawCopy("same");
  // TH1F* eff = (TH1F*) recoEtaPtProj->ProjectionX();
  // eff->Scale(1./recoEtaPtProj->GetYaxis()->GetNbins());
  Plotter::SetHistAxes(eff,"#font[12]{p}_{T} (GeV/#font[12]{c})","#varepsilon");
  Plotter::SetHist(eff,"",20,kBlack,1.2);
  eff->DrawCopy();
  eff_wPv->DrawCopy("same");

  TPaveText * pave = new TPaveText();
  Plotter::SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.97, 0.7,0.95);
  pave->AddText("ALICE, Work in Progress");
  pave->AddText("pp, 13.6 TeV");
  pave->AddText(Form("%s",finalNames[part].Data()));
  pave->Draw("same");

  TLegend *leg= Plotter::CreateLegend(0.55, 0.92, 0.55, 0.65,0.05);
  leg->AddEntry(eff,"all generated","p");
  leg->AddEntry(eff_wPv,"generated only with PV","p");
  leg->Draw();

  TH1F * eff_eta = (TH1F*) recoHist->Project3D("y");
  eff_eta->SetName("eff_eta");
  TH1F * eff_eta_wPv = (TH1F*) recoHist->Project3D("y");
  eff_eta_wPv->SetName("eff_eta_wPv");
  TH1F * genEtaProj = (TH1F*) genHist->ProjectionY();
  genEtaProj->SetName("genEtaProj");
  TH1F * genEtaProj_wPv = (TH1F*) genHist_wPV->ProjectionY();
  genEtaProj_wPv->SetName("genPtProj_wPv");
  eff_eta_wPv->RebinX(2);
  eff_eta->RebinX(2);
  eff_eta->Divide(genEtaProj);
  eff_eta_wPv->Divide(genEtaProj_wPv);

  TCanvas * canEtaEff = Plotter::CreateCanvas("canEtaEff");
  Plotter::SetHistAxes(eff_eta,"#eta","#varepsilon");
  Plotter::SetHist(eff_eta,"",20,kBlack,1.2);
  Plotter::SetHist(eff_eta_wPv,"",24,kBlack,1.2);
  eff_eta->GetYaxis()->SetRangeUser(0,0.2999);
  eff_eta->DrawCopy("");
  eff_eta_wPv->DrawCopy("same");
  pave->Draw("same");
  leg->Draw();

}
