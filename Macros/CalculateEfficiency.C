#include "Plotter.h"

void CalculateEfficiency(Int_t part = 4){

  TString name[]={"K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus","Pion"};
  TString finalNames[]={"K_{S}^{0}","#Lambda","#bar{#Lambda}","#Xi^{-}","#Xi^{+}","#Omega^{-}","#Omega^{+})","#pi^{+}+#pi^{-}"};

  TFile * fFile = new TFile(Form("../data/AnalysisResults_MC_20_09_Cascades.root"));

  TH2F * genHist  = (TH2F *) fFile->Get(Form("correlate-strangeness/Generated/h%s",name[part].Data()));
  genHist->Sumw2();
  TH2F * genHist_wPV_yCut  = (TH2F *) fFile->Get(Form("correlate-strangeness/GeneratedWithPV/h%s_MidYVsMult",name[part].Data()));
  genHist_wPV_yCut->Sumw2();
  TH2F * genHist_wPV  = (TH2F *) fFile->Get(Form("correlate-strangeness/GeneratedWithPV/h%s",name[part].Data()));
  genHist_wPV->Sumw2();
  TH3F * recoHist  = (TH3F *) fFile->Get(Form("correlate-strangeness/h%sEtaVsPtVsPhi",name[part].Data()));
  recoHist->Sumw2();
  TH3F * recoHist_Bckg  = (TH3F *) fFile->Get(Form("correlate-strangeness/h%sEtaVsPtVsPhiBg",name[part].Data()));
  recoHist_Bckg->Sumw2();

  TH3F * recoHist_SpectrY  = (TH3F *) fFile->Get(Form("correlate-strangeness/h3d%sSpectrumY",name[part].Data()));
  recoHist_SpectrY->Sumw2();

  TH2F * recoEtaPtProj = (TH2F*) recoHist->Project3D("yx");
  recoEtaPtProj->SetName(Form("%s",name[part].Data()));
  TH2F * recoEtaPtProj_Bckg = (TH2F*) recoHist_Bckg->Project3D("yx");

  TH2F * recoEtaPtProj_woBckg = (TH2F *) recoEtaPtProj->Clone();
  recoEtaPtProj_woBckg->SetName(Form("%s_woBckg",name[part].Data()));
  recoEtaPtProj_woBckg->Add(recoEtaPtProj_Bckg,-1);

  TH2F * eff_woBckg_wPV = (TH2F *) recoEtaPtProj_woBckg->Clone();
  eff_woBckg_wPV->SetName(Form("%s_woBckg_wPV",name[part].Data()));
  eff_woBckg_wPV->Divide(genHist_wPV);
  TH2F * eff_wPV = (TH2F *) recoEtaPtProj->Clone();
  eff_wPV->SetName(Form("%s_wPV",name[part].Data()));
  eff_wPV->Divide(genHist_wPV);

  recoEtaPtProj->RebinY(2);
  genHist->RebinY(2);
  genHist_wPV->RebinY(2);
  recoEtaPtProj_woBckg->RebinY(2);
  recoEtaPtProj->Divide(genHist);
  recoEtaPtProj_woBckg->Divide(genHist);
  TCanvas * can = Plotter::CreateCanvas("c");
  Plotter::Set2DHistAxes(recoEtaPtProj,"p_{T}","#eta","eff","");
  Plotter::Set2DHistAxes(recoEtaPtProj_woBckg,"p_{T}","#eta","eff","");

  recoEtaPtProj->DrawCopy("surf2");

  TCanvas * can2 = Plotter::CreateCanvas("c2");
  recoEtaPtProj_woBckg->DrawCopy("surf2");

  TCanvas * can3 = Plotter::CreateCanvas("c3");
  Plotter::Set2DHistAxes(eff_wPV,"p_{T}","#eta","eff","");
  eff_wPV->DrawCopy("surf2");

  TCanvas * can4 = Plotter::CreateCanvas("c4");
  Plotter::Set2DHistAxes(eff_woBckg_wPV,"p_{T}","#eta","eff","");
  eff_woBckg_wPV->DrawCopy("surf2");

  //1d Projections
  recoHist->GetYaxis()->SetRangeUser(-0.8,0.8);
  recoHist_Bckg->GetYaxis()->SetRangeUser(-0.8,0.8);
  genHist->GetYaxis()->SetRangeUser(-0.8,0.8);

  TH1F * pt_eff = (TH1F*) recoHist->Project3D("x");
  pt_eff->SetName(Form("pt_%s",name[part].Data()));

  TH1F * pt_eff_wPv = (TH1F*) recoHist->Project3D("x");
  pt_eff_wPv->SetName(Form("pt_%s_wPv",name[part].Data()));

  TH1F * pt_eff_woBckg = (TH1F*) recoHist->Project3D("x");
  pt_eff_woBckg->Add((TH1F*) recoHist_Bckg->Project3D("x"),-1);
  pt_eff_woBckg->SetName(Form("pt_%s_woBckg",name[part].Data()));

  TH1F * pt_eff_woBckg_wPV = (TH1F*) recoHist->Project3D("x");
  pt_eff_woBckg_wPV->Add((TH1F*) recoHist_Bckg->Project3D("x"),-1);
  pt_eff_woBckg_wPV->SetName(Form("pt_%s_woBckg_wPV",name[part].Data()));


  TH1F * genPtProj = (TH1F*) genHist->ProjectionX();
  genPtProj->SetName("genPtProj");
  genHist_wPV->GetYaxis()->SetRangeUser(-0.8,0.8);
  TH1F * genPtProj_wPv = (TH1F*) genHist_wPV->ProjectionX();
  genPtProj_wPv->SetName("genPtProj_wPv");

  TH1F * genPtProj_wPv_yCut = (TH1F*) genHist_wPV_yCut->ProjectionX();
  genPtProj_wPv_yCut->SetName("genPtProj_wPv_yCut");

  pt_eff->Divide(genPtProj);
  pt_eff_wPv->Divide(genPtProj_wPv);
  pt_eff_woBckg->Divide(genPtProj);
  pt_eff_woBckg_wPV->Divide(genPtProj_wPv);

  TCanvas * canPtEff = Plotter::CreateCanvas("canPtEff");

  Plotter::SetHist(pt_eff_wPv,"",24,kBlack,1.2);
  Plotter::SetHistAxes(pt_eff,"#font[12]{p}_{T} (GeV/#font[12]{c})","#varepsilon");
  Plotter::SetHist(pt_eff,"",20,kBlack,1.2);
  Plotter::SetHist(pt_eff_woBckg,"",20,kBlue+1,1.2);
  Plotter::SetHist(pt_eff_woBckg_wPV,"",24,kBlue+1,1.2);
  pt_eff->DrawCopy();
  pt_eff_wPv->DrawCopy("same");
  pt_eff_woBckg->DrawCopy("same");
  pt_eff_woBckg_wPV->DrawCopy("same");

  TPaveText * pave = new TPaveText();
  Plotter::SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.97, 0.7,0.95);
  pave->AddText("ALICE, Work in Progress");
  pave->AddText("pp, 13.6 TeV");
  pave->AddText(Form("%s, |#eta| < 0.8",finalNames[part].Data()));
  pave->Draw("same");

  TLegend *leg= Plotter::CreateLegend(0.55, 0.92, 0.45, 0.65,0.05);
  leg->AddEntry(pt_eff,"Gen: all","p");
  leg->AddEntry(pt_eff_wPv,"Gen: only with PV","p");
  leg->AddEntry(pt_eff_woBckg,"Gen: all, Rec: bckg. subtracted","p");
  leg->AddEntry(pt_eff_woBckg_wPV,"Gen: only with PV, Rec: bckg. subtracted","p");
  leg->Draw();

  TH1F * eta_eff = (TH1F*) recoHist->Project3D("y");
  eta_eff->SetName(Form("eta_%s",name[part].Data()));

  TH1F * eta_eff_wPv = (TH1F*) recoHist->Project3D("y");
  eta_eff_wPv->SetName(Form("eta_%s_wPV",name[part].Data()));

  TH1F * eta_eff_woBckg = (TH1F*) recoHist->Project3D("y");
  eta_eff_woBckg->Add((TH1F*) recoHist_Bckg->Project3D("y"),-1);
  eta_eff_woBckg->SetName(Form("eta_%s_woBckg",name[part].Data()));

  TH1F * eta_eff_woBckg_wPV = (TH1F*) recoHist->Project3D("y");
  eta_eff_woBckg_wPV->Add((TH1F*) recoHist_Bckg->Project3D("y"),-1);
  eta_eff_woBckg_wPV->SetName(Form("eta_%s_woBckg_wPV",name[part].Data()));

  TH1F * genEtaProj = (TH1F*) genHist->ProjectionY();
  genEtaProj->SetName("genEtaProj");
  TH1F * genEtaProj_wPv = (TH1F*) genHist_wPV->ProjectionY();
  genEtaProj_wPv->SetName("genPtProj_wPv");
  eta_eff_wPv->RebinX(2);
  eta_eff->RebinX(2);
  eta_eff_woBckg->RebinX(2);
  eta_eff_woBckg_wPV->RebinX(2);

  eta_eff->Divide(genEtaProj);
  eta_eff_wPv->Divide(genEtaProj_wPv);
  eta_eff_woBckg->Divide(genEtaProj);
  eta_eff_woBckg_wPV->Divide(genEtaProj_wPv);

  TCanvas * canEtaEff = Plotter::CreateCanvas("canEtaEff");
  Plotter::SetHistAxes(eta_eff,"#eta","#varepsilon");
  Plotter::SetHist(eta_eff,"",20,kBlack,1.2);
  Plotter::SetHist(eta_eff_wPv,"",24,kBlack,1.2);
  Plotter::SetHist(eta_eff_woBckg,"",20,kBlue+1,1.2);
  Plotter::SetHist(eta_eff_woBckg_wPV,"",24,kBlue+1,1.2);
  eta_eff->GetYaxis()->SetRangeUser(0,0.2999);
  eta_eff->DrawCopy("");
  eta_eff_wPv->DrawCopy("same");
  eta_eff_woBckg->DrawCopy("same");
  eta_eff_woBckg_wPV->DrawCopy("same");
  pave->Draw("same");
  leg->Draw();

  recoHist_SpectrY->GetZaxis()->SetRange(2,2);
  TH1F * pt_eff_ycut_woBckg_wPV = (TH1F *) recoHist_SpectrY->Project3D("x");
  pt_eff_ycut_woBckg_wPV->SetName(Form("pt_%s_ycut_woBckg_wPV",name[part].Data()));
  recoHist_SpectrY->GetZaxis()->SetRange(1,1);
  TH1F * pt_Bckg = (TH1F *) recoHist_SpectrY->Project3D("x");
  pt_Bckg->SetName("pt_Bckg");
  recoHist_SpectrY->GetZaxis()->SetRange(3,3);
  pt_Bckg->Add((TH1F *) recoHist_SpectrY->Project3D("x"));

  pt_eff_ycut_woBckg_wPV->Add(pt_Bckg,-1);
  pt_eff_ycut_woBckg_wPV->Divide(genPtProj_wPv_yCut);

  TFile *fileNew = new TFile(Form("../data/Efficiency/Eff_%s.root",name[part].Data()), "RECREATE");
  eta_eff->Write();
  eta_eff_wPv->Write();
  eta_eff_woBckg->Write();
  eta_eff_woBckg_wPV->Write();
  pt_eff->Write();
  pt_eff_wPv->Write();
  pt_eff_woBckg->Write();
  pt_eff_woBckg_wPV->Write();
  recoEtaPtProj_woBckg->Write();
  eff_wPV->Write();
  eff_woBckg_wPV->Write();
  recoEtaPtProj->Write();
  pt_eff_ycut_woBckg_wPV->Write();
  fileNew->Close();

}
