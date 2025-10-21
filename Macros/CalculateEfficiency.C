#include "Plotter.h"
#include "Definitions.h"

void CalculateEfficiency(Int_t part = 4){

  Int_t rebinX = 2;
  Int_t rebinY = 8;
  TString name[]={"Trigger","K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus","Pion"};
  TString finalNames[]={"charged hadrons","K_{S}^{0}","#Lambda","#bar{#Lambda}","#Xi^{-}","#Xi^{+}","#Omega^{-}","#Omega^{+}","#pi^{+}+#pi^{-}"};

  TFile * fFile = new TFile(Form("../data/AnalysisResults_Eff_LHC25b4b6.root"));

  TH3F * genHist  = (TH3F *) fFile->Get(Form("h-strange-correlation/GeneratedWithPV/h%s",name[part].Data()));//_id31232
  genHist->RebinY(rebinY);
  genHist->RebinX(rebinX);
  genHist->Sumw2();
  TH3F * recoHist;
  if(part==0)  recoHist = (TH3F *) fFile->Get(Form("h-strange-correlation/h%sPrimaryEtaVsPt",name[part].Data()));//_requireTrigger_3
  else recoHist = (TH3F *) fFile->Get(Form("h-strange-correlation/h%sEtaVsPtVsPhi",name[part].Data()));
  recoHist->Sumw2();
  recoHist->RebinY(rebinY);
  recoHist->RebinX(rebinX);
  TH3F * recoHist_Bckg;
  if(part>0)  {
    recoHist_Bckg= (TH3F *) fFile->Get(Form("h-strange-correlation/h%sEtaVsPtVsPhiBg",name[part].Data()));//_requireTrigger_3
    recoHist_Bckg->Sumw2();
    recoHist_Bckg->RebinY(rebinY);
    recoHist_Bckg->RebinX(rebinX);
  }
  TFile *fileNew = new TFile(Form("../data/Efficiency/Eff_LHC25b4b6_%s.root",name[part].Data()), "RECREATE");
  if(part==0) {
    TH3F * hist3Deff = (TH3F *) recoHist->Clone();
    hist3Deff->SetName(Form("hEfficiency%sMultTmp",name[part].Data()));
    hist3Deff->Divide(genHist);
    const Int_t nBinsX = hist3Deff->GetXaxis()->GetNbins();
    double binsx[nBinsX+1];
    for (int i =0; i<hist3Deff->GetXaxis()->GetNbins(); i++) {
      binsx[i]=hist3Deff->GetXaxis()->GetBinLowEdge(i+1);
      if(i==hist3Deff->GetXaxis()->GetNbins()-1)
        binsx[i+1]=hist3Deff->GetXaxis()->GetBinUpEdge(i+1);
    }
    const Int_t nBinsY = hist3Deff->GetYaxis()->GetNbins()+2;
    double binsy[nBinsY+1];
    binsy[0] = hist3Deff->GetYaxis()->GetBinLowEdge(1)-0.2;
    for (int i =1; i<hist3Deff->GetYaxis()->GetNbins()+1; i++) {
      binsy[i]=hist3Deff->GetYaxis()->GetBinLowEdge(i);
      // cout << binsy[i] << endl;
      if(i==hist3Deff->GetYaxis()->GetNbins()){
        binsy[i+1]=hist3Deff->GetYaxis()->GetBinUpEdge(i);
        binsy[i+2]=hist3Deff->GetYaxis()->GetBinUpEdge(i)+0.2;
        // cout << binsy[i+1] << endl;
        // cout << binsy[i+2] << endl;
      }
      
    }
    const Int_t nBinsZ = hist3Deff->GetZaxis()->GetNbins()+2;
    double binsz[nBinsZ+1];
    binsz[0] = hist3Deff->GetZaxis()->GetBinLowEdge(1)-10;
    for (int i =1; i<hist3Deff->GetZaxis()->GetNbins()+1; i++) {
      binsz[i]=hist3Deff->GetZaxis()->GetBinLowEdge(i);
      if(i==hist3Deff->GetZaxis()->GetNbins()){
        binsz[i+1]=hist3Deff->GetZaxis()->GetBinUpEdge(i);
        binsz[i+2]=hist3Deff->GetZaxis()->GetBinUpEdge(i)+10;
      }
    }
    TH3F * histEffExtended = new TH3F(Form("hEfficiency%sMult",name[part].Data()),"",nBinsX,binsx,nBinsY,binsy,nBinsZ,binsz);
    for (int i =1; i<histEffExtended->GetXaxis()->GetNbins()+1; i++) {
      for (int j =1; j<histEffExtended->GetYaxis()->GetNbins()+1; j++) {
        for (int k=1; k<histEffExtended->GetZaxis()->GetNbins()+1; k++) {
          if(j==1&&k==1){
            histEffExtended->SetBinContent(i,j,k,hist3Deff->GetBinContent(i,j,k));
            histEffExtended->SetBinError(i,j,k,hist3Deff->GetBinError(i,j,k));
          }
          else if(j==1&&k==histEffExtended->GetZaxis()->GetNbins()){
            histEffExtended->SetBinContent(i,j,k,hist3Deff->GetBinContent(i,j,k-2));
            histEffExtended->SetBinError(i,j,k,hist3Deff->GetBinError(i,j,k-2));
          }
          else if(k==1&&histEffExtended->GetYaxis()->GetNbins()){
            histEffExtended->SetBinContent(i,j,k,hist3Deff->GetBinContent(i,j-2,k));
            histEffExtended->SetBinError(i,j,k,hist3Deff->GetBinError(i,j-2,k));
          }
          else if(j==1){
            histEffExtended->SetBinContent(i,j,k,hist3Deff->GetBinContent(i,j,k-1));
            histEffExtended->SetBinError(i,j,k,hist3Deff->GetBinError(i,j,k-1));
          }
          else if(k==1){
            histEffExtended->SetBinContent(i,j,k,hist3Deff->GetBinContent(i,j-1,k));
            histEffExtended->SetBinError(i,j,k,hist3Deff->GetBinError(i,j-1,k));
          }
          else if(j==histEffExtended->GetYaxis()->GetNbins()&&k==histEffExtended->GetZaxis()->GetNbins()){
            histEffExtended->SetBinContent(i,j,k,hist3Deff->GetBinContent(i,j-2,k-2));
            histEffExtended->SetBinError(i,j,k,hist3Deff->GetBinError(i,j-2,k-2));
          }
          else if (j==histEffExtended->GetYaxis()->GetNbins()){
            histEffExtended->SetBinContent(i,j,k,hist3Deff->GetBinContent(i,j-2,k-1));
            histEffExtended->SetBinError(i,j,k,hist3Deff->GetBinError(i,j-2,k-1));
          }
          else if (k==histEffExtended->GetZaxis()->GetNbins()){
            histEffExtended->SetBinContent(i,j,k,hist3Deff->GetBinContent(i,j-1,k-2));
            histEffExtended->SetBinError(i,j,k,hist3Deff->GetBinError(i,j-1,k-2));
          }
          else{
            histEffExtended->SetBinContent(i,j,k,hist3Deff->GetBinContent(i,j-1,k-1));
            histEffExtended->SetBinError(i,j,k,hist3Deff->GetBinError(i,j-1,k-1));
          }
        }
      }
    }

    TH3F * hist3DeffUncert = (TH3F *) histEffExtended->Clone();
    hist3DeffUncert->SetName(Form("hEfficiencyUncertainty%sMult",name[part].Data()));
    for(int i =1; i<hist3DeffUncert->GetXaxis()->GetNbins()+1; i++){
      for(int j =1; j< hist3DeffUncert->GetYaxis()->GetNbins()+1; j++){
        for(int k=1; k<hist3DeffUncert->GetZaxis()->GetNbins()+1; k++){
          hist3DeffUncert->SetBinContent(i,j,k, histEffExtended->GetBinError(i,j,k));
        }
      }
    }
    histEffExtended->Write();
    hist3DeffUncert->Write();
  }
  TH2F * genEtaPtProj = (TH2F*)genHist ->Project3D("yx");
  genEtaPtProj->SetName("genEtaPtProj");
  TH2F * recoEtaPtProj = (TH2F*) recoHist->Project3D("yx");
  recoEtaPtProj->SetName(Form("hEfficiency%s",name[part].Data()));
  if(part>0) {
    TH2F * recoEtaPtProj_Bckg = (TH2F*) recoHist_Bckg->Project3D("yx");
    recoEtaPtProj_Bckg->SetName(Form("hEfficiency_bckg%s",name[part].Data()));
    recoEtaPtProj->Add(recoEtaPtProj_Bckg,-1);
  }

  recoEtaPtProj->Divide(genEtaPtProj);
  TCanvas * can = CreateCanvas("c");
  gPad->SetTheta(45);
  gPad->SetPhi(40);
  gPad->GetFrame()->SetLineColor(0);
  Set2DHistAxes(recoEtaPtProj,"#font[52]{p}_{T} (GeV/#font[52]{c})","#eta","#varepsilon","");
  // if(part==0) recoEtaPtProj->GetXaxis()->SetRangeUser(2,50);
  recoEtaPtProj->GetZaxis()->SetTitleOffset(1);
  recoEtaPtProj->GetXaxis()->SetTitleOffset(1.5);
  recoEtaPtProj->GetYaxis()->SetTitleOffset(1.5);
  // recHistClone->DrawCopy("lego2z");
  recoEtaPtProj->GetXaxis()->SetRangeUser(0.21,50);
  recoEtaPtProj->DrawCopy("surf2 fb");
  recoEtaPtProj->DrawCopy("surf same fb");

  TPaveText * pave2D = new TPaveText();
  SetPaveText(pave2D,42,0.05, 0, 0,12,0,0.02,0.45, 0.8,0.97);
  pave2D->AddText(Form("ALICE, Work in progress"));
  pave2D->AddText("pp, 5.36 TeV");
  pave2D->AddText("LHC25b4b6");
  pave2D->AddText(Form("%s",finalNames[part].Data()));
  pave2D->Draw();

  TH2F * effUncert = (TH2F*)recoEtaPtProj->Clone();
  effUncert->SetName(Form("hEfficiencyUncertainty%s",name[part].Data()));
  for (int x=1; x<effUncert->GetXaxis()->GetNbins()+1; x++) {
    for (int y=1; y<effUncert->GetYaxis()->GetNbins()+1; y++) {
      effUncert->SetBinContent(x,y,recoEtaPtProj->GetBinError(x,y));
    }
  }
  TCanvas * canUncert = CreateCanvas("cUncert");
  gPad->SetTheta(45);
  gPad->SetPhi(40);
  gPad->GetFrame()->SetLineColor(0);
  effUncert->GetZaxis()->SetTitleOffset(1);
  effUncert->GetXaxis()->SetTitleOffset(1.5);
  effUncert->GetYaxis()->SetTitleOffset(1.5);
  // recHistClone->DrawCopy("lego2z");
  effUncert->GetZaxis()->SetTitle("#varepsilon_{uncertanity}");
  effUncert->GetXaxis()->SetRangeUser(0.21,50);
  effUncert->DrawCopy("surf2 fb");
  effUncert->DrawCopy("surf same fb");

  //1d Projections
  recoHist->GetYaxis()->SetRangeUser(-0.8,0.8);
  // if(part==0)recoHist->GetXaxis()->SetRangeUser(2,50);
  genHist->GetYaxis()->SetRangeUser(-0.8,0.8);
  // if(part==0)genHist->GetXaxis()->SetRangeUser(2,50);
  
  TH1F * pt_eff = (TH1F*) recoHist->Project3D("x");
  pt_eff->SetName(Form("pt_%s",name[part].Data()));
  TH1F * genPtProj = (TH1F*) genHist->Project3D("x");
  genPtProj->SetName("genPtProj");
  pt_eff->Divide(genPtProj);
  
  TCanvas * canPtEff;
  if(part==0) canPtEff = CreateCanvas("canPtEff",600,750,true);
  else canPtEff = CreateCanvas("canPtEff");
  TPad *padRatioPtEff;
  if(part==0){ 
    padRatioPtEff = new TPad("padRatioPtEff","padRatioPtEff",0.001,0.001,0.999,0.3);
    padRatioPtEff->SetMargin(0.12,0.02,0.25,0.01);
    padRatioPtEff->Draw();
    canPtEff->GetPadSave()->cd();
  }
  SetHistAxes(pt_eff,"#font[52]{p}_{T} (GeV/#font[52]{c})","#varepsilon");
  SetHist(pt_eff,"",20,kBlack,1.2);
  pt_eff->DrawCopy();
  TPaveText * pave = new TPaveText();
  SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.97, 0.7,0.95);
  pave->AddText("ALICE, Work in Progress");
  pave->AddText("pp, 5.36 TeV");
  pave->AddText("LHC25b4b6");
  pave->AddText(Form("%s",finalNames[part].Data()));
  pave->Draw("same");

  TLegend *legEff = new TLegend(0.6,0.6,0.9,0.9);
  legEff->SetBorderSize(0);
  legEff->SetFillStyle(0);

  if(part==0) {
    TH1F * effMult[nMultBins-1];
    TH1F * effMultRatio[nMultBins-1];
    TH1F * gen1Dmult[nMultBins-1];

    TH2F * eff2DMult[nMultBins-1];
    TH2F * gen2Dmult[nMultBins-1];
    TCanvas * canEffMult[nMultBins-1];
    for(Int_t iMult=0; iMult< nMultBins-1; iMult++){
      // if(iMult==1){
        // recoHist->GetZaxis()->SetRange(iMult, iMult+1);
        // genHist->GetZaxis()->SetRange(iMult, iMult+1);
      // }else{
        recoHist->GetZaxis()->SetRange(iMult+1, iMult+1);
        genHist->GetZaxis()->SetRange(iMult+1, iMult+1);
      // }
      effMult[iMult] = (TH1F *)recoHist->Project3D("x");
      effMult[iMult]->SetName(Form("eff%s%d",name[part].Data(), iMult));
      gen1Dmult[iMult] = (TH1F *)genHist->Project3D("x");
      gen1Dmult[iMult]->SetName(Form("gen1Dmult%s%d",name[part].Data(), iMult));
      effMult[iMult]->Write();
      effMult[iMult]->Divide(gen1Dmult[iMult]);
      SetHist(effMult[iMult], "", 24, colorsMultiplicity[iMult], 1.2);
      canPtEff->GetPadSave()->cd();
      effMult[iMult]->DrawCopy("same");
      legEff->AddEntry(effMult[iMult],Form("%s",multiplicityPave[iMult+1].Data()),"pl");
      effMult[iMult]->Write();
      if(iMult==nMultBins-2) {
        legEff->Draw("same");
      }
      padRatioPtEff->cd();
      effMultRatio[iMult] = (TH1F*) effMult[iMult]->Clone();
      effMultRatio[iMult]->SetName(Form("effMultRatio%i",iMult));
      effMultRatio[iMult] ->Divide(pt_eff);
      if(iMult==0){
        SetHistAxesSmallPad(effMultRatio[iMult], "#font[52]{p}_{T} (GeV/#font[52]{c})", "ratio over MB");
        effMultRatio[iMult]->DrawCopy();
      }else effMultRatio[iMult]->DrawCopy("same");

      canEffMult[iMult]=CreateCanvas(Form("eff2dMult%i",iMult));
      eff2DMult[iMult] = (TH2F *)recoHist->Project3D("yx");
      eff2DMult[iMult]->SetName(Form("eff2DMult%s%d",name[part].Data(), iMult));
      gen2Dmult[iMult] = (TH2F *)genHist->Project3D("yx");
      gen2Dmult[iMult]->SetName(Form("gen2Dmult%s%d",name[part].Data(), iMult));
      eff2DMult[iMult]->Divide(gen2Dmult[iMult]);
      eff2DMult[iMult]->SetTitle(Form("%i mult bin",iMult));
      eff2DMult[iMult]->DrawCopy("surf2");
    }
  }

  TH1F * eta_eff = (TH1F*) recoHist->Project3D("y");
  eta_eff->SetName(Form("eta_%s",name[part].Data()));
  //
  TH1F * genEtaProj = (TH1F*) genHist->Project3D("y");
  genEtaProj->SetName("genEtaProj");
  eta_eff->Divide(genEtaProj);
  //
  TCanvas * canEtaEff = CreateCanvas("canEtaEff");
  SetHistAxes(eta_eff,"#eta","#varepsilon");
  SetHist(eta_eff,"",20,kBlack,1.2);
  eta_eff->DrawCopy("");
  // if(part==0) pave->AddText(Form("%s, 2 < #font[52]{p}_{T} < 50 GeV/#font[52]{c}",finalNames[part].Data()));
  pave->Draw("same");

  
  eta_eff->Write();
  pt_eff->Write();
  recoEtaPtProj->Write();
  effUncert->Write();
  if(part==0){
    
  }
  // fileNew->Close();

}
