#include <TPaveText.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <THn.h>
#include <TLegend.h>
#include "Plotter.h"

void ScaleMixing(TH2F * histMix){
    Int_t nPhiBins = histMix->GetYaxis()->GetNbins();
    Int_t bin0 =histMix->GetXaxis()->FindBin(0.);
    Double_t scale=0.;

    for(Int_t iBinPhi=nPhiBins/2; iBinPhi<nPhiBins; iBinPhi++){
        scale+=histMix->GetBinContent(bin0,iBinPhi+1);
    }
    scale=scale/(nPhiBins/2);

    if(scale!=0) histMix->Scale(0.5/scale);
    //cout << "scale     "<< scale << endl;

};

int IsAtEdge(TH1F* hist, float value){
  int returnValue = -1;
  for(Int_t ii=1; ii<hist->GetNbinsX()+1; ii++){
    if(TMath::Abs(hist->GetBinLowEdge(ii)-value)<1e-6) returnValue = ii; // returns bin if at edge
  }
  return returnValue;
};

void Convolution(Int_t mult_i = 3, Int_t iPtTrigger = 3){


  TFile * fFile = new TFile (Form("../data/AnalysisResults_ForMCclosure_14_08.root"));
  TDirectoryFile *my_dir = (TDirectoryFile*)fFile->Get("correlate-strangeness_id14337/ClosureTest");
  TFile * fFile2 = new TFile (Form("../data/AnalysisResultsMCClosure_GenExpanded.root"));

  TString particle[8]={"Pion","K0Short","Lambda","AntiLambda","XiMinus","XiPlus","OmegaMinus","OmegaPlus"};
  TString histTrigger = "fHist3dTrigger";
  TString particle1[8]={"Pion","K_{S}^{0}","#Lambda","#bar{#Lambda}","Xi^{-}","#Xi^{+}","#Omega^{-}","#Omega^{+}"};
  TString histAsso[8] = {"hPion","hK0Short","hLambda","hAntiLambda","hXiMinus","hXiPlus","hOmegaMinus","hOmegaPlus"};
  TString histmix[8] = {"fHist5d2pcMixedPions","fHist5d2pcMixedK0Short","fHist5d2pcMixedLambda","fHist5d2pcMixedAntiLambda","fHist5d2pcMixedXiMinus","fHist5d2pcMixedXiPlus","fHist5d2pcMixedOmegaMinus","fHist5d2pcMixedOmegaPlus"};
  TString histsame[8] = {"fHist5d2pcPions","fHist5d2pcK0Short","fHist5d2pcLambda","fHist5d2pcAntiLambda","fHist5d2pcXiMinus","fHist5d2pcXiPlus","fHist5d2pcOmegaMinus","fHist5d2pcOmegaPlus"};
  Int_t ptTriggerEdge[] = {2,4,6,10,100};
  TString ptTriggerBinRange[4] = {"2-4","4-6","6-10","10-100"};
  TString ptBinRange[9] = {"0-1","1-2","2-3","3-4","4-6","6-8","8-10","10-12","12-15"};
  TString region[]={"Signal","RightBg","LeftBg"};
  const Int_t particles = 1;
  const Int_t nPtBins = 9;
  TH1 * fHistTrigger;
  TH1 * fHistAsso[particles];
  TH1F * fHistConvolution[particles][nPtBins];
  TH2F * fHistConvolutionMix[particles][nPtBins];
  TH2F * fHistConvolutionSame[particles][nPtBins];
  TH2F * fHistCorrected[particles][nPtBins];
  Double_t ptBins[nPtBins+1]={0.,1.,2.,3.,4.,6.,8.,10.,12.,15.};

  TFile * fFileNew = TFile::Open (Form("../data/Convolution_PtTrigger%d.root",iPtTrigger),"RECREATE");

  TH3D* trigger = (TH3D*) my_dir->Get("hTrigger");
  //trigger->Write();

  fHistTrigger = (TH1F*) trigger->ProjectionY("triggerEta",trigger->GetXaxis()->FindBin(ptTriggerEdge[iPtTrigger]), trigger->GetXaxis()->FindBin(ptTriggerEdge[iPtTrigger+1]), 1, trigger->GetZaxis()->GetNbins());
  for(int iPtAsso=0; iPtAsso<nPtBins; iPtAsso++){
    //fHistTrigger -> Write();

    for(int iPart=1; iPart<particles+1; iPart++){

      // TCanvas *c1 = new TCanvas("c1", "", 800,600);

      THnF* mix = (THnF*) fFile2->Get(Form("Same_%s__Expanded",particle[iPart].Data()));
      fHistConvolutionMix[iPart-1][iPtAsso] = (TH2F*) mix->Projection(0,1);
      fHistConvolutionMix[iPart-1][iPtAsso]->SetName(Form("Mixed_%s_Expanded_%i",particle[iPart].Data(),iPtAsso));
      //to advoid bining edge warning
      for(int i=1; i<fHistConvolutionMix[iPart-1][iPtAsso]->GetXaxis()->GetNbins()+1; i++){
        for(int j=1; j<fHistConvolutionMix[iPart-1][iPtAsso]->GetYaxis()->GetNbins()+1; j++){
          fHistConvolutionMix[iPart-1][iPtAsso]->SetBinContent(i,j,0);
        }
      }
      int Nbins2d = fHistConvolutionMix[iPart-1][iPtAsso]->GetYaxis()->GetNbins();
      // fHistConvolutionMix[iPart-1][iPtAsso]->SetName(Form("fHistConvolution2d%s_pt%s_PtTriggerbin%s", particle[iPart].Data(),ptBinRange[iPtAsso].Data(),ptTriggerBinRange[iPtTrigger].Data()));
      TH3F* asso = (TH3F*) my_dir->Get(histAsso[iPart].Data());
      int assoBinPtLeft = asso->GetXaxis()->FindBin(ptBins[iPtAsso]);
      int assoBinPtRight = asso->GetXaxis()->FindBin(ptBins[iPtAsso+1]);
      fHistAsso[iPart-1] = asso->ProjectionY(histAsso[iPart].Data(),assoBinPtLeft,assoBinPtRight,1, asso->GetZaxis()->GetNbins());
      // //fHistAsso[iPart-1] -> Write();
      //
      fHistConvolution[iPart-1][iPtAsso] = new TH1F(Form("fHistConvolution%i",iPtAsso),"fHistConvolution" , 64,-1.60,1.60);
      //
      for(Int_t ib1=1; ib1<fHistTrigger->GetNbinsX()+1; ib1++){
        for(Int_t ib2=1; ib2<fHistAsso[iPart-1]->GetNbinsX()+1; ib2++){
          float thisDeltaEta = fHistTrigger->GetBinCenter(ib1)-fHistAsso[iPart-1]->GetBinCenter(ib2);
          float multiplication = fHistTrigger->GetBinContent(ib1)*fHistAsso[iPart-1]->GetBinContent(ib2);

          //Determine if bin matches an edge precisely
          int isThisAtEdge = IsAtEdge(fHistConvolution[iPart-1][iPtAsso], thisDeltaEta);

          if(isThisAtEdge>-1){
            //divide into adjacent bins if necessary
            fHistConvolution[iPart-1][iPtAsso]->Fill(thisDeltaEta+1e-3, 0.5*multiplication);
            fHistConvolution[iPart-1][iPtAsso]->Fill(thisDeltaEta-1e-3, 0.5*multiplication);
            int etabin1 = fHistConvolutionMix[iPart-1][iPtAsso]->GetXaxis()->FindBin(thisDeltaEta+1e-3);
            int etabin2 = fHistConvolutionMix[iPart-1][iPtAsso]->GetXaxis()->FindBin(thisDeltaEta-1e-3);
            for(Int_t ib3 = 1; ib3<Nbins2d+1; ib3++){
              double originalValue1 = fHistConvolutionMix[iPart-1][iPtAsso]->GetBinContent(etabin1, ib3);
              fHistConvolutionMix[iPart-1][iPtAsso]->SetBinContent(etabin1, ib3, originalValue1+0.5*multiplication);
            }
            for(Int_t ib3 = 1; ib3<Nbins2d+1; ib3++){
              double originalValue2 = fHistConvolutionMix[iPart-1][iPtAsso]->GetBinContent(etabin2, ib3);
              fHistConvolutionMix[iPart-1][iPtAsso]->SetBinContent(etabin2, ib3, originalValue2+0.5*multiplication);
            }
          }
          else{
            fHistConvolution[iPart-1][iPtAsso]->Fill(thisDeltaEta, multiplication);
            int etabin = fHistConvolutionMix[iPart-1][iPtAsso]->GetXaxis()->FindBin(thisDeltaEta);
            for(Int_t ib3 = 1; ib3<Nbins2d+1; ib3++){
              double originalValue = fHistConvolutionMix[iPart-1][iPtAsso]->GetBinContent(etabin, ib3);
              fHistConvolutionMix[iPart-1][iPtAsso]->SetBinContent(etabin, ib3, originalValue+multiplication);
            }
          }
        }
      }

      //fHistConvolutionMix[iPart-1][iPtAsso]->RebinY(2);
      //fHistConvolutionMix[iPart-1][iPtAsso]->RebinX(2);
      ScaleMixing(fHistConvolutionMix[iPart-1][iPtAsso]);
      // //fHistConvolutionMix[iPart-1][iPtAsso]->SetTitle(Form("fHistConvolution2d%s_pt%s_PtTriggerbin%s", particle[iPart].Data(),ptBinRange[iPtAsso].Data(),ptTriggerBinRange[iPtTrigger].Data()));
      fHistConvolutionMix[iPart-1][iPtAsso]->Write(Form("fHistConvolution2d%s_pt%i_PtTriggerbin%i", particle[iPart].Data(),iPtAsso,iPtTrigger));
      //
      // THnF* same = (THnF*) fFile2->Get(Form("Same_%s_Gen",particle[iPart].Data()));
      //
      // same->GetAxis(2)->SetRange(iPtAsso+1,iPtAsso+1);
      // same->GetAxis(3)->SetRange(iPtTrigger+1,iPtTrigger+1);
      // fHistConvolutionSame[iPart-1][iPtAsso] = (TH2F *) same->Projection(0,1);
      // fHistConvolutionSame[iPart-1][iPtAsso]->Scale(1./fHistConvolutionSame[iPart-1][iPtAsso]->GetXaxis()->GetBinWidth(2));
      // fHistConvolutionSame[iPart-1][iPtAsso]->Scale(1./fHistConvolutionSame[iPart-1][iPtAsso]->GetYaxis()->GetBinWidth(2));
      // //fHistConvolutionSame[iPart-1][iPtAsso]->RebinY(2);
      // //fHistConvolutionSame[iPart-1][iPtAsso]->RebinX(2);
      // fHistConvolutionSame[iPart-1][iPtAsso]->Write(Form("SamePro_%s_Pt%d_PtTrigger%s",particle[iPart].Data(),iPtAsso,ptTriggerBinRange[iPtTrigger].Data()));
      // fHistConvolutionSame[iPart-1][iPtAsso]->Scale(fHistConvolutionSame[iPart-1][iPtAsso]->GetXaxis()->GetBinWidth(2));
      // fHistConvolutionSame[iPart-1][iPtAsso]->Scale(fHistConvolutionSame[iPart-1][iPtAsso]->GetYaxis()->GetBinWidth(2));
      // fHistConvolutionSame[iPart-1][iPtAsso]->Divide(fHistConvolutionMix[iPart-1][iPtAsso]);
      // //fHistConvolutionSame[iPart-1][iPtAsso]->SetTitle(Form("fHistCorrected_%s_pt%s_ptTrigger%s",particle[iPart].Data(),ptBinRange[iPtAsso].Data(),ptTriggerBinRange[iPtTrigger].Data()));
      // fHistConvolutionSame[iPart-1][iPtAsso]->SetName(Form("fHistCorrected_%s_pt%s_ptTrigger%s",particle[iPart].Data(),ptBinRange[iPtAsso].Data(),ptTriggerBinRange[iPtTrigger].Data()));
      // fHistConvolutionSame[iPart-1][iPtAsso]->Scale(1./fHistConvolutionSame[iPart-1][iPtAsso]->GetXaxis()->GetBinWidth(2));
      // fHistConvolutionSame[iPart-1][iPtAsso]->Scale(1./fHistConvolutionSame[iPart-1][iPtAsso]->GetYaxis()->GetBinWidth(2));
      // fHistConvolutionSame[iPart-1][iPtAsso]->Write();

      //TPaveText *pave = new TPaveText();
      //Plotter::SetPaveText(pave,42,0.04, 0, 0,15,0,0.15,0.35,0.65,0.80);
      //pave->SetTextSize(0.04);
      //pave->AddText(Form("%s",particle1[iPart-1].Data()));
      //pave->AddText(Form("Multiplicity: %s",multBinRange[mult_i-1].Data()));
      //pave->AddText(Form("PtAsso: %d-%d GeV/c",iPtAsso-1,iPtAsso));
      //pave->AddText(Form("PtTrigger: %s GeV/c",ptTriggerBinRange[iPtTrigger-1].Data()));
      //pave->Draw("same");
      //c1->SetTitle(Form("%s_PtAsso%d_PtTrigger%s_Mult%s",particle[iPart-1].Data(),iPtAsso,ptTriggerBinRange[iPtTrigger-1].Data(),multBinRange[mult_i-1].Data()));
      //c1->Write(Form("%s_PtAsso%d_PtTrigger%s_Mult%s",particle[iPart-1].Data(),iPtAsso,ptTriggerBinRange[iPtTrigger-1].Data(),multBinRange[mult_i-1].Data()));
      //c1->SaveAs(Form("../data/%s/Plots/Convolution/%s_PtAsso%d_PtTrigger%s_Mult%s.pdf",mode[modenum].Data(),particle[iPart-1].Data(),iPtAsso,ptTriggerBinRange[iPtTrigger-1].Data(),multBinRange[mult_i-1].Data()));

      // delete c1;
      // delete fHistConvolution[iPart-1][iPtAsso];
      // delete fHistAsso[iPart-1];
      // delete fHistConvolutionSame[iPart-1][iPtAsso];
      // delete fHistConvolutionMix[iPart-1][iPtAsso];

    }
  }
  fFileNew->Close();

}
