#include "Plotter.h"
void ScaleMixing(TH2F * histMix);

void MixingCorrection(TString particle="K0Short",TString region = "RightBg",Int_t mult_i = 2){
  TFile * fFile = new TFile ("../data/LongTrain_2Aug/AnalysisResults_expanded.root");
  // TString multBinsName[]={"0_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_60Mult","60_70Mult","70_80Mult","80_90Mult","90_100Mult"};
  TString multBinsName[]={"0_1Mult","1_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_70Mult","70_100Mult"};
  // TString multBinsName[]={"MB"};

  THnF* same = (THnF*) fFile->Get(Form("Same_%s_%s_Expanded",particle.Data(),region.Data()));
  // same->GetAxis(5)->SetRange(mult_i+1,mult_i+1);

  // TFile * fFileMix = new TFile ("../data/Convolution_MBruns.root");
  THnF* mix = (THnF*) fFile->Get(Form("Mixed_%s_%s_Expanded",particle.Data(),region.Data()));
  // if(!mix) cout << "no mix!!" << endl;
  // mix->GetAxis(5)->SetRange(mult_i+1,mult_i+1);

  const Int_t nPtBins = 9;
  const Int_t nPtTriggBins = 4;
  TH2F * fHistCorrected[nPtBins][nPtTriggBins];
  TH2F * fHistSame[nPtBins][nPtTriggBins];
  TH2F * fHistMix[nPtBins][nPtTriggBins];
  TH2F * fHistMixProjection[nPtBins][nPtTriggBins][10];
  TH2F * fHistSameProjection[nPtBins][nPtTriggBins][10];

  TFile * fFileNew = TFile::Open (Form("../data/LongTrain_2Aug/%s/MixCorrected_%s_%s_%s.root",particle.Data(),region.Data(),particle.Data(),multBinsName[mult_i].Data()),"RECREATE");
  for (Int_t iPtTrigg = 0; iPtTrigg < nPtTriggBins; iPtTrigg++) {
    same->GetAxis(3)->SetRange(iPtTrigg+1,iPtTrigg+1);
    mix->GetAxis(3)->SetRange(iPtTrigg+1,iPtTrigg+1);
    for (Int_t iPt = 0; iPt < nPtBins; iPt++) {

      same->GetAxis(4)->SetRange(1,10);
      mix->GetAxis(4)->SetRange(1,10);

      same->GetAxis(2)->SetRange(iPt+1,iPt+1);
      mix->GetAxis(2)->SetRange(iPt+1,iPt+1);

      fHistSame[iPt][iPtTrigg] = (TH2F *) same->Projection(0,1);
      fHistSame[iPt][iPtTrigg]->SetName(Form("fHistSame_%s_pt%d_ptTrigg%d",particle.Data(),iPt,iPtTrigg));
      // fHistSame[iPt][iPtTrigg]->RebinY(2);
      // fHistSame[iPt][iPtTrigg]->RebinX(2);
      fHistSame[iPt][iPtTrigg]->Scale(1./fHistSame[iPt][iPtTrigg]->GetXaxis()->GetBinWidth(2));
      fHistSame[iPt][iPtTrigg]->Scale(1./fHistSame[iPt][iPtTrigg]->GetYaxis()->GetBinWidth(2));
      fHistSame[iPt][iPtTrigg]->Write();

      fHistMix[iPt][iPtTrigg] = (TH2F *) mix->Projection(0,1);
      // fHistMix[iPt][iPtTrigg] = (TH2F *)fFileMix->Get(Form("fHistConvolution2d%s_pt%i_PtTriggerbin%i",particle.Data(),iPt,iPtTrigg));
      fHistMix[iPt][iPtTrigg]->SetName(Form("fHistMixed_%s_pt%d_ptTrigg%d",particle.Data(),iPt,iPtTrigg));
      // fHistMix[iPt][iPtTrigg]->RebinY(2);
      // fHistMix[iPt][iPtTrigg]->RebinX(2);
      ScaleMixing(fHistMix[iPt][iPtTrigg]);
      fHistMix[iPt][iPtTrigg]->Write();

      for (Int_t iPvz = 0; iPvz < 10; iPvz++) {
        same->GetAxis(4)->SetRange(iPvz+1,iPvz+1);
        mix->GetAxis(4)->SetRange(iPvz+1,iPvz+1);
        // Int_t iPvz=0;
        fHistSameProjection[iPt][iPtTrigg][iPvz] = (TH2F *) same->Projection(0,1);
        fHistSameProjection[iPt][iPtTrigg][iPvz]->SetName(Form("fHistSame_%s_pt%d%d_ptTrigg%d",particle.Data(),iPvz,iPt,iPtTrigg));
        // fHistSameProjection[iPt][iPtTrigg][iPvz]->RebinY(2);
        // fHistSameProjection[iPt][iPtTrigg][iPvz]->RebinX(2);

        fHistMixProjection[iPt][iPtTrigg][iPvz] = (TH2F *) mix->Projection(0,1);
        // fHistMixProjection[iPt][iPtTrigg][iPvz] = (TH2F *)fFileMix->Get(Form("fHistConvolution2d%s_pt%i_PtTriggerbin%i",particle.Data(),iPt,iPtTrigg));

        fHistMixProjection[iPt][iPtTrigg][iPvz]->SetName(Form("fHistMix_%s_pt%d%d_ptTrigg%d",particle.Data(),iPvz,iPt,iPtTrigg));
        fHistMixProjection[iPt][iPtTrigg][iPvz]->Write();
        // fHistMixProjection[iPt][iPtTrigg][iPvz]->RebinY(2);
        // fHistMixProjection[iPt][iPtTrigg][iPvz]->RebinX(2);

        ScaleMixing(fHistMixProjection[iPt][iPtTrigg][iPvz]);
        fHistSameProjection[iPt][iPtTrigg][iPvz]->Divide(fHistMixProjection[iPt][iPtTrigg][iPvz]);
        fHistSameProjection[iPt][iPtTrigg][iPvz]->Write();

        if(iPvz==0){
          fHistCorrected[iPt][iPtTrigg]=(TH2F*)fHistSameProjection[iPt][iPtTrigg][iPvz]->Clone();
          fHistCorrected[iPt][iPtTrigg]->SetName(Form("fHistCorrected_%s_pt%d_ptTrigg%d",particle.Data(),iPt,iPtTrigg));
        }else fHistCorrected[iPt][iPtTrigg]->Add(fHistSameProjection[iPt][iPtTrigg][iPvz]);

      }

      // fHistCorrected[iPt][iPtTrigg]->Scale(1./fHistCorrected[iPt][iPtTrigg]->GetXaxis()->GetBinWidth(2));
      // fHistCorrected[iPt][iPtTrigg]->Scale(1./fHistCorrected[iPt][iPtTrigg]->GetYaxis()->GetBinWidth(2));

      Plotter::Set2DHistAxes(fHistCorrected[iPt][iPtTrigg],"#Delta#eta","#Delta#varphi","#frac{d^{2}N}{d#Delta#eta d#Delta#varphi}",Form("h-%s, p_{T}^{assoc} bin: %d, p_{T}^{trigg} bin: %d",particle.Data(),iPt,iPtTrigg));
      fHistCorrected[iPt][iPtTrigg]->Write();
    }

  }
  fFileNew->Close();
}
//____________________________________________________________________
void ScaleMixing(TH2F * histMix){
    Int_t nPhiBins = histMix->GetYaxis()->GetNbins();
    Int_t bin0 =histMix->GetXaxis()->FindBin(0.);
    Double_t scale=0.;

    for(Int_t iBinPhi=nPhiBins/2; iBinPhi<nPhiBins; iBinPhi++){
        scale+=histMix->GetBinContent(bin0,iBinPhi+1);
    }
    scale=scale/(nPhiBins/2);

    if(scale!=0) histMix->Scale(1./scale);
    cout << "scale     "<< scale << endl;

};
