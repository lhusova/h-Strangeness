#include "Plotter.h"
void ScaleMixing(TH2F * histMix);

void MixingCorrection(TString particle="K0Short",TString region = "Signal"){
  TFile * fFile = new TFile ("../data/AnalysisResults_Hyperloop_09_05.root");

  THnF* same = (THnF*) fFile->Get(Form("correlate-strangeness/sameEvent/%s/%s",region.Data(),particle.Data()));
  THnF* mix = (THnF*) fFile->Get(Form("correlate-strangeness/mixedEvent/%s/%s",region.Data(),particle.Data()));

  const Int_t nPtBins = 7;
  TH2F * fHistCorrected[nPtBins];
  TH2F * fHistSame[nPtBins];
  TH2F * fHistMix[nPtBins];
  TH2F * fHistMixProjection[nPtBins][10];
  TH2F * fHistSameProjection[nPtBins][10];

  TFile * fFileNew = TFile::Open (Form("../data/MixCorrected_%s_%s.root",region.Data(),particle.Data()),"RECREATE");
  for (Int_t iPt = 0; iPt < nPtBins; iPt++) {

    same->GetAxis(2)->SetRange(iPt+1,iPt+1);
    mix->GetAxis(2)->SetRange(iPt+1,iPt+1);

    fHistSame[iPt] = (TH2F *) same->Projection(0,1);
    fHistSame[iPt]->SetName(Form("fHistSame_%s_pt%d",particle.Data(),iPt));
    fHistSame[iPt]->Write();

    fHistMix[iPt] = (TH2F *) mix->Projection(0,1);
    fHistMix[iPt]->SetName(Form("fHistMixed_%s_pt%d",particle.Data(),iPt));
    ScaleMixing(fHistMix[iPt]);
    fHistMix[iPt]->Write();

    for (Int_t iPvz = 0; iPvz < 10; iPvz++) {
      same->GetAxis(3)->SetRange(iPvz+1,iPvz+1);
      mix->GetAxis(3)->SetRange(iPvz+1,iPvz+1);

      fHistSameProjection[iPt][iPvz] = (TH2F *) same->Projection(0,1);
      fHistSameProjection[iPt][iPvz]->SetName(Form("fHistCorrected_%s_pt%d%d",particle.Data(),iPvz,iPt));

      fHistMixProjection[iPt][iPvz] = (TH2F *) mix->Projection(0,1);
      fHistMixProjection[iPt][iPvz]->SetName(Form("fHistMix_%s_pt%d%d",particle.Data(),iPvz,iPt));

      ScaleMixing(fHistMixProjection[iPt][iPvz]);
      fHistSameProjection[iPt][iPvz]->Divide(fHistMixProjection[iPt][iPvz]);

      if(iPvz==0){
        fHistCorrected[iPt]=(TH2F*)fHistSameProjection[iPt][iPvz]->Clone();
        fHistCorrected[iPt]->SetName(Form("fHistCorrected_%s_pt%d",particle.Data(),iPt));
      }else fHistCorrected[iPt]->Add(fHistSameProjection[iPt][iPvz]);

    }

    fHistCorrected[iPt]->Scale(1./fHistCorrected[iPt]->GetXaxis()->GetBinWidth(2));
    fHistCorrected[iPt]->Scale(1./fHistCorrected[iPt]->GetYaxis()->GetBinWidth(2));

    Plotter::Set2DHistAxes(fHistCorrected[iPt],"#Delta#eta","#Delta#varphi","#frac{d^{2}N}{d#Delta#eta d#Delta#varphi}",Form("h-%s, p_{T} bin: %d",particle.Data(),iPt));
    fHistCorrected[iPt]->Write();

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
