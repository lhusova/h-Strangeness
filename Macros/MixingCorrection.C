#include "Plotter.h"
void ScaleMixing(TH2F * histMix);

void MixingCorrection(TString particle="OmegaPlus"){
  TFile * fFile = new TFile ("../data/AnalysisResults_Hyperloop_16_04.root");

  TH3F* same = (TH3F*) fFile->Get(Form("correlate-strangeness/sameEvent/Hadron%s",particle.Data()));
  TH3F* mix = (TH3F*) fFile->Get(Form("correlate-strangeness/mixedEvent/Hadron%s",particle.Data()));

  const Int_t nPtBins = 10;
  TH2F * fHistCorrected[nPtBins];
  TH2F * fHistMixProjection[nPtBins];

  TFile * fFileNew = TFile::Open (Form("../data/MixCorrected_%s.root",particle.Data()),"RECREATE");
  for (Int_t iPt = 0; iPt < nPtBins; iPt++) {

    same->GetZaxis()->SetRange(iPt+1,iPt+1);
    // mix->GetZaxis()->SetRange(iPt+1,iPt+1);

    fHistCorrected[iPt] = (TH2F *) same->Project3D("xy");
    fHistCorrected[iPt]->SetName(Form("fHistCorrected_%s_pt%d",particle.Data(),iPt));
    fHistCorrected[iPt]->RebinX(2);
    fHistCorrected[iPt]->RebinY(2);

    fHistMixProjection[iPt] = (TH2F *) mix->Project3D("xy");
    fHistMixProjection[iPt]->SetName(Form("fHistMix_%s_pt%d",particle.Data(),iPt));

    fHistMixProjection[iPt]->RebinX(2);
    fHistMixProjection[iPt]->RebinY(2);

    ScaleMixing(fHistMixProjection[iPt]);

    fHistCorrected[iPt]->Divide(fHistMixProjection[iPt]);
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
