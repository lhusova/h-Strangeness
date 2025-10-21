#include "Plotter.h"
#include "Definitions.h"
void ScaleMixing(TH2F * histMix, Bool_t fitMaximum);

void MixingCorrection(TString particle="Hadron",TString region = "Signal",Int_t mult_i = 3, Bool_t fitMaximum =false, Bool_t correctMixing = true, Bool_t mcTrue = true, Bool_t integrMixing=false, Bool_t ptTrigDep = false, Bool_t ptAssoc =false){
  TString mixCorr = "_Mixing";
  if(!correctMixing) mixCorr ="_noMixing";
  if(mcTrue) mixCorr = "";
  TString mixingDivision ="";
  if(integrMixing) mixingDivision+="_integrMixing";
  if(!integrMixing) mixingDivision+= "_multMix";
  if(ptTrigDep) mixingDivision+="_pTtriggDep";
  if(ptAssoc)mixingDivision+="_pTAssocDep";
  // if(mcTrue||!correctMixing)mixingDivision = "";

  TFile * fFile = new TFile (Form("../data/Expanded/AnalysisResults_h%sgen%s_LHC25b4b_multEff.root",particle.Data(),mixCorr.Data()));
  
  THnF* same = (THnF*) fFile->Get(Form("Same_%s_%s_Expanded",particle.Data(),region.Data()));
  if(mult_i>0) same->GetAxis(5)->SetRange(mult_i,mult_i);
  
  TFile * fFileNew = TFile::Open (Form("../data/MixCorrected/%s/Corrected%s_MCgen%s_%s_%s_LHC25b4b_multEff.root",particle.Data(),mixCorr.Data(),mixingDivision.Data(),region.Data(),multiplicityNames[mult_i].Data()),"RECREATE");
  
  TFile * fFileMix;
  THnF* mix;
  TH2F * fHistMix;
  if(correctMixing&&!mcTrue) {
    mix = (THnF*) fFile->Get(Form("Mixed_%s_%s_Expanded",particle.Data(),region.Data()));
    if(!integrMixing && mult_i>0) mix->GetAxis(5)->SetRange(mult_i,mult_i);
    if(!ptTrigDep|| !ptAssoc){
      cout << "!ptTrigDep" << endl;
      fHistMix = (TH2F *) mix->Projection(0,1);
      fHistMix->SetName(Form("fHistMixed_%s",particle.Data()));
      fHistMix->RebinX(2);
      fHistMix->RebinY(2);
      ScaleMixing(fHistMix,fitMaximum);
    }
    
  }
  TH2F * fHistCorrected[nPtBins][nPtTriggBins];
  TH2F * fHistSame[nPtBins][nPtTriggBins];

  TH2F * fHistMixProjection[nPtBins][nPtTriggBins];
  TH2F * fHistSameProjection[nPtBins][nPtTriggBins][1];


  if(correctMixing&&!mcTrue&&!ptTrigDep) {
    fHistMix->Write();
    fHistMix->ProjectionX()->Write();
    }
  for (Int_t iPtTrigg = 0; iPtTrigg < nPtTriggBins; iPtTrigg++) {
    same->GetAxis(3)->SetRange(iPtTrigg+1,iPtTrigg+1);
    if(mcTrue) fFileMix = new TFile (Form("../data/Convolution_hh_PtTrigger%i.root",iPtTrigg));
    
    if(correctMixing&&ptTrigDep) {
      mix->GetAxis(3)->SetRange(iPtTrigg+1,iPtTrigg+1);
    } 
    for (Int_t iPt = 0; iPt < nPtBins; iPt++) {
      same->GetAxis(2)->SetRange(iPt+1,iPt+1);
      if(correctMixing&&ptAssoc) {
        mix->GetAxis(2)->SetRange(iPt+1,iPt+1);
      }

      fHistSame[iPt][iPtTrigg] = (TH2F *) same->Projection(0,1);
      fHistSame[iPt][iPtTrigg]->SetName(Form("fHistSame_%s_pt%d_ptTrigg%d",particle.Data(),iPt,iPtTrigg));
      fHistSame[iPt][iPtTrigg]->RebinX(2);
      fHistSame[iPt][iPtTrigg]->RebinY(2);
      fHistSame[iPt][iPtTrigg]->Scale(1./fHistSame[iPt][iPtTrigg]->GetXaxis()->GetBinWidth(2));
      fHistSame[iPt][iPtTrigg]->Scale(1./fHistSame[iPt][iPtTrigg]->GetYaxis()->GetBinWidth(2));
      fFileNew->cd();
      fHistSame[iPt][iPtTrigg]->Write();

      if(correctMixing&&(mcTrue||ptTrigDep||ptAssoc)){
        if(mcTrue) fHistMixProjection[iPt][iPtTrigg] = (TH2F *)fFileMix->Get(Form("fHistConvolution2d%s_pt%i_PtTriggerbin%i", particle.Data(),iPt,iPtTrigg));
        if(ptTrigDep||ptAssoc)fHistMixProjection[iPt][iPtTrigg] = (TH2F *)mix->Projection(0,1);
        fHistMixProjection[iPt][iPtTrigg]->SetName(Form("fHistMixed_%s_pt%d_ptTrigg%d",particle.Data(),iPt,iPtTrigg));
        fHistMixProjection[iPt][iPtTrigg]->RebinX(2);
        fHistMixProjection[iPt][iPtTrigg]->RebinY(2);
        ScaleMixing(fHistMixProjection[iPt][iPtTrigg],fitMaximum);
        fHistMixProjection[iPt][iPtTrigg]->Write();
      }
      // for (Int_t iPvz = 0; iPvz < 10; iPvz++) {
      //   same->GetAxis(4)->SetRange(iPvz+1,iPvz+1);
      //   mix->GetAxis(4)->SetRange(iPvz+1,iPvz+1);
        Int_t iPvz=0; 
        fHistSameProjection[iPt][iPtTrigg][iPvz] = (TH2F *) same->Projection(0,1);
        fHistSameProjection[iPt][iPtTrigg][iPvz]->SetName(Form("fHistSame_%s_pt%d%d_ptTrigg%d",particle.Data(),iPvz,iPt,iPtTrigg));
        fHistSameProjection[iPt][iPtTrigg][iPvz]->RebinX(2);
        fHistSameProjection[iPt][iPtTrigg][iPvz]->RebinY(2);
        fHistSameProjection[iPt][iPtTrigg][iPvz]->Write();
        if(correctMixing&&(mcTrue||ptTrigDep||ptAssoc))fHistSameProjection[iPt][iPtTrigg][iPvz]->Divide(fHistMixProjection[iPt][iPtTrigg]); 
        else if(correctMixing) fHistSameProjection[iPt][iPtTrigg][iPvz]->Divide(fHistMix);

        if(iPvz==0){
          fHistCorrected[iPt][iPtTrigg]=(TH2F*)fHistSameProjection[iPt][iPtTrigg][iPvz]->Clone();
          fHistCorrected[iPt][iPtTrigg]->SetName(Form("fHistCorrected_%s_pt%d_ptTrigg%d",particle.Data(),iPt,iPtTrigg));
        }else fHistCorrected[iPt][iPtTrigg]->Add(fHistSameProjection[iPt][iPtTrigg][iPvz]);

      // fHistCorrected[iPt][iPtTrigg]->Scale(1./fHistCorrected[iPt][iPtTrigg]->GetXaxis()->GetBinWidth(2));
      // fHistCorrected[iPt][iPtTrigg]->Scale(1./fHistCorrected[iPt][iPtTrigg]->GetYaxis()->GetBinWidth(2));

      Set2DHistAxes(fHistCorrected[iPt][iPtTrigg],"#Delta#eta","#Delta#varphi","#frac{d^{2}N}{d#Delta#eta d#Delta#varphi}",Form("h-%s, p_{T}^{assoc} bin: %d, p_{T}^{trigg} bin: %d",particle.Data(),iPt,iPtTrigg));
      fHistCorrected[iPt][iPtTrigg]->Write();
    }

  }
  // fFileNew->Close();
}
//____________________________________________________________________
void ScaleMixing(TH2F * histMix, Bool_t fitMaximum){
    Int_t nPhiBins = histMix->GetYaxis()->GetNbins();
    Int_t bin0 =histMix->GetXaxis()->FindBin(0.);
    Double_t scale=0.;
    if(fitMaximum){
      TF1 * leftFit = new TF1("leftFit","[0]*x+[1]",-1.5,0);
      TF1 * rightFit = new TF1("rightFit","[0]*x+[1]",0.,1.5);
      Double_t scaletmp = 0.;
      for(int i=1; i< nPhiBins+1; i++){
        histMix->GetYaxis()->SetRange(i,i);
        TH1F * tmp = (TH1F*)histMix->ProjectionX();
        tmp->SetName(Form("tmp%i",i));
        tmp->Fit(leftFit,"R");
        tmp->Fit(rightFit,"R");
        tmp->Write();
        rightFit->Write();
        leftFit->Write();
        scaletmp+=(leftFit->GetParameter(1)+rightFit->GetParameter(1))/2;
      }
      histMix->GetYaxis()->SetRange(1,nPhiBins);
      scale = scaletmp/nPhiBins;
    }else{
      for(Int_t iBinPhi=nPhiBins/2; iBinPhi<nPhiBins; iBinPhi++){
        scale+=histMix->GetBinContent(bin0,iBinPhi+1);
      }
      scale=scale/(nPhiBins/2);
    }
    if(scale!=0) histMix->Scale(1./scale);
    cout << "scale     "<< scale << endl;

};
