#include "Plotter.h"
#include "Definitions.h"
#include "TFile.h"

using namespace std;

TH1F * GetBackgroundHist(TH1F* hist,Bool_t fixBins = false, int nBins = 4);
TH1F * GetBackgroundHistVonMises(TH1F* hist, Int_t ptAssoc, Double_t parameters[6]);
TH1F * GetBackgroundLongRange(TH1F* hist,TH2F * hist2d);
TH1F * GetBackgroundFlow(TH1F* hist,Int_t part, Double_t meanPtTrigg, Double_t ptAssoc, TFile *  fileFlow, TFile *  fileFlow_v3,Bool_t include_v3);
Int_t GetNbinsFilled(TH1F * hist);
Int_t GetMinimumBin(TH1F* hist, Double_t minval); 
double CalculateMinimumError(TF1 * meanFunction,TF1 * minusFunction,TF1 * plusFunction,float minimum,float maximum);
double fitPedestalOnly(double *x,double *par) {
   double arg = 0;
   if ((x[0]>-1.&&x[0]<1.)||(x[0]>TMath::Pi()/2)) arg = 0;
   else arg =par[0];
   return arg;
}

void CalculateYield(Int_t part=2,Int_t ptTriggBin = 1, Int_t multClass = 1, Int_t bckg = 0){
  Bool_t doDensityInseadOfYield = kTRUE;
  TString desity = "";
  if(doDensityInseadOfYield){
    desity = "_density";
  }
  Bool_t doZYAMSyst = kFALSE;
  Int_t particleType =part;
  if(part==2)particleType=3;
  if(part==3)particleType=5;
  if(part==4)particleType=7;
  if(part==5)particleType=8;

  TFile *fileFlow = new TFile("../data/Flow/prapared_flow_v2.root");
  TFile *fileFlow_v3 = new TFile("../data/Flow/prapared_flow_v3.root");

  TFile * fFile[3];
  TFile * fFile2[3];

  Int_t nFile  = 1;
  if(part>0&&part<4) nFile = sizeof(InvMassRanges) / sizeof(TString);

  for (Int_t i = 0; i < nFile; i++) {
    fFile[i] = new TFile(Form("../data/MixCorrected/%s/Corrected_Mixing_integrMixing_%s_Corrected_LHC24ap_pass1_medium.root",nameSave[part].Data(),multiplicityNames[multClass].Data()));
    // fFile[i] = new TFile(Form("../FromKai/RootFiles/Shoving/MixCorrected/MixCorrectedMC_3.root"));
    // if(part>0&&part<4) fFile2[i] = new TFile(Form("../data/MixCorrected/%s/Corrected_Mixing_multMix_%s_%s_LHC24ap_pass1_medium.root",name[particleType+1].Data(),InvMassRanges[i].Data(),multiplicityNames[multClass].Data()));
  }
  TFile * fFileTrigger = new TFile("../data/DataFromTrain/AnalysisResults_hXi_raw_corrected_LHC24ap_pass1_medium.root");
  TH2F* histTriggers;
  if(part<2) histTriggers = (TH2F*) fFileTrigger->Get("correlate-strangeness_Trigg_1/sameEvent/TriggerParticlesV0");
  else if(part<4)histTriggers = (TH2F*) fFileTrigger->Get("h-strange-correlation/sameEvent/TriggerParticlesCascade");
  else if(part<5) histTriggers = (TH2F*) fFileTrigger->Get("correlate-strangeness/sameEvent/TriggerParticlesPion");
  else histTriggers = (TH2F*) fFileTrigger->Get("h-strange-correlation/sameEvent/TriggerParticlesHadron");
  
  if(multClass>0){
    // if(multClass==1)histTriggers->GetYaxis()->SetRange(multClass,multClass+1);
    // else 
    histTriggers->GetYaxis()->SetRange(multClass,multClass);
  }
  TH1F * hist1DTriggers = (TH1F *) histTriggers->ProjectionX();
  Float_t nTrigg=1;
  Double_t nTriggErr=1;
  // if(ptTriggBin<3)
  nTrigg=hist1DTriggers->IntegralAndError(hist1DTriggers->FindBin(ptTriggBins[ptTriggBin]),hist1DTriggers->FindBin(ptTriggBins[ptTriggBin+1]-0.1),nTriggErr);
  // else nTrigg= hist1DTriggers->Integral(hist1DTriggers->FindBin(ptTriggBins[ptTriggBin]),hist1DTriggers->FindBin(50-0.1)) + hist1DTriggers->GetBinContent(hist1DTriggers->GetXaxis()->GetNbins()+1);
  cout << nTrigg << endl;

  TH2F * fHistCorrected[nPtBins][nRegions];
  TH1F * fHistCorrectedProjection[nPtBins];
  TH1F * fHistBckg[nPtBins];
  Double_t yield_UE[nPtBins];
  Double_t yield_UE_error[nPtBins];
  Double_t yields[2][nPtBins];
  Double_t yields_err[2][nPtBins];
  cout << nPtBins << " _____ " << ptBins[nPtBins] << endl;
  TH1D * fHistNear = new TH1D("fHistNear","",nPtBins,ptBins);
  TH1D * fHistAway = new TH1D("fHistAway","",nPtBins,ptBins);
  TH1D * fHistUE = new TH1D("fHistUE","",nPtBins,ptBins);
  TH1F * fHistTrigg;
  Int_t minimalPtBin =0;
  int nBinsForZYAM = 5;
  TFile * fFileNew = TFile::Open (Form("../data/Yields/%s/Yields%s_Corrected_LHC24ap_pass1_medium_%s_%s_ptTrigg%d.root",nameSave[part].Data(),desity.Data(),multiplicityNames[multClass].Data(),nameBackground[bckg].Data(),ptTriggBin),"RECREATE");
  hist1DTriggers->Write();

  Double_t parameters[6]={ 1.9,0.0299079,10,0.02,0.,0.3};


  for (Int_t iPt = minimalPtBin; iPt < nPtBins; iPt++) {
    for (Int_t iFile = 0; iFile < nFile; iFile++) { // side-band subtraction
      fHistCorrected[iPt][iFile] = (TH2F*)fFile[iFile]->Get(Form("fHistCorrected_%s_pt%d_ptTrigg%d",nameSave[part].Data(),iPt,ptTriggBin));
      // if(part>0&&part<4)fHistCorrected[iPt][iFile]->Add ((TH2F*)fFile2[iFile]->Get(Form("fHistCorrected_%s_pt%d_ptTrigg%d",name[particleType+1].Data(),iPt,ptTriggBin)));
      // if(iFile>0) fHistCorrected[iPt][0]->Add(fHistCorrected[iPt][iFile],-1);
    }
    // if(bckg==1&&iPt==minimalPtBin){
    //   fHistUE->SetBinContent(iPt+1,0.001/(ptBins[iPt+1]-ptBins[iPt]));
    //   fHistNear->SetBinContent(iPt+1,0.001/(ptBins[iPt+1]-ptBins[iPt]));
    //   fHistAway->SetBinContent(iPt+1,0.001/(ptBins[iPt+1]-ptBins[iPt]));
    //   continue;
    // }
    fHistCorrected[iPt][0]->Write();
    // fHistCorrected[iPt][0]->GetXaxis()->SetRangeUser(-1.1,1.1);
    // cout << "bin content :" << fHistCorrected[iPt][0]->GetBinContent(20,10)<<" err 2d " << fHistCorrected[iPt][0]->GetBinError(20,10) << endl;
    // fHistCorrected[iPt][0]->Scale(1./nTrigg);
    // cout << "bin content after scaling :" << fHistCorrected[iPt][0]->GetBinContent(20,10) << " err 2d after scaling " << fHistCorrected[iPt][0]->GetBinError(20,10) << endl;

    fHistCorrectedProjection[iPt] = (TH1F *) fHistCorrected[iPt][0]->ProjectionY(Form("phiProj_withBckg_pT%d",iPt),5,21,"e");
    if(!doDensityInseadOfYield)
      fHistCorrectedProjection[iPt]->Scale(1./fHistCorrectedProjection[iPt]->GetXaxis()->GetBinWidth(2));
    if(iPt==minimalPtBin){
      fHistTrigg = (TH1F *) fHistCorrectedProjection[iPt]->Clone();
      fHistTrigg->SetName("fHistTrigg");
      for (int i = 1; i < fHistTrigg->GetXaxis()->GetNbins()+1; i++) {
        fHistTrigg->SetBinContent(i,nTrigg);
        fHistTrigg->SetBinError(i,nTriggErr);
      }
      fHistTrigg->Write();
    }
    fHistCorrectedProjection[iPt]->Divide(fHistTrigg);
    fHistCorrectedProjection[iPt]->SetMarkerStyle(20);
    if(doDensityInseadOfYield){
      float etaWidth = fHistCorrected[iPt][0]->GetXaxis()->GetBinUpEdge(21)-fHistCorrected[iPt][0]->GetXaxis()->GetBinLowEdge(5);
      fHistCorrectedProjection[iPt]->Scale(1./etaWidth);
    }
    // fHistCorrectedProjection[iPt]->Scale(fHistCorrected[iPt][0]->GetYaxis()->GetBinWidth(2));


    SetHistAxes(fHistCorrectedProjection[iPt],"#Delta#varphi","#frac{1}{N_{trigg}} #frac{dN}{d#Delta#varphi}");

    // fHistCorrectedProjection[iPt]->SetName(Form("phiProj_withBckg_pT%d",iPt));
    fHistCorrectedProjection[iPt]->Write();
    if(doZYAMSyst){
      if(iPt<3)nBinsForZYAM=2;
      else if(iPt>2&&iPt<5)nBinsForZYAM=3;
      else if(iPt>4&&iPt<8)nBinsForZYAM=5;
      else if(iPt>=8)nBinsForZYAM=6;
  }
    if(bckg==0)fHistBckg[iPt] = GetBackgroundHist(fHistCorrectedProjection[iPt],false,nBinsForZYAM);
    else if(bckg==1)fHistBckg[iPt] = GetBackgroundLongRange(fHistCorrectedProjection[iPt],fHistCorrected[iPt][0]);
    else if (bckg==2)fHistBckg[iPt] = GetBackgroundFlow(fHistCorrectedProjection[iPt],part,hist1DTriggers->GetMean(),(ptBins[iPt]+ptBins[iPt+1])/2,fileFlow,fileFlow_v3,kFALSE);
    else if (bckg==3)fHistBckg[iPt] = GetBackgroundFlow(fHistCorrectedProjection[iPt],part,hist1DTriggers->GetMean(),(ptBins[iPt]+ptBins[iPt+1])/2,fileFlow,fileFlow_v3,kTRUE);
    else if(bckg==4){
      fHistBckg[iPt] = GetBackgroundHistVonMises(fHistCorrectedProjection[iPt],iPt,parameters);
    }
    else {
      cout << "ERROR: Underlying Event not defined!!!!" << endl;
      return;
    }
    fHistBckg[iPt]->SetName(Form("underlying_ev_pT%d",iPt));
    fHistBckg[iPt]->SetMarkerStyle(20);
    fHistCorrectedProjection[iPt]->Add(fHistBckg[iPt],-1);
    fHistBckg[iPt]->Write();

    // if(doDensityInseadOfYield){
    //   float phWidth = fHistBckg[iPt]->GetXaxis()->GetBinUpEdge(36)-fHistBckg[iPt]->GetXaxis()->GetBinLowEdge(1);
    //   cout << " phi width " << phWidth << endl;
    //   fHistBckg[iPt]->Scale(1./phWidth);
    // }

    if(bckg==1) {
      yield_UE[iPt]=fHistBckg[iPt]->IntegralAndError(1,18,yield_UE_error[iPt],"width")*2;
      yield_UE[iPt] = yield_UE[iPt] / GetNbinsFilled(fHistBckg[iPt]);
      yield_UE[iPt] = yield_UE[iPt] * 18;
    }
    else {
      yield_UE[iPt] = 0;
      yield_UE_error[iPt]=0;
      for (int i=1; i<fHistBckg[iPt]->GetXaxis()->GetNbins()+1; i++) {
        if(doDensityInseadOfYield){
          if(i< fHistCorrectedProjection[iPt]->FindBin(1.1))
            continue;
          if(i>fHistCorrectedProjection[iPt]->FindBin(1.7))
            continue;
          yield_UE[iPt]+=fHistBckg[iPt]->GetBinContent(i) + fHistCorrectedProjection[iPt]->GetBinContent(i);
          yield_UE_error[iPt]=TMath::Sqrt(yield_UE_error[iPt]*yield_UE_error[iPt]+fHistBckg[iPt]->GetBinError(i)*fHistBckg[iPt]->GetBinError(i)+2*fHistBckg[iPt]->GetBinError(i)*yield_UE_error[iPt]);
        }else{
          yield_UE[iPt]+=fHistBckg[iPt]->GetBinContent(i)*fHistBckg[iPt]->GetBinWidth(i);
          yield_UE_error[iPt]=TMath::Sqrt(yield_UE_error[iPt]*yield_UE_error[iPt]+fHistBckg[iPt]->GetBinError(i)*fHistBckg[iPt]->GetBinError(i)+2*fHistBckg[iPt]->GetBinError(i)*yield_UE_error[iPt]);
        }
      }
      cout << yield_UE[iPt] << endl;
      if(doDensityInseadOfYield){
        float phWidth = fHistBckg[iPt]->GetXaxis()->GetBinUpEdge(fHistBckg[iPt]->FindBin(1.7))-fHistBckg[iPt]->GetXaxis()->GetBinLowEdge(fHistBckg[iPt]->FindBin(1.1));
        yield_UE[iPt]=yield_UE[iPt]/phWidth;
        yield_UE_error[iPt]=yield_UE_error[iPt]/phWidth;
      }else{
        yield_UE_error[iPt]=yield_UE_error[iPt]*fHistBckg[iPt]->GetBinWidth(1);
      }
      cout << yield_UE_error[iPt] << endl;
    }
    //yield_UE[iPt]=fHistBckg[iPt]->IntegralAndError(1,36,yield_UE_error[iPt],"width");//(fHistBckg[iPt]->GetXaxis()->GetFirst(),fHistBckg[iPt]->GetXaxis()->GetLast(),yield_UE_error[iPt],"width");
    yield_UE[iPt]=yield_UE[iPt]/50;
    yield_UE_error[iPt]=yield_UE_error[iPt]/50;
    cout << (ptBins[iPt+1]-ptBins[iPt]) << endl;
    fHistUE->SetBinContent(iPt+1,yield_UE[iPt]/(ptBins[iPt+1]-ptBins[iPt]));
    if(bckg==1)fHistUE->SetBinError(iPt+1,(yield_UE_error[iPt]*2)/(ptBins[iPt+1]-ptBins[iPt]));
    else fHistUE->SetBinError(iPt+1,yield_UE_error[iPt]/(ptBins[iPt+1]-ptBins[iPt]));

    
    fHistCorrectedProjection[iPt]->SetName(Form("phiProj_noBckg_pT%d",iPt));
    fHistCorrectedProjection[iPt]->Write();

    if(doDensityInseadOfYield){
      float phWidth = fHistCorrectedProjection[iPt]->GetXaxis()->GetBinUpEdge(fHistCorrectedProjection[iPt]->FindBin(1))-fHistCorrectedProjection[iPt]->GetXaxis()->GetBinLowEdge(fHistCorrectedProjection[iPt]->FindBin(-1));
      yields[0][iPt]=fHistCorrectedProjection[iPt]->IntegralAndError(fHistCorrectedProjection[iPt]->FindBin(-1),fHistCorrectedProjection[iPt]->FindBin(1),yields_err[0][iPt],"")/phWidth;
    }else
      yields[0][iPt]=fHistCorrectedProjection[iPt]->IntegralAndError(fHistCorrectedProjection[iPt]->FindBin(-1),fHistCorrectedProjection[iPt]->FindBin(1),yields_err[0][iPt],"width");
    yields[1][iPt]=fHistCorrectedProjection[iPt]->IntegralAndError(fHistCorrectedProjection[iPt]->FindBin(TMath::Pi()-TMath::Pi()/2),fHistCorrectedProjection[iPt]->FindBin(TMath::Pi()+TMath::Pi()/2),yields_err[1][iPt],"width");
    fHistNear->SetBinContent(iPt+1,yields[0][iPt]/(ptBins[iPt+1]-ptBins[iPt]));
    fHistNear->SetBinError(iPt+1,yields_err[0][iPt]/(ptBins[iPt+1]-ptBins[iPt]));
    fHistAway->SetBinContent(iPt+1,yields[1][iPt]/(ptBins[iPt+1]-ptBins[iPt]));
    fHistAway->SetBinError(iPt+1,yields_err[1][iPt]/(ptBins[iPt+1]-ptBins[iPt]));
  }

  if(doDensityInseadOfYield)
    SetHistAxes(fHistNear,"#font[52]{p}^{assoc}_{T} (GeV/#font[52]{c})","1/#font[52]{N}_{Trigg}1/#Delta#eta#Delta#varphi d#font[52]{N}/d#font[52]{p}_{T}");
  else 
    SetHistAxes(fHistNear,"#font[52]{p}^{assoc}_{T} (GeV/#font[52]{c})","1/#font[52]{N}_{Trigg}d#font[52]{N}/d#font[52]{p}_{T}");
  SetHist(fHistNear,"",20,kRed+1,1.);

  SetHistAxes(fHistAway,"#font[52]{p}^{assoc}_{T} (GeV/#font[52]{c})","1/#font[52]{N}_{Trigg}d#font[52]{N}/d#font[52]{p}_{T}");
  SetHist(fHistAway,"",20,kBlue+1,1.);

  if(doDensityInseadOfYield)
    SetHistAxes(fHistUE,"#font[52]{p}^{assoc}_{T} (GeV/#font[52]{c})","1/#font[52]{N}_{Trigg}1/#Delta#eta#Delta#varphi d#font[52]{N}/d#font[52]{p}_{T}");
  else  
    SetHistAxes(fHistUE,"#font[52]{p}^{assoc}_{T} (GeV/#font[52]{c})","1/#font[52]{N}_{Trigg}d#font[52]{N}/d#font[52]{p}_{T}");
  SetHist(fHistUE,"",20,kGreen+1,1.);

  TCanvas *can = CreateCanvas("c");
  gStyle->SetErrorX(0.01);
  // can->GetPadSave()->SetLogy();
  fHistNear->GetYaxis()->SetRangeUser(0.5*fHistUE->GetMinimum(),1.5*fHistNear->GetMaximum());
  fHistNear->DrawCopy("");
  fHistAway->DrawCopy("p same");
  fHistUE->DrawCopy("p same");

  hist1DTriggers->Write();
  fHistNear->Write();
  fHistAway->Write();
  fHistUE->Write();

  TPaveText *pave = new TPaveText();
  SetPaveText(pave,42,0.05, 0, 0,33,0,0.55,0.95, 0.65,0.95);
  pave->AddText("ALICE, Work in Progress");
  pave->AddText("pp, 5.36 TeV");
  pave->AddText(Form("%s",multiplicityPave[multClass].Data()));
  pave->AddText(Form("h-%s",finalNames[part].Data()));
  pave->AddText(Form("%g < #font[52]{p}^{trigg}_{T} < %g GeV/#font[52]{c}",ptTriggBins[ptTriggBin],ptTriggBins[ptTriggBin+1]));
  pave->AddText("|#Delta#eta| < 1.1");
  pave->Draw("same");

  TLegend *leg = CreateLegend(0.15, 0.45, 0.15, 0.45,0.05);
  leg->AddEntry(fHistNear,"Near-side, |#Delta#varphi|<#pi/2","pl");
  leg->AddEntry(fHistAway,"Away-side, |#Delta#varphi-#pi|<#pi/2","pl");
  leg->AddEntry(fHistUE,"Underlying event #times 1/50","pl");
  leg->Draw("same");

  // can->SaveAs(Form("../Plots/Yields/%s/Yield_ppRefEff_MultMixing_%s_%s_fullrangePeak_11_%s_ptTrigg%d.pdf",nameSave[part].Data(),nameSave[part].Data(),multiplicityNames[multClass].Data(),nameBackground[bckg].Data(),ptTriggBin));
  // fFileNew->Close();
}
//____________________________________________________________________
TH1F * GetBackgroundHist(TH1F* hist, Bool_t fixBins, Int_t nBins){

  // Int_t bins[8]={1,2,3,4,33,34,35,36};
  Int_t bins[4]={1,18,19,36};
  Double_t value =0;
  Double_t err=0;
  Double_t minimum =0;
  cout<< "  nBins ZYAM " <<nBins << endl;
  for (size_t i = 0; i < nBins; i++) {
    if(fixBins){
      value+=hist->GetBinContent(bins[i]);
      err+=TMath::Power(hist->GetBinError(bins[i]),2);
    }else{
      minimum = hist->GetMinimum(minimum);
      value+= minimum ;//hist->GetBinContent(bins[i]);
      cout <<"minmum  " << GetMinimumBin(hist,minimum) << ": " <<minimum<< " err___ " << hist->GetBinError(GetMinimumBin(hist,minimum)) << endl;
      err+=TMath::Power(hist->GetBinError(GetMinimumBin(hist,minimum)),2);
    }
    
  }
  err=TMath::Sqrt(err);
  value=value/nBins;
  err=err/nBins;
  cout << err << endl;
  cout << value << endl;

  TH1F* bckg = (TH1F*) hist->Clone();
  for (size_t iPhi = 1; iPhi < hist->GetXaxis()->GetNbins()+1; iPhi++) {
    bckg->SetBinContent(iPhi,value);
    bckg->SetBinError(iPhi,err);
  }
  return bckg;

}
//____________________________________________________________________
TH1F * GetBackgroundLongRange(TH1F* hist,TH2F * hist2d){

  hist2d->GetXaxis()->SetRange(1,4);
  TH1F * projLeft = (TH1F*)hist2d->ProjectionY();

  hist2d->GetXaxis()->SetRange(22,25);
  TH1F * projRight = (TH1F*)hist2d->ProjectionY();

  TH1F * histBckg = (TH1F*) projLeft->Clone();
  histBckg->Add(projRight);

  // TH1F * fHistZyam = GetBackgroundHist(hist);
  double scale = 0;
  double scaleDenom = 0;
  Int_t filled = 0;
  for (Int_t b = 19; b <= 36; ++b) {
    if (histBckg->GetBinContent(b) != 0) {
      scaleDenom += histBckg->GetBinContent(b);
      scale+= hist->GetBinContent(b);
      filled++;
    }
  }
  if (filled > 0) scale = scale / scaleDenom;
  else scale = 1;

  histBckg->Scale(scale);
  // for (size_t i = 16; i < histBckg->GetXaxis()->GetNbins(); i++) {
  //   histBckg->SetBinContent(i+1,fHistZyam->GetBinContent(i+1));
  //   histBckg->SetBinError(i+1,fHistZyam->GetBinError(i+1));
  // }
   // histBckg->Scale(hist->GetBinContent(17)/histBckg->GetBinContent(17));

  return histBckg;
}
//____________________________________________________________________
TH1F * GetBackgroundFlow(TH1F* hist, Int_t part, Double_t meanPtTrigg, Double_t ptAssoc, TFile *  fileFlow, TFile *  fileFlow_v3, Bool_t include_v3){

  // Int_t bins[8]={1,2,3,4,33,34,35,36};
  Int_t bins[6]={1,2,3,16,17,18};
  Double_t value =0;

  for (size_t i = 0; i < 6; i++) {
    value+=hist->GetBinContent(bins[i]);
  }
  value=value/6;

  TGraphErrors * trigg = (TGraphErrors * ) fileFlow->Get("v2_charged");
  TGraphErrors * assoc;
  if(part==0) assoc = (TGraphErrors * ) fileFlow->Get("v2_K0_Combined");
  else if(part==1) assoc = (TGraphErrors * ) fileFlow->Get("v2_Lambda_Combined");
  else if(part==4) assoc = (TGraphErrors * ) fileFlow->Get("v2_Pion_Combined");
  else {
    cout << "ERROR!! Cascade v2 not available!!! ";
    return 0x0;
  }

  TF1 * bcgkFunction;
  if(include_v3){
    TGraphErrors * trigg_v3 = (TGraphErrors * ) fileFlow_v3->Get("v3_charged");
    TGraphErrors * assoc_3;
    if(part==0) assoc_3 = (TGraphErrors * ) fileFlow_v3->Get("v3_K0_Combined");
    else if(part==1) assoc_3 = (TGraphErrors * ) fileFlow_v3->Get("v3_Lambda_Combined");
    else if(part==4) assoc_3 = (TGraphErrors * ) fileFlow_v3->Get("v3_Pion_Combined");
    else {
      cout << "ERROR!! Cascade v2 not available!!! ";
      return 0x0;
    }

    bcgkFunction = new TF1("bcgkFunction","[0]*(1+2*[1]*([2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)))",-TMath::Pi()/2,3/2*TMath::Pi());
    bcgkFunction->SetParameter(0,value);

    bcgkFunction->SetParameter(1,trigg->Eval(meanPtTrigg));
    cout << assoc->Eval(ptAssoc) << endl;
    bcgkFunction->SetParameter(2,assoc->Eval(ptAssoc));
    bcgkFunction->SetParameter(3,assoc_3->Eval(ptAssoc));

  }else{
    bcgkFunction = new TF1("bcgkFunction","[0]*(1+2*[1]*[2]*TMath::Cos(2*x))",-TMath::Pi()/2,3/2*TMath::Pi());
    cout << value << endl;
    cout << trigg->Eval(meanPtTrigg) << endl;
    bcgkFunction->SetParameter(0,value);

    bcgkFunction->SetParameter(1,trigg->Eval(meanPtTrigg));
    cout << assoc->Eval(ptAssoc) << endl;
    bcgkFunction->SetParameter(2,assoc->Eval(ptAssoc));
  }


  TH1F * bckg = (TH1F *)hist->Clone();

  for (Int_t i = 0; i < bckg->GetXaxis()->GetNbins(); i++) {
    bckg->SetBinContent(i+1,bcgkFunction->Eval(bckg->GetBinCenter(i+1)));
  }
  return bckg;

}
Int_t GetNbinsFilled(TH1F * hist) {
  Int_t count = 0;
  for (Int_t i = 1; i <= 18; ++i) {
    if (hist->GetBinContent(i) != 0) {
      count++;
    }
  }
  cout << count << endl;
  return count;
}
TH1F * GetBackgroundHistVonMises(TH1F* hist, Int_t ptAssoc, Double_t parameters[6]){
  cout << "--------- iPt " << ptAssoc<<endl;
  TF1 * fitPedestal = new TF1("fitPedestal",fitPedestalOnly,-TMath::Pi()/2,3*TMath::Pi()/2,1);
  TCanvas * can = new TCanvas(Form("canVonMises_pT%d",ptAssoc));
  
  TH1F * histClone = (TH1F *) hist->Clone();
  histClone->SetName(Form("histClone_pT%d",ptAssoc));
  histClone->DrawCopy("e");
  // histClone->Fit(fitPedestal,"R");
  // if(parameters[5]<0)parameters[5]=0.5;
  if(parameters[2]<0)parameters[2]=9.;
  if(parameters[1]<0)parameters[1]=1.;
  if(ptAssoc==3){
    parameters[0]=0.00746312;
    parameters[1]=0.0104298;
    parameters[2]=10;
    parameters[3]=0.03;
    parameters[5]=0.15;
  }
  TF1 * funcMises = new TF1("funcMises", "[0]+[1]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x))+[3]*TMath::Gaus(x,[4],[5])",-TMath::Pi()/2,TMath::Pi()/2);//+ [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))
  funcMises->SetParNames("pedestal","norm_near","kappa_near","norm Gauss","mean","sigma");//,"norm_away","kappa_away");
  funcMises->SetParameter(0,parameters[0]);
  // funcMises->SetParLimits(0,0.,1e10);
  funcMises->SetParameter(1,parameters[1]);
  funcMises->SetParameter(2,parameters[2]);
  funcMises->SetParameter(3,parameters[3]);
  // funcMises->SetParLimits(3,-1,1e10);
  funcMises->FixParameter(4,0);
  funcMises->SetParameter(5,0.5);
  // funcMises->SetParameter(0,0);
  funcMises->SetLineColor(kGray+2);
  cout << "Fitting Von Mises function to the data..." << endl;
  
  histClone->Fit(funcMises,"R");
  // if(funcMises->GetParameter(2)<=0){
  //   // while(funcMises->GetParameter(2)<=0){
  //     funcMises->SetParameter(1,2);
  //     funcMises->SetParameter(2,5);
  //     funcMises->SetParameter(3,0.03);
  //     histClone->Fit(funcMises,"R");
  //   // }
  // }
  // if(std::abs(funcMises->GetParameter(5))<0.3||funcMises->GetParameter(3)>1000){
  // //   while(funcMises->GetParameter(5)<=0||funcMises->GetParameter(3)>1000){
  //   // if(ptAssoc<3){
  //     // funcMises->SetParameter(1,2);
  //     // funcMises->SetParameter(2,5);
  //   // }
  // if(!(ptAssoc==6||ptAssoc==8||ptAssoc==9)){
      // funcMises->SetParameter(5,0.3);
  //     funcMises->SetParameter(2,40);
  // // }
      // funcMises->SetParameter(3,0.02);
      TFitResultPtr frp(0);
      frp = histClone->Fit(funcMises,"RS");
  //   }
      // TFitResult * ft =  pt.Get();
      // if ( frp >= 0 ) {
		    // frp->Print();
      // }else    
  //     if(ptAssoc<2)  {
      if(!frp->IsValid()){//||std::abs(funcMises->GetParameter(5))<0.1||funcMises->GetParameter(3)<0||funcMises->GetParameter(3)>1000||funcMises->GetParameter(2)<0){//||funcMises->GetParameter(3)<0||funcMises->GetParameter(3)>1000||funcMises->GetParameter(2)<0){//||std::abs(funcMises->GetParameter(5))<0.1||funcMises->GetParameter(3)<0){
        while(!frp->IsValid()){//||std::abs(funcMises->GetParameter(5))<0.1||funcMises->GetParameter(3)>1000||funcMises->GetParameter(2)<0){//||funcMises->GetParameter(3)<0||funcMises->GetParameter(3)>1000||funcMises->GetParameter(2)<0){//||std::abs(funcMises->GetParameter(5))<0.2||funcMises->GetParameter(3)<0||funcMises->GetParameter(0)<0){
  //     // // // //   if(ptAssoc<8)  {
  //     // // //   
          // if(ptAssoc==2)
          funcMises->SetParameter(1,0.05);
          funcMises->SetParameter(2,25);
  //     // //     // if(ptAssoc==2)
      funcMises->SetParameter(5,0.3);
  // //     if(ptAssoc<6){
          funcMises->SetParameter(3,0.06);
          frp = histClone->Fit(funcMises,"RS");
  //       }
      }
    }
  
  // if((funcMises->GetParameter(3)<0||funcMises->GetParameter(3)>1000)){
  //   while(funcMises->GetParameter(3)<0||funcMises->GetParameter(3)>1000||!frp->IsValid()){
  //     funcMises->SetParameter(2,50);
  //     funcMises->SetParameter(1,0.01);
  //     funcMises->SetParameter(3,0.02);
  //     // funcMises->SetParameter(5,0.1);
  //     frp = histClone->Fit(funcMises,"RS");
  //   }
  // }
  //     if(ptAssoc<2&&std::abs(funcMises->GetParameter(5))>100){
  //   while(std::abs(funcMises->GetParameter(5))>100){
  //     funcMises->SetParameter(2,0.1);
  //     // funcMises->SetParameter(3,0.02);
  //     funcMises->SetParameter(5,0.1);
  //     frp = histClone->Fit(funcMises,"RS");
  //   }
  // }
  TF1 * pedestalFunction = new TF1("pedestalFunction","[0]",-TMath::Pi()/2,TMath::Pi()/2);
  pedestalFunction->SetParameter(0,funcMises->GetParameter(0));
  pedestalFunction->SetLineColor(kGreen+1);
  pedestalFunction->Draw("same");
  TF1 * nearFunction = new TF1("nearFunction","[0]/(2*TMath::Pi()*TMath::BesselI0([1]))*TMath::Exp([1]*TMath::Cos(x))",-TMath::Pi()/2,TMath::Pi()/2);
  nearFunction->SetParameter(0,funcMises->GetParameter(1));
  nearFunction->SetParameter(1,funcMises->GetParameter(2));
  nearFunction->SetLineColor(kRed+1);
  nearFunction->Draw("same");
  TF1 * nearGauss = new TF1("nearGauss","[0]*TMath::Gaus(x,[1],[2])",-TMath::Pi()/2,TMath::Pi()/2);
  nearGauss->SetParameter(0,funcMises->GetParameter(3));
  nearGauss->SetParameter(1,funcMises->GetParameter(4));
  nearGauss->SetParameter(2,funcMises->GetParameter(5));
  nearGauss->SetLineColor(kBlue+1);
  nearGauss->Draw("same");
  TF1 * awayFunction = new TF1("awayFunction","[0]/(2*TMath::Pi()*TMath::BesselI0([1]))*TMath::Exp([1]*TMath::Cos(x))",-TMath::Pi()/2,3*TMath::Pi()/2);
  awayFunction->SetParameter(0,funcMises->GetParameter(3));
  awayFunction->SetParameter(1,funcMises->GetParameter(4));
  awayFunction->SetLineColor(kBlue+1);
  // awayFunction->Draw("same");
  funcMises->Write();
  TF1 * funcMisesErrPos = new TF1("funcMisesErrPos", "[0]+[1]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x))+[3]*TMath::Gaus(x,[4],[5])",-TMath::Pi()/2,TMath::Pi()/2);//+ [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))
  funcMisesErrPos->FixParameter(0,funcMises->GetParameter(0)+funcMises->GetParError(0));
  funcMisesErrPos->FixParameter(1,funcMises->GetParameter(1)+funcMises->GetParError(1));
  funcMisesErrPos->FixParameter(2,funcMises->GetParameter(2)+funcMises->GetParError(2));
  funcMisesErrPos->FixParameter(3,funcMises->GetParameter(3)+funcMises->GetParError(3));
  funcMisesErrPos->FixParameter(4,funcMises->GetParameter(4)+funcMises->GetParError(4));
  funcMisesErrPos->FixParameter(5,funcMises->GetParameter(5)+funcMises->GetParError(5));
  funcMisesErrPos->Write();
  TF1 * funcMisesErrNeg = new TF1("funcMisesErrNeg", "[0]+[1]/(2*TMath::Pi()*TMath::BesselI0([2]))*TMath::Exp([2]*TMath::Cos(x))+[3]*TMath::Gaus(x,[4],[5])",-TMath::Pi()/2,TMath::Pi()/2);//+ [3]/(2*TMath::Pi()*TMath::BesselI0([4]))*TMath::Exp([4]*TMath::Cos(x))
  funcMisesErrNeg->FixParameter(0,funcMises->GetParameter(0)-funcMises->GetParError(0));
  funcMisesErrNeg->FixParameter(1,funcMises->GetParameter(1)-funcMises->GetParError(1));
  funcMisesErrNeg->FixParameter(2,funcMises->GetParameter(2)-funcMises->GetParError(2));
  funcMisesErrNeg->FixParameter(3,funcMises->GetParameter(3)-funcMises->GetParError(3));
  funcMisesErrNeg->FixParameter(4,funcMises->GetParameter(4)-funcMises->GetParError(4));
  funcMisesErrNeg->FixParameter(5,funcMises->GetParameter(5)-funcMises->GetParError(5));
  funcMisesErrNeg->Write();
  double err1 = CalculateMinimumError(funcMises,funcMisesErrNeg,funcMisesErrPos,-TMath::Pi()/2,-1);
  double err2 = CalculateMinimumError(funcMises,funcMisesErrNeg,funcMisesErrPos,1,TMath::Pi()/2);
  TH1F* bckg = (TH1F*) hist->Clone();
  cout<< (funcMises->GetMinimum(-TMath::Pi()/2,-1)+funcMises->GetMinimum(1,TMath::Pi()/2))/2 <<  endl;
  cout<< "Error.  " << TMath::Sqrt(err1*err1+err2*err2)/2 <<  endl;
  cout<< "Prew.  " << funcMises->GetParError(TMath::Sqrt(err1*err1+err2*err2)/2) <<  endl;
  for (size_t iPhi = 1; iPhi < hist->GetXaxis()->GetNbins()+1; iPhi++) {
    bckg->SetBinContent(iPhi,(funcMises->GetMinimum(-TMath::Pi()/2,-1)+funcMises->GetMinimum(1,TMath::Pi()/2))/2);
    bckg->SetBinError(iPhi,TMath::Sqrt(err1*err1+err2*err2)/2);
  }
  parameters[0]=funcMises->GetParameter(0);
  parameters[1]=funcMises->GetParameter(1);
  parameters[2]=funcMises->GetParameter(2);
  parameters[3]=funcMises->GetParameter(3);
  parameters[5]=funcMises->GetParameter(5);  
  return bckg;

}
double CalculateMinimumError(TF1 * meanFunction,TF1 * minusFunction,TF1 * plusFunction,float minimum,float maximum){
  double minMean = meanFunction->GetMinimum(minimum, maximum);
  double minMeanX = meanFunction->GetMinimumX(minimum, maximum);

  // Find minimum values for the error-shifted functions
  double minMinus = minusFunction->Eval(minMeanX);
  double minPlus  = plusFunction->Eval(minMeanX);

  // Calculate the absolute differences
  double errMinus = fabs(minMean - minMinus);
  double errPlus  = fabs(minMean - minPlus);

  // Return the largest difference as the error
  return std::max(errMinus, errPlus);
}

Int_t GetMinimumBin(TH1F* hist, Double_t minval) 
{ 
  Int_t bin, binx;
  Int_t locm;
  Int_t xfirst  = hist->GetXaxis()->GetFirst();
  Int_t xlast   = hist->GetXaxis()->GetLast();
  Double_t minimum, value;
  minimum = FLT_MAX;
  for (binx=xfirst;binx<=xlast;binx++) {
    value = hist->GetBinContent(binx);
    if (value < minimum && value>minval) {
      minimum = value;
      bin  = binx;
    }
        
   }
   return bin;
}
