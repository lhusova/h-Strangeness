#include "Riostream.h"
#include "TH1F.h"

void ErrRatioCorr(TH1F* hNum, TH1F* hDenom, TH1F* hRatio, Int_t FullCorr){
  //FullCorr == 1 means ro = 1; //full correlation                                                           
  //FullCorr == 0 means ro = sigmaDenom/sigmaNum ; //data sample at numerator is a subset of denominator     
  //FullCorr == 2 :another possibility (followed by Fiorella's according to Roger Barlow): I think this is an approximation of Barlow's prescription, and my procedure is better. I found out that errors calculated in this way vary within +-10% from those calculated by me.

  Float_t Err1=0;
  Float_t Err2=0;
  Float_t ErrC=0;
  Float_t Err=0;

  for (Int_t b=1; b<=hNum->GetNbinsX();b++){
    if (hNum->GetBinContent(b)==0 ||hDenom->GetBinContent(b)==0){
      hRatio->SetBinError(b,0);
      continue;
    }
    Err1=pow(hNum->GetBinError(b)/hNum->GetBinContent(b),2);
    Err2=pow(hDenom->GetBinError(b)/hDenom->GetBinContent(b),2);
    if (FullCorr==0){
      if (hDenom->GetBinError(b)<hNum->GetBinError(b))      ErrC=pow(hDenom->GetBinError(b),2)/(hNum->GetBinContent(b)*hDenom->GetBinContent(b)); //Num is a subsample of denom
      else       ErrC=pow(hNum->GetBinError(b),2)/(hNum->GetBinContent(b)*hDenom->GetBinContent(b)); //Denom is a subsample of num
    }
    else if (FullCorr==1){
      ErrC=hDenom->GetBinError(b) * hNum->GetBinError(b)/(hNum->GetBinContent(b)*hDenom->GetBinContent(b));
    }
    Err=sqrt(Err1+Err2-2*ErrC); //are we sure there's the 2?
    if (Err1+Err2-ErrC<0) {
      cout << "Error not defined! (NaN)" << endl;
      Err = sqrt(Err1+Err2);
    }
    hRatio->SetBinError(b,Err*hRatio->GetBinContent(b));
    if (FullCorr==2){
      hRatio->SetBinError(b, sqrt(TMath::Abs(pow(hNum->GetBinError(b),2) - pow(hDenom->GetBinError(b),2)) ) /hDenom->GetBinContent(b));
    }
    //  cout << "bin: " << b << " err 1: " << Err*hRatio->GetBinContent(b) << " err 2: " <<  hRatio->GetBinError(b) << " 2/1: " << hRatio->GetBinError(b) / Err*hRatio->GetBinContent(b) << endl;
  }
}

//Nota:

//*************************************************************
//R = a/A where a is a subsample of A: ro = Sigma_A/Sigma_a <1
//(Sigma_R/R)^2 = (Sigma_a/a)^2 + (Sigma_A/A)^2 -2 (Sigma_A)^2 / (aA)

//R' = 1/R = A/a
//Sigma_R' = (Sigma_R/R)*R' = as above but a is now called Denom and A is now called Num
