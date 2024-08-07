#include "Plotter.h"
#include <TPad.h>

#if !defined(R__ALPHA) && !defined(R__SOLARIS) && !defined(R__ACC) && !defined(R__FBSD)
NamespaceImp(Plotter)
#endif
namespace Plotter{}
//_________________________________
void Plotter::SetPaveText(TPaveText * pave, Font_t tfont,Float_t tsize, Style_t fstyle, Color_t fcolor,Short_t align,Int_t bordersize,Double_t x1ndc, Double_t x2ndc, Double_t y1ndc,Double_t y2ndc){
    pave->SetTextFont(tfont);
    pave->SetTextSize(tsize);
    pave->SetFillStyle(fstyle);
    pave->SetFillColor(fcolor);
    pave->SetTextAlign(align);
    pave->SetBorderSize(bordersize);
    pave->SetX1NDC(x1ndc);
    pave->SetX2NDC(x2ndc);
    pave->SetY1NDC(y1ndc);
    pave->SetY2NDC(y2ndc);
}
//_________________________________
TCanvas * Plotter::CreateCanvas(TString padname,Int_t width,Int_t height,Bool_t twoPads){
  TCanvas * can = new TCanvas(padname.Data(),padname.Data(),width,height);
  TPad * pad = 0x0;
  if(twoPads){
    pad = new TPad(padname.Data(),padname.Data(),0.001,0.3,0.999,0.999);
    pad ->SetMargin(0.12,0.02,0.01,0.05);
  }else{
    pad = new TPad(padname.Data(),padname.Data(),0.001,0.001,0.999,0.999);
    pad ->SetMargin(0.12,0.02,0.12,0.05);
  }
  pad ->Draw();
  if(!twoPads) pad->cd();
  pad->SetTicky();
  pad->SetTickx();
  can-> SetPadSave(pad);

  return can;
}
//_________________________________
void Plotter::SetGraphAxes(TGraph *gr, TString xTitle, TString yTitle) {
  gr->GetXaxis()->SetTitle(Form("%s",xTitle.Data()));
  gr->GetYaxis()->SetTitle(Form("%s",yTitle.Data()));
  gr->GetXaxis()->SetTitleSize(0.055);
  gr->GetYaxis()->SetTitleSize(0.055);
  gr->GetYaxis()->SetTitleOffset(1.);
  gr->GetXaxis()->SetLabelSize(0.045);
  gr->GetYaxis()->SetLabelSize(0.045);
}
//_________________________________
void Plotter::SetHistAxes(TH1 *hist, TString xTitle, TString yTitle) {
  hist->GetXaxis()->SetTitle(Form("%s",xTitle.Data()));
  hist->GetYaxis()->SetTitle(Form("%s",yTitle.Data()));
  hist->GetXaxis()->CenterTitle(kFALSE);
  hist->GetXaxis()->SetTitleSize(0.055);
  hist->GetYaxis()->SetTitleSize(0.055);
  hist->GetYaxis()->SetTitleOffset(1.);
  hist->GetXaxis()->SetLabelSize(0.045);
  hist->GetYaxis()->SetLabelSize(0.045);
  hist->SetStats(kFALSE);
}
//_________________________________
void Plotter::SetFunctionAxes(TF1 *func, TString xTitle, TString yTitle) {
  func->GetXaxis()->SetTitle(Form("%s",xTitle.Data()));
  func->GetYaxis()->SetTitle(Form("%s",yTitle.Data()));
  func->GetXaxis()->SetTitleSize(0.055);
  func->GetYaxis()->SetTitleSize(0.055);
  func->GetYaxis()->SetTitleOffset(1.);
  func->GetXaxis()->SetLabelSize(0.045);
  func->GetYaxis()->SetLabelSize(0.045);
  func->SetTitle("");
}
//_________________________________
void Plotter::SetHistAxesSmallPad(TH1 *hist, TString xTitle, TString yTitle) {
  hist->GetXaxis()->SetTitle(Form("%s",xTitle.Data()));
  hist->GetYaxis()->SetTitle(Form("%s",yTitle.Data()));
  hist->GetXaxis()->CenterTitle(kFALSE);
  hist->GetXaxis()->SetTitleSize(0.11);
  hist->GetYaxis()->SetNdivisions(505);
  hist->GetYaxis()->SetTitleSize(0.11);
  hist->GetYaxis()->SetTitleOffset(0.5);
  hist->GetXaxis()->SetLabelSize(0.11);
  hist->GetYaxis()->SetLabelSize(0.11);
  hist->SetStats(kFALSE);
}
//_________________________________
void Plotter::SetGraph(TGraph *gr, TString title,Int_t mStyle, Color_t col,Float_t mSize,Float_t alpha, Int_t line_s, Int_t line_w){

  gr->SetTitle(title.Data());
  gr->SetMarkerStyle(mStyle);
  gr->SetMarkerColor(col);
  gr->SetLineColor(col);
  gr->SetMarkerSize(mSize);
  gr->SetFillStyle(1001);
  gr->SetFillColorAlpha(col, alpha);
  gr->SetLineStyle(line_s);
  gr->SetLineWidth(line_w);

}
//_________________________________
void Plotter::SetHist(TH1 *hist, TString title,Int_t mStyle, Color_t col,Float_t mSize, Bool_t fill, Int_t line_s, Int_t line_w){

  hist->SetTitle(title.Data());
  hist->SetMarkerStyle(mStyle);
  hist->SetMarkerColor(col);
  hist->SetLineColor(col);
  hist->SetMarkerSize(mSize);
  if(fill){
    hist->SetFillStyle(1001);
    hist->SetFillColorAlpha(col, 0.5);
  }
  hist->SetLineStyle(line_s);
  hist->SetLineWidth(line_w);

}
//_________________________________
TF1 * Plotter::DrawUnity(Color_t col, Float_t par, Float_t min, Float_t max, Int_t style){
  TF1 * fUnity =  new TF1("fUnity","[0]",min,max);
  fUnity->SetParameter(0,par);
  fUnity->SetLineColor(col);
  fUnity->SetLineWidth(2);
  fUnity->SetLineStyle(style);
  fUnity->Draw("same");
  return fUnity;
}
//_________________________________
TF1 * Plotter::DrawLin(Color_t col, Float_t parConst, Float_t parLin, Float_t min, Float_t max, Int_t style){
  TF1 * fLin =  new TF1("fUnity","[0]+[1]*x",min,max);
  fLin->SetParameter(0,parConst);
  fLin->SetParameter(1,parLin);
  fLin->SetLineColor(col);
  fLin->SetLineWidth(2);
  fLin->SetLineStyle(style);
  fLin->Draw("same");
  return fLin;
}
//_________________________________
TLegend * Plotter::CreateLegend(Float_t x1pos,Float_t x2pos,Float_t y1pos, Float_t y2pos,Float_t textSize){
  TLegend *lg = new TLegend(x1pos,y1pos,x2pos,y2pos);
  lg->SetTextSize(textSize);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);

  return lg;
}
//_________________________________
void Plotter::Set2DHistAxes(TH2 *hist, TString xTitle, TString yTitle, TString zTitle, TString title){
  hist->GetXaxis()->SetTitle(Form("%s",xTitle.Data()));
  hist->GetXaxis()->CenterTitle();
  hist->GetYaxis()->SetTitle(Form("%s",yTitle.Data()));
  hist->GetYaxis()->CenterTitle();
  hist->GetZaxis()->SetTitle(Form("%s",zTitle.Data()));
  hist->SetTitle(Form("%s",title.Data()));
  hist->GetXaxis()->SetTitleSize(0.055);
  hist->GetYaxis()->SetTitleSize(0.055);
  hist->GetZaxis()->SetTitleSize(0.055);
  hist->GetXaxis()->SetTitleOffset(1.);
  hist->GetYaxis()->SetTitleOffset(1.);
  hist->GetZaxis()->SetTitleOffset(1.5);
  hist->GetXaxis()->SetLabelSize(0.045);
  hist->GetYaxis()->SetLabelSize(0.045);
  hist->GetZaxis()->SetLabelSize(0.045);
  hist->SetStats(kFALSE);
}
