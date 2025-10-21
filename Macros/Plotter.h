#include <TPaveText.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TLegend.h>
#include <TColor.h>

    void SetPaveText(TPaveText * pave, Font_t tfont,Float_t tsize, Style_t fstyle, Color_t fcolor,Short_t align,Int_t bordersize,Double_t x1ndc, Double_t x2ndc, Double_t y1ndc,Double_t y2ndc){
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
    TCanvas * CreateCanvas(TString padname="pad",Int_t width = 700,Int_t height=550,Bool_t twoPads = kFALSE, Bool_t threePads = kFALSE, Bool_t fourPads = kFALSE){
      TCanvas * can = new TCanvas(padname.Data(),padname.Data(),width,height);
      TPad * pad = 0x0;
      if(twoPads){
        pad = new TPad(padname.Data(),padname.Data(),0.001,0.3,0.999,0.999);
        pad ->SetMargin(0.12,0.02,0.01,0.05);
      }else if(threePads){
        pad = new TPad(padname.Data(),padname.Data(),0.001,0.3,0.75,0.999);
        pad ->SetMargin(0.12,0.01,0.01,0.01);
      }else if(fourPads){
        pad = new TPad(padname.Data(),padname.Data(),0.001,0.001,0.295,0.999);
        pad ->SetMargin(0.2,0.01,0.12,0.01);
      }else{
        pad = new TPad(padname.Data(),padname.Data(),0.001,0.001,0.999,0.999);
        pad ->SetMargin(0.12,0.02,0.12,0.05);
      }
      pad ->Draw();
      if(!twoPads&&!threePads&&!fourPads) pad->cd();
      pad->SetTicky();
      pad->SetTickx();
      can-> SetPadSave(pad);

      return can;
    }
    void SetGraphAxes(TGraph *gr, TString xTitle, TString yTitle) {
      gr->GetXaxis()->SetTitle(Form("%s",xTitle.Data()));
      gr->GetYaxis()->SetTitle(Form("%s",yTitle.Data()));
      gr->GetXaxis()->SetTitleSize(0.055);
      gr->GetYaxis()->SetTitleSize(0.055);
      gr->GetYaxis()->SetTitleOffset(1.);
      gr->GetXaxis()->SetLabelSize(0.045);
      gr->GetYaxis()->SetLabelSize(0.045);
    }
    void SetGraph(TGraph *gr, TString title,Int_t mStyle, Color_t col,Float_t mSize,Float_t alpha=0.5, Int_t line_s=1, Int_t line_w=1){
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
    void SetHistAxes(TH1 *hist, TString xTitle, TString yTitle){
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
    void SetHistAxesSmallPad(TH1 *hist, TString xTitle, TString yTitle){
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
    void SetGraphAxesSmallPad(TGraph *gr, TString xTitle, TString yTitle){
      gr->GetXaxis()->SetTitle(Form("%s",xTitle.Data()));
      gr->GetYaxis()->SetTitle(Form("%s",yTitle.Data()));
      gr->GetXaxis()->CenterTitle(kFALSE);
      gr->GetXaxis()->SetTitleSize(0.11);
      gr->GetYaxis()->SetNdivisions(505);
      gr->GetXaxis()->SetNdivisions(505);
      gr->GetYaxis()->SetTitleSize(0.11);
      gr->GetYaxis()->SetTitleOffset(0.5);
      gr->GetXaxis()->SetLabelSize(0.11);
      gr->GetYaxis()->SetLabelSize(0.11);
    }
    void SetHist(TH1 *hist, TString title,Int_t mStyle, Color_t col,Float_t mSize, Bool_t fill=kFALSE, Int_t line_s=1, Int_t line_w=1){
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
    TF1 * DrawConstant(Color_t col, Float_t par = 1,Float_t min=0.2, Float_t max=10, Int_t style = 9){
      TF1 * fConstant =  new TF1("fConstant","[0]",min,max);
      fConstant->SetParameter(0,par);
      fConstant->SetLineColor(col);
      fConstant->SetLineWidth(3);
      fConstant->SetLineStyle(style);
      fConstant->Draw("same");
      return fConstant;
    }
    TF1 * DrawLin(Color_t col, Float_t parConst, Float_t parLin, Float_t min, Float_t max, Int_t style = 9){
      TF1 * fLin =  new TF1("fLin","[0]+[1]*x",min,max);
      fLin->SetParameter(0,parConst);
      fLin->SetParameter(1,parLin);
      fLin->SetLineColor(col);
      fLin->SetLineWidth(2);
      fLin->SetLineStyle(style);
      fLin->Draw("same");
      return fLin;
    }
    TLegend * CreateLegend(Float_t x1pos,Float_t x2pos,Float_t y1pos, Float_t y2pos,Float_t textSize){
      TLegend *lg = new TLegend(x1pos,y1pos,x2pos,y2pos);
      lg->SetTextSize(textSize);
      lg->SetBorderSize(0);
      lg->SetFillStyle(0);

      return lg;
    }
    void Set2DHistAxes(TH2 *hist, TString xTitle, TString yTitle, TString zTitle, TString title){
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
    void SetFunctionAxes(TF1 *func, TString xTitle, TString yTitle){
      func->GetXaxis()->SetTitle(Form("%s",xTitle.Data()));
      func->GetYaxis()->SetTitle(Form("%s",yTitle.Data()));
      func->GetXaxis()->SetTitleSize(0.055);
      func->GetYaxis()->SetTitleSize(0.055);
      func->GetYaxis()->SetTitleOffset(1.);
      func->GetXaxis()->SetLabelSize(0.045);
      func->GetYaxis()->SetLabelSize(0.045);
      func->SetTitle("");
    }
    void SetColorPalete(const Int_t nClasses, Int_t colors[]) {
      const int NRGBs = 5;
      double stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
      double red[NRGBs]   = { 0.00, 0.00, 0.9*0.87, 1.00, 0.51 };
      double green[NRGBs] = { 0.00, 0.81, 0.9*1.00, 0.20, 0.00 };
      double blue[NRGBs]  = { 0.51, 0.9*1.00, 0.12, 0.00, 0.00 };

      int FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, nClasses+1);
      for(int i=0; i< nClasses ; i++)
      {
          colors[i] = FI + nClasses - i;

      }
    }
