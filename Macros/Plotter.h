#include <TPaveText.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TLegend.h>

namespace Plotter {
    void SetPaveText(TPaveText * pave, Font_t tfont,Float_t tsize, Style_t fstyle, Color_t fcolor,Short_t align,Int_t bordersize,Double_t x1ndc, Double_t x2ndc, Double_t y1ndc,Double_t y2ndc);
    TCanvas * CreateCanvas(TString padname="pad",Int_t width = 700,Int_t height=550,Bool_t twoPads = kFALSE);
    void SetGraphAxes(TGraph *gr, TString xTitle, TString yTitle) ;
    void SetGraph(TGraph *gr, TString title,Int_t mStyle, Color_t col,Float_t mSize,Float_t alpha=0.5, Int_t line_s=1, Int_t line_w=1);
    void SetHistAxes(TH1 *hist, TString xTitle, TString yTitle) ;
    void SetHistAxesSmallPad(TH1 *hist, TString xTitle, TString yTitle);
    void SetHist(TH1 *hist, TString title,Int_t mStyle, Color_t col,Float_t mSize, Bool_t fill=kFALSE, Int_t line_s=1, Int_t line_w=1);
    TF1 * DrawUnity(Color_t col, Float_t par = 1,Float_t min=0.2, Float_t max=10, Int_t style = 9);
    TF1 * DrawLin(Color_t col, Float_t parConst, Float_t parLin, Float_t min, Float_t max, Int_t style = 9);
    TLegend * CreateLegend(Float_t x1pos,Float_t x2pos,Float_t y1pos, Float_t y2pos,Float_t textSize);
    void Set2DHistAxes(TH2 *hist, TString xTitle, TString yTitle, TString zTitle, TString title);
    void SetFunctionAxes(TF1 *func, TString xTitle, TString yTitle);
};
