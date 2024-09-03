#include "Riostream.h"
#include "TTimer.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TH3F.h>
#include "TNtuple.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TProfile.h"
#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TFile.h>
#include <TLine.h>
#include <TSpline.h>
#include "TFitResult.h"
#include "TGraphAsymmErrors.h"

const Int_t nPtTriggBins = 4;
const Int_t nMultBins = 9;

Float_t LimSupSpectra = 0.99;
Float_t LimInfSpectra = 0.4 * 1e-5; // 0.4 * 1e-6 for jet
Float_t xTitle = 15;
Float_t xOffset = 4;
Float_t yTitle = 30;
Float_t yOffset = 2;

Float_t xLabel = 30;
Float_t yLabel = 30;
Float_t xLabelOffset = 0.05;
Float_t yLabelOffset = 0.01;

Float_t tickX = 0.025;
Float_t tickY = 0.03;

Float_t LimSupMultRatio = 2.5; // 1.5 for jet
Float_t UpRangePt[4] = {15, 10, 10, 6};
Float_t LimInfMultRatio = 0; // 0.5 for jet
Float_t YoffsetSpectraRatio = 1.1;
Float_t xTitleR = 35;
Float_t xOffsetR = 1;
Float_t yTitleR = 30;
Float_t yOffsetR = 2;

Float_t xLabelR = 25;
Float_t yLabelR = 30;
Float_t xLabelOffsetR = 0.02;
Float_t yLabelOffsetR = 0.035;

Float_t tickXRatio = 0.04;
Float_t tickYRatio = 0.04;

Int_t ColorFit[] = {634, 797, 815, 429, 867, 601};
Int_t ColorModel[] = {kBlue + 1, kGreen + 2};
Int_t ColorMult[] = {634, 628, 807, kOrange - 4, 797, 815, 418, 429, 867, 856, 601, kViolet, kPink + 9, kPink + 1, 1};
// Float_t SizeMult[] = {2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8, 2, 2, 2.8, 2.5, 2.8};
//  Float_t SizeMultRatio[] = {1, 1, 1.8, 1.5, 1.8, 1, 1, 1.8, 1.5, 1.8, 1, 1, 1.8, 1.5, 1.8};
//  Int_t MarkerMult[] = {20, 21, 33, 34, 29, 20, 21, 33, 34, 29, 20, 21, 33, 34, 29};
Float_t SizeMult[] = {1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8, 1.8};
Float_t SizeMultRatio[] = {1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2};
Int_t MarkerMult[] = {21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21};
// Float_t ScaleFactor[] = {16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1};
Float_t ScaleFactor[] = {256, 128, 64, 32, 16, 8, 4, 2, 1};
Int_t ColorMB = 1;
Float_t SizeMB = 2;
Int_t MarkerMB = 24;
Float_t ScaleFactorMB = pow(2, 8);

TString particleName[] = {"K0s", "Lam", "Xi", "Omega", "Pion"};
Float_t particleMass[] = {0.497, 1.115, 1.321, 1.672, 0.1396};
TString finalNames[] = {"K_{S}^{0}", "(#Lambda+#bar{#Lambda})", "(#Xi^{+}+#Xi^{-})", "(#Omega^{+}+#Omega^{-})", "#pi^{+}+#pi^{-}"};
TString multiplicityPave[] = {"MB", "0-1%", "1-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-70%", "70-100%"};
TString multiplicityNames[] = {"minBias", "0_1Mult", "1_10Mult", "10_20Mult", "20_30Mult", "30_40Mult", "40_50Mult", "50_70Mult", "70_100Mult"};
TString multiplicityNamesShort[] = {"MB", "0_1Mult", "1_10Mult", "10_20Mult", "20_30Mult", "30_40Mult", "40_50Mult", "50_70Mult", "70_100Mult"};
TString namesRegions[3] = {"fHistNear", "fHistAway", "fHistUE"};
TString namesRegionsShort[3] = {"Near", "Away", "UE"};
TString paveRegions[3] = {"Near-side", "Away-side", "Underlying event"};
TString PhiRegions[3] = {"|#Delta#varphi| < #pi/2", "#pi/2 < |#Delta#varphi| < 3/2#pi", "-#pi/2 < |#Delta#varphi| < 3/2#pi"};
Color_t colRegions[3][3] = {{kRed + 1, kRed - 4, kRed + 2}, {kBlue, kAzure + 7, kBlue + 1}, {kGreen + 2, kGreen + 1, kGreen + 3}};
Int_t markers[4] = {20, 21, 29, 33};
Double_t ptTriggBins[] = {2., 4., 6., 10., 50};

Int_t Color[] = {634, 628, 797, 815, 418, 429, 867, 601, 1};
Int_t ColorPt[] = {634, 628, 807, kOrange - 4, 797, 815, 418, 429, 867, 856, 601, kViolet, kPink + 9, kPink + 1, 1};
Int_t ColorPart[] = {634, kBlue + 2};

TString TitleXPt = "#it{p}_{T} (GeV/#it{c})";
TString TitleYYield = "1/#it{N}_{trigg} d#it{N}/d#it{p}_{T} (GeV/#it{c})^{-1}";
TString TitleYYieldPtInt = "1/#it{N}_{trigg} d#it{N}/d#it{y}";
TString TitleYYieldPtIntToMB = "(1/#it{N}_{evt} d#it{N}/d#it{y}) / (1/#it{N}_{evt} d#it{N}/d#it{y})_{MB}";
TString TitleXMult = "#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5}";
TString TitleXMultToMB = "#LTd#it{N}_{ch}/d#it{#eta}#GT_{|#it{#eta}|<0.5} / #LTd#it{N}_{ch}/d#it{#eta}#GT^{MB}_{|#it{#eta}|<0.5}";

void StylePad(TPad *pad, Float_t LMargin, Float_t RMargin, Float_t TMargin, Float_t BMargin)
{
  pad->SetFillColor(0);
  pad->SetTickx(1);
  pad->SetTicky(1);
  pad->SetLeftMargin(LMargin);
  pad->SetRightMargin(RMargin);
  pad->SetTopMargin(TMargin);
  pad->SetBottomMargin(BMargin);
}

void StyleHisto(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString TitleX, TString TitleY, TString title)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(1.5);
  histo->GetXaxis()->SetTitle(TitleX);
  histo->GetXaxis()->SetTitleSize(0.04);
  histo->GetXaxis()->SetTitleOffset(1.2);
  histo->GetYaxis()->SetTitle(TitleY);
  histo->GetYaxis()->SetTitleSize(0.04);
  histo->GetYaxis()->SetTitleOffset(1.3);
  histo->SetTitle(title);
}

void StyleHistoYield(TH1F *histo, Float_t Low, Float_t Up, Int_t color, Int_t style, TString TitleX, TString TitleY, TString title, Float_t mSize, Float_t xOffset, Float_t yOffset)
{
  histo->GetYaxis()->SetRangeUser(Low, Up);
  histo->SetLineColor(color);
  histo->SetMarkerColor(color);
  histo->SetMarkerStyle(style);
  histo->SetMarkerSize(mSize);
  histo->GetXaxis()->SetTitle(TitleX);
  histo->GetXaxis()->SetTitleSize(0.05);
  histo->GetXaxis()->SetLabelSize(0.05);
  histo->GetXaxis()->SetTitleOffset(xOffset);
  histo->GetYaxis()->SetTitle(TitleY);
  histo->GetYaxis()->SetTitleSize(0.05);
  histo->GetYaxis()->SetTitleOffset(yOffset); // 1.2
  histo->GetYaxis()->SetLabelSize(0.05);
  histo->SetTitle(title);
}

void SetFont(TH1F *histo)
{
  histo->GetXaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetYaxis()->SetLabelFont(43);
}
void SetTickLength(TH1F *histo, Float_t TickLengthX, Float_t TickLengthY)
{
  histo->GetXaxis()->SetTickLength(TickLengthX);
  histo->GetYaxis()->SetTickLength(TickLengthY);
}

void SetHistoTextSize(TH1F *histo, Float_t XSize, Float_t XLabelSize, Float_t XOffset, Float_t XLabelOffset, Float_t YSize, Float_t YLabelSize, Float_t YOffset, Float_t YLabelOffset)
{
  histo->GetXaxis()->SetTitleSize(XSize);
  histo->GetXaxis()->SetLabelSize(XLabelSize);
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetXaxis()->SetLabelOffset(XLabelOffset);
  histo->GetYaxis()->SetTitleSize(YSize);
  histo->GetYaxis()->SetLabelSize(YLabelSize);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetYaxis()->SetLabelOffset(YLabelOffset);
}

void StyleCanvas(TCanvas *canvas, Float_t TopMargin, Float_t BottomMargin, Float_t LeftMargin, Float_t RightMargin)
{
  canvas->SetFillColor(0);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  gPad->SetTopMargin(TopMargin);
  gPad->SetLeftMargin(LeftMargin);
  gPad->SetBottomMargin(BottomMargin);
  gPad->SetRightMargin(RightMargin);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(0);
  gStyle->SetLegendFont(42);
}