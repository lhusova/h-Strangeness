#include "Plotter.h"
#include "TFile.h"
#include "Style.h"
#include "ErrRatioCorr.C"

void PlotRegionComparison(Int_t iPart = 0)
{

  if (iPart > 4)
  {
    cout << "Error!!! Assoc particle not defined !!!" << endl;
    return;
  }

  TString particleName[] = {"K0s", "Lam", "Xi", "Omega", "Pion"};
  Float_t particleMass[] = {0.497, 1.115, 1.321, 1.672, 0.1396};
  TString finalNames[] = {"K_{S}^{0}", "(#Lambda+#bar{#Lambda})", "(#Xi^{+}+#Xi^{-})", "(#Omega^{+}+#Omega^{-})", "#pi^{+}+#pi^{-}"};
  // TString multiplicityNames[]={"0_10Mult","10_20Mult","20_30Mult","30_40Mult","40_50Mult","50_60Mult","60_70Mult","70_80Mult","80_90Mult","90_100Mult"};
  // TString multiplicityPave[]={"0-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%","80-90%","90-100%"};
  TString multiplicityNames[] = {"minBias", "0_1Mult", "1_10Mult", "10_20Mult", "20_30Mult", "30_40Mult", "40_50Mult", "50_70Mult", "70_100Mult"};
  TString multiplicityPave[] = {"MB", "0-1%", "1-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-70%", "70-100%"};
  TString namesRegions[3] = {"fHistNear", "fHistAway", "fHistUE"};
  TString paveRegions[3] = {"Near-side", "Away-side", "Underlying event"};
  TString paveRegionsShort[3] = {"NS", "AS", "UE"};
  TString PhiRegions[3] = {"|#Delta#varphi| < #pi/2", "#pi/2 < |#Delta#varphi| < 3/2#pi", "-#pi/2 < |#Delta#varphi| < 3/2#pi"};
  Color_t colRegions[3][3] = {{kRed + 1, kRed - 4, kRed + 2}, {kBlue, kAzure + 7, kBlue + 1}, {kGreen + 2, kGreen + 1, kGreen + 3}};
  Int_t markers[4] = {20, 21, 29, 33};
  Double_t ptTriggBins[] = {2., 4., 6., 10., 100};
  // Double_t ptTriggBins[] = {0., 1., 2., 3., 100};

  TFile *file[10][4];
  TH1F *fProj[10][4][3];
  TCanvas *can[10][4][3];
  TCanvas *canMult[4][3];
  TPaveText *pave[10][4][3];
  TLegend *leg = Plotter::CreateLegend(0.65, 0.95, 0.25, 0.5, 0.05);
  Double_t yield, erryield, erryieldSist, err;

  for (Int_t iPtTrigg = 0; iPtTrigg < 3; iPtTrigg++)
  {
    for (Int_t iMult = 0; iMult < 9; iMult++)
    {
      for (Int_t iReg = 0; iReg < 3; iReg++)
      {
        file[iMult][iPtTrigg] = new TFile(Form("../../K0_Yields/Yields_longTrain_%s_%s_fullrangePeak_11_flat_ptTrigg%d.root", particleName[iPart].Data(), multiplicityNames[iMult].Data(), iPtTrigg));
        fProj[iMult][iPtTrigg][iReg] = (TH1F *)file[iMult][iPtTrigg]->Get(namesRegions[iReg].Data());
        if (iReg == 2)
          fProj[iMult][iPtTrigg][iReg]->Scale(50);
        // fProj[iMult][iPtTrigg][iReg]->GetYaxis()->SetMaxDigits(2);
        // fProj[iMult][iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0, fProj[iMult][iPtTrigg][iReg]->GetMaximum() * 1.2);
        // fProj[iMult][iPtTrigg][iReg]->GetXaxis()->SetRangeUser(0, 15);
      }
    }
  }

  // PLOT: YIELDS VS PT in MULT CLASSES
  TH1F *fProjScaled[10][4][3];
  TH1F *fProjSistScaled[10][4][3];
  TH1F *fHistSpectrumStatMultRatio[9];
  TH1F *fHistSpectrumSistMultRatio[9];
  TCanvas *canvasPtSpectra;
  Float_t LLUpperPad = 0.33;
  Float_t ULLowerPad = 0.33;
  TPad *pad1;
  TPad *padL1;
  TString sScaleFactorFinal[9];
  Float_t ScaleFactorFinal[9];

  for (Int_t iPtTrigg = 0; iPtTrigg < 3; iPtTrigg++)
  {
    canvasPtSpectra = new TCanvas(Form("canvasPtSpectra_pttrigg%i", iPtTrigg), Form("canvasPtSpectra_pttrigg%i", iPtTrigg), 700, 900);
    pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
    padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);

    StylePad(pad1, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
    StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B

    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendBorderSize(0);

    TLegend *legendAllMult = new TLegend(0.22, 0.03, 0.93, 0.28);
    legendAllMult->SetHeader("");
    legendAllMult->SetNColumns(2);
    legendAllMult->SetFillStyle(0);
    TLegendEntry *lheaderAllMult = (TLegendEntry *)legendAllMult->GetListOfPrimitives()->First();
    lheaderAllMult->SetTextSize(0.04);

    TLegend *LegendTitle = new TLegend(0.54, 0.65, 0.95, 0.92);
    LegendTitle->SetFillStyle(0);
    LegendTitle->SetTextAlign(33);
    LegendTitle->SetTextSize(0.04);
    LegendTitle->AddEntry("", "#bf{ALICE Work in progress}", "");
    LegendTitle->AddEntry("", "Run 3, pp #sqrt{#it{s}} = 13.6 TeV", "");
    LegendTitle->AddEntry("", Form("h-%s, %.1f < #it{p}_{T, trigg} < %.1f GeV/#it{c}", finalNames[iPart].Data(), ptTriggBins[iPtTrigg], ptTriggBins[iPtTrigg + 1]), "");
    // LegendTitle->AddEntry("", Form("%.1f < p_{T, trigg} < %.1f", ptTriggBins[iPtTrigg], ptTriggBins[iPtTrigg + 1]), "");
    LegendTitle->AddEntry("", Form("|#it{#eta}_{trigg}| < 0.8, |#it{#eta}_{%s}| < 0.8", finalNames[iPart].Data()), "");
    // LegendTitle->AddEntry("", "|#Delta#it{#eta}| < 1.1", "");

    LimSupSpectra = 9999.99;
    LimInfSpectra = 0.4 * 1e-6;
    TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 15.5);
    for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
      hDummy->SetBinContent(i, 1e-12);
    canvasPtSpectra->cd();
    SetFont(hDummy);
    StyleHistoYield(hDummy, LimInfSpectra, LimSupSpectra, 1, 1, TitleXPt, TitleYYield, "", 1, 1.15, 1.6);
    SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
    SetTickLength(hDummy, tickX, tickY);
    hDummy->GetXaxis()->SetRangeUser(0, 8.5);
    pad1->Draw();
    pad1->cd();
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    hDummy->Draw("same");

    for (Int_t iReg = 0; iReg < 3; iReg++)
    {
      // LegendTitle->AddEntry("", Form("%s", paveRegions[iReg].Data()), "");

      if (iReg == 2)
        ScaleFactorMB = pow(2, 10);
      else
        ScaleFactorMB = pow(2, 8);
      Int_t iMultEff = -1;
      for (Int_t iMult = 0; iMult < 9; iMult++)
      {
        if (iMult != 0 && iMult != 1 && iMult != 8)
          continue; // keep only 0-1%, 70-100% and 0-100%
        iMultEff++;
        ScaleFactorFinal[iMult] = ScaleFactor[iMult];
        if (iMult == 0)
          ScaleFactorFinal[iMult] = ScaleFactorMB;
        fProjScaled[iMult][iPtTrigg][iReg] = (TH1F *)fProj[iMult][iPtTrigg][iReg]->Clone(Form("fHistSpectrumStatScaled_%s_pttrigg%i", multiplicityNames[iMult].Data(), iPtTrigg));
        fProjSistScaled[iMult][iPtTrigg][iReg] = (TH1F *)fProj[iMult][iPtTrigg][iReg]->Clone(Form("fHistSpectrumSistScaled_%s_pttrigg%i", multiplicityNames[iMult].Data(), iPtTrigg));
        fProjScaled[iMult][iPtTrigg][iReg]->Scale(ScaleFactorFinal[iMult]);
        fProjSistScaled[iMult][iPtTrigg][iReg]->Scale(ScaleFactorFinal[iMult]);
        for (Int_t b = 1; b <= fProj[iMult][iPtTrigg][iReg]->GetNbinsX(); b++)
        {
          // cout << "bin " << b << " " << fProj[iMult][iPtTrigg][iReg]->GetBinContent(b) << "+-" << fProj[iMult][iPtTrigg][iReg]->GetBinError(b) << endl;
          //  cout << "bin " << b << " " << fProjSist[iMult][iPtTrigg][iReg]->GetBinContent(b) << "+-" << fProjSist[iMult][iPtTrigg][iReg]->GetBinError(b) << endl;
          // cout << "bin " << b << " " << fProjScaled[iMult][iPtTrigg][iReg]->GetBinContent(b) << "+-" << fProjScaled[iMult][iPtTrigg][iReg]->GetBinError(b) << endl;
        }
        fProjScaled[iMult][iPtTrigg][iReg]->SetMarkerColor(colRegions[iReg][iMultEff]);
        fProjScaled[iMult][iPtTrigg][iReg]->SetLineColor(colRegions[iReg][iMultEff]);
        fProjScaled[iMult][iPtTrigg][iReg]->SetMarkerStyle(MarkerMult[iMult]);
        fProjScaled[iMult][iPtTrigg][iReg]->SetMarkerSize(SizeMult[iMult]);
        fProjSistScaled[iMult][iPtTrigg][iReg]->SetMarkerColor(colRegions[iReg][iMultEff]);
        fProjSistScaled[iMult][iPtTrigg][iReg]->SetLineColor(colRegions[iReg][iMultEff]);
        fProjSistScaled[iMult][iPtTrigg][iReg]->SetMarkerStyle(MarkerMult[iMult]);
        fProjSistScaled[iMult][iPtTrigg][iReg]->SetMarkerSize(SizeMult[iMult]);
        if (iMult != 0)
        {
          fProjScaled[iMult][iPtTrigg][iReg]->Draw("same ex0");
          fProjSistScaled[iMult][iPtTrigg][iReg]->SetFillStyle(0);
          fProjSistScaled[iMult][iPtTrigg][iReg]->Draw("same e2");
        }
        sScaleFactorFinal[iMult] = Form(" (x2^{%i})", int(log2(ScaleFactorFinal[iMult])));
        if (ScaleFactorFinal[iMult] == 1)
          sScaleFactorFinal[iMult] = "";
        else if (ScaleFactorFinal[iMult] == 2)
          sScaleFactorFinal[iMult] = " (x2)";
        if (iMult!=0) legendAllMult->AddEntry(fProjScaled[iMult][iPtTrigg][iReg], Form("%s %s", paveRegionsShort[iReg].Data(), multiplicityPave[iMult].Data()) + sScaleFactorFinal[iMult] + " ", "pef");
      } // end loop on mult
    }
    LegendTitle->Draw("");
    legendAllMult->Draw("");

    Int_t ChosenMult = 0; // MB
    TString TitleYSpectraRatio = Form("Ratio to %s", multiplicityPave[ChosenMult].Data());
    TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 10000, 0, 15.5);
    SetFont(hDummyRatio);
    StyleHistoYield(hDummyRatio, LimInfMultRatio, LimSupMultRatio, 1, 1, TitleXPt, TitleYSpectraRatio, "", 1, 1.15, YoffsetSpectraRatio);
    SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
    SetTickLength(hDummyRatio, tickXRatio, tickYRatio);
    hDummyRatio->GetXaxis()->SetRangeUser(0, 15.5);
    canvasPtSpectra->cd();
    padL1->Draw();
    padL1->cd();
    // gPad->SetLogy();
    hDummyRatio->Draw("same");

    for (Int_t iReg = 0; iReg < 3; iReg++)
    {
      Int_t iMultEff = -1;
      for (Int_t iMult = 0; iMult < 9; iMult++)
      {
        if (iMult != 0 && iMult != 1 && iMult != 8)
          continue; // keep only 0-1%, 70-100% and 0-100%
        iMultEff++;  
        fHistSpectrumStatMultRatio[iMult] = (TH1F *)fProj[iMult][iPtTrigg][iReg]->Clone(Form("fHistSpectrumStatRatio_%s_pttrigg%i_%s", multiplicityNames[iMult].Data(), iPtTrigg, namesRegions[iReg].Data()));
        fHistSpectrumSistMultRatio[iMult] = (TH1F *)fProj[iMult][iPtTrigg][iReg]->Clone(Form("fHistSpectrumSistRatio_%s_pttrigg%i_%s", multiplicityNames[iMult].Data(), iPtTrigg, namesRegions[iReg].Data()));
        fHistSpectrumStatMultRatio[iMult]->Divide(fProj[ChosenMult][iPtTrigg][iReg]);
        fHistSpectrumSistMultRatio[iMult]->Divide(fProj[ChosenMult][iPtTrigg][iReg]);
        ErrRatioCorr(fProj[iMult][iPtTrigg][iReg], fProj[ChosenMult][iPtTrigg][iReg], fHistSpectrumStatMultRatio[iMult], 0);
        ErrRatioCorr(fProj[iMult][iPtTrigg][iReg], fProj[ChosenMult][iPtTrigg][iReg], fHistSpectrumSistMultRatio[iMult], 0);
        fHistSpectrumStatMultRatio[iMult]->GetYaxis()->SetRangeUser(LimInfMultRatio, LimSupMultRatio);
        fHistSpectrumSistMultRatio[iMult]->GetYaxis()->SetRangeUser(LimInfMultRatio, LimSupMultRatio);
        fHistSpectrumStatMultRatio[iMult]->SetMarkerColor(colRegions[iReg][iMultEff]);
        fHistSpectrumStatMultRatio[iMult]->SetLineColor(colRegions[iReg][iMultEff]);
        fHistSpectrumStatMultRatio[iMult]->SetMarkerStyle(MarkerMult[iMult]);
        fHistSpectrumStatMultRatio[iMult]->SetMarkerSize(SizeMultRatio[iMult]);
        fHistSpectrumSistMultRatio[iMult]->SetMarkerColor(colRegions[iReg][iMultEff]);
        fHistSpectrumSistMultRatio[iMult]->SetLineColor(colRegions[iReg][iMultEff]);
        fHistSpectrumSistMultRatio[iMult]->SetMarkerStyle(MarkerMult[iMult]);
        fHistSpectrumSistMultRatio[iMult]->SetMarkerSize(SizeMultRatio[iMult]);
        if (iMult != ChosenMult)
        {
          fHistSpectrumStatMultRatio[iMult]->Draw("same ex0");
          fHistSpectrumSistMultRatio[iMult]->SetFillStyle(0);
          fHistSpectrumSistMultRatio[iMult]->Draw("same e2");
        }
      } // end loop on mult
    }
    TString stringoutpdf = Form("../../Ptspectrum_%s_AllRegions_ptTrigg%d", particleName[iPart].Data(), iPtTrigg);
    canvasPtSpectra->SaveAs(stringoutpdf + ".pdf");
    canvasPtSpectra->SaveAs(stringoutpdf + ".png");
    canvasPtSpectra->SaveAs(stringoutpdf + ".eps");
  } // end loop on pttrigg
}
