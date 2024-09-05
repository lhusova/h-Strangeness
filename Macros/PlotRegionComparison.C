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

  TFile *file[nMultBins][nPtTriggBins];
  TFile *fileSyst;
  TH1F *fProj[nMultBins][nPtTriggBins][nRegions];
  TH1F *fProjRelSyst[nMultBins][nPtTriggBins][nRegions];
  TH1F *fProjSyst[nMultBins][nPtTriggBins][nRegions];

  fileSyst = new TFile("../../Uncertainty.root", "");
  if (!fileSyst)
  {
    cout << "File not found" << endl;
    return;
  }

  for (Int_t iPtTrigg = 0; iPtTrigg < nPtTriggBins; iPtTrigg++)
  {
    for (Int_t iMult = 0; iMult < nMultBins; iMult++)
    {
      for (Int_t iReg = 0; iReg < nRegions; iReg++)
      {
        file[iMult][iPtTrigg] = new TFile(Form("../../K0_Yields_Kai/Yields_%s_%s_fullrangePeak_11_flat_ptTrigg%d.root", particleName[iPart].Data(), multiplicityNamesShort[iMult].Data(), iPtTrigg));
        if (!file[iMult][iPtTrigg])
        {
          cout << "File not found" << endl;
          return;
        }
        fProj[iMult][iPtTrigg][iReg] = (TH1F *)file[iMult][iPtTrigg]->Get(namesRegions[iReg].Data());
        if (!fProj[iMult][iPtTrigg][iReg])
        {
          cout << "Histogram not found" << endl;
          return;
        }

        fProj[iMult][iPtTrigg][iReg]->SetName(Form("fhist_%s_%s_%d", namesRegionsShort[iReg].Data(), multiplicityNames[iMult].Data(), iPtTrigg));
        fProjSyst[iMult][iPtTrigg][iReg] = (TH1F *)fProj[iMult][iPtTrigg][iReg]->Clone(Form("fhistSyst_%s_%s_%d", namesRegionsShort[iReg].Data(), multiplicityNames[iMult].Data(), iPtTrigg));
        fProjRelSyst[iMult][iPtTrigg][iReg] = (TH1F *)fileSyst->Get(Form("fhistsyst_%s_%s_pttrigger%d", namesRegionsShort[iReg].Data(), multiplicityNamesShort[0].Data(), iPtTrigg));
        if (!fProjRelSyst[iMult][iPtTrigg][iReg])
        {
          cout << "Histogram uncertainty not found" << endl;
          return;
        }
        if (iReg == 2)
        {
          fProj[iMult][iPtTrigg][iReg]->Scale(50);
          fProjSyst[iMult][iPtTrigg][iReg]->Scale(50);
        }
        for (Int_t b = 1; b <= fProj[iMult][iPtTrigg][iReg]->GetNbinsX(); b++)
        {
          fProjSyst[iMult][iPtTrigg][iReg]->SetBinError(b, fProjRelSyst[iMult][iPtTrigg][iReg]->GetBinContent(b) * fProjSyst[iMult][iPtTrigg][iReg]->GetBinContent(b));
        }
      }
    }
  }

  // PLOT: YIELDS VS PT in MULT CLASSES
  TH1F *fProjScaled[nMultBins][nPtTriggBins][nRegions];
  TH1F *fProjSistScaled[nMultBins][nPtTriggBins][nRegions];
  TH1F *fHistSpectrumStatMultRatio[nMultBins];
  TH1F *fHistSpectrumSistMultRatio[nMultBins];
  TCanvas *canvasPtSpectra;
  Float_t LLUpperPad = 0.33;
  Float_t ULLowerPad = 0.33;
  TPad *pad1;
  TPad *padL1;
  TString sScaleFactorFinal[nMultBins];
  Float_t ScaleFactorFinal[nMultBins];

  for (Int_t iPtTrigg = 0; iPtTrigg < nPtTriggBins; iPtTrigg++)
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
    LimSupSpectra = 99.99;
    if (iPtTrigg == 0)
    {
      LimInfSpectra = 0.4 * 1e-8;
      LimSupSpectra = 999.99;
    }
    else if (iPtTrigg == 1)
    {
      LimInfSpectra = 0.4 * 1e-7;
      LimSupSpectra = 9999.99;
    }
    else if (iPtTrigg == 2)
    {
      LimInfSpectra = 0.4 * 1e-7;
      LimSupSpectra = 9999.99;
    }
    else
    {
      LimInfSpectra = 0.4 * 1e-6;
      LimSupSpectra = 999999.99;
    }
    TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, 15.5);
    for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
      hDummy->SetBinContent(i, 1e-12);
    canvasPtSpectra->cd();
    SetFont(hDummy);
    StyleHistoYield(hDummy, LimInfSpectra, LimSupSpectra, 1, 1, TitleXPt, TitleYYield, "", 1, 1.15, 1.6);
    SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
    SetTickLength(hDummy, tickX, tickY);
    hDummy->GetXaxis()->SetRangeUser(0, UpRangePt[iPtTrigg]);
    pad1->Draw();
    pad1->cd();
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    hDummy->Draw("same");

    for (Int_t iReg = 0; iReg < nRegions; iReg++)
    {
      if (iReg == 2)
        ScaleFactorMB = pow(2, 10);
      else
        ScaleFactorMB = pow(2, 8);
      Int_t iMultEff = -1;
      for (Int_t iMult = 0; iMult < nMultBins; iMult++)
      {
        if (iMult != 0 && iMult != 1 && iMult != 8)
          continue; // keep only 0-1%, 70-100% and 0-100%
        iMultEff++;
        ScaleFactorFinal[iMult] = ScaleFactor[iMult];
        if (iMult == 0)
          ScaleFactorFinal[iMult] = ScaleFactorMB;
        fProjScaled[iMult][iPtTrigg][iReg] = (TH1F *)fProj[iMult][iPtTrigg][iReg]->Clone(Form("fHistSpectrumStatScaled_%s_pttrigg%i", multiplicityNames[iMult].Data(), iPtTrigg));
        fProjSistScaled[iMult][iPtTrigg][iReg] = (TH1F *)fProjSyst[iMult][iPtTrigg][iReg]->Clone(Form("fHistSpectrumSistScaled_%s_pttrigg%i", multiplicityNames[iMult].Data(), iPtTrigg));
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
        if (iMult != 0)
          legendAllMult->AddEntry(fProjScaled[iMult][iPtTrigg][iReg], Form("%s %s", namesRegionsShort[iReg].Data(), multiplicityPave[iMult].Data()) + sScaleFactorFinal[iMult] + " ", "pef");
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
    hDummyRatio->GetXaxis()->SetRangeUser(0, UpRangePt[iPtTrigg]);
    canvasPtSpectra->cd();
    padL1->Draw();
    padL1->cd();
    // gPad->SetLogy();
    hDummyRatio->Draw("same");

    for (Int_t iReg = 0; iReg < nRegions; iReg++)
    {
      Int_t iMultEff = -1;
      for (Int_t iMult = 0; iMult < nMultBins; iMult++)
      {
        if (iMult != 0 && iMult != 1 && iMult != 8)
          continue; // keep only 0-1%, 70-100% and 0-100%
        iMultEff++;
        fHistSpectrumStatMultRatio[iMult] = (TH1F *)fProj[iMult][iPtTrigg][iReg]->Clone(Form("fHistSpectrumStatRatio_%s_pttrigg%i_%s", multiplicityNames[iMult].Data(), iPtTrigg, namesRegions[iReg].Data()));
        fHistSpectrumSistMultRatio[iMult] = (TH1F *)fProjSyst[iMult][iPtTrigg][iReg]->Clone(Form("fHistSpectrumSistRatio_%s_pttrigg%i_%s", multiplicityNames[iMult].Data(), iPtTrigg, namesRegions[iReg].Data()));
        fHistSpectrumStatMultRatio[iMult]->Divide(fProj[ChosenMult][iPtTrigg][iReg]);
        fHistSpectrumSistMultRatio[iMult]->Divide(fProjSyst[ChosenMult][iPtTrigg][iReg]);
        ErrRatioCorr(fProj[iMult][iPtTrigg][iReg], fProj[ChosenMult][iPtTrigg][iReg], fHistSpectrumStatMultRatio[iMult], 0);
        ErrRatioCorr(fProjSyst[iMult][iPtTrigg][iReg], fProjSyst[ChosenMult][iPtTrigg][iReg], fHistSpectrumSistMultRatio[iMult], 0);
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
