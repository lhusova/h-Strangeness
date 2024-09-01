#include "Plotter.h"
#include "TFile.h"
#include "Style.h"
#include "ErrRatioCorr.C"

void IntegratedYieldVsMult(Int_t iPart = 0, Int_t iReg = 0)
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
  TString multiplicityNamesShort[] = {"MB", "0_1Mult", "1_10Mult", "10_20Mult", "20_30Mult", "30_40Mult", "40_50Mult", "50_70Mult", "70_100Mult"};
  TString multiplicityPave[] = {"MB", "0-1%", "1-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-70%", "70-100%"};
  TString namesRegions[3] = {"fHistNear", "fHistAway", "fHistUE"};
  TString namesRegionsShort[3] = {"Near", "Away", "UE"};
  TString paveRegions[3] = {"Near-side", "Away-side", "Underlying event"};
  TString PhiRegions[3] = {"|#Delta#varphi| < #pi/2", "#pi/2 < |#Delta#varphi| < 3/2#pi", "-#pi/2 < |#Delta#varphi| < 3/2#pi"};
  Color_t colRegions[3][5] = {{kRed - 7, kRed - 4, kRed + 1, kRed + 2, kRed + 4}, {kAzure + 6, kAzure + 7, kBlue, kBlue + 1, kBlue + 3}, {kSpring + 7, kGreen + 1, kGreen + 2, kGreen + 3, kGreen + 5}};
  Int_t markers[4] = {20, 21, 29, 33};
  Double_t ptTriggBins[] = {2., 4., 6., 10., 50};
  // Double_t ptTriggBins[] = {0., 1., 2., 3., 100};

  TFile *file[10][4];
  TFile *fileSyst;
  TH1F *fProj[10][4][3];
  TH1F *fProjRelSyst[10][4][3];
  TH1F *fProjSyst[10][4][3];
  TCanvas *can[10][4][3];
  TCanvas *canMult[4][3];
  TPaveText *pave[10][4][3];
  TLegend *leg = Plotter::CreateLegend(0.65, 0.95, 0.25, 0.5, 0.05);
  Double_t yield, erryield, erryieldSist, err;

  // const char *labels[11] = {"90-100%","80-90%","70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","0-10%","MB"};
  const char *labels[9] = {"70-100%", "50-70%", "40-50%", "30-40%", "20-30%", "10-20%", "1-10%", "0-1%", "MB"};
  TH1F *histYield[4][3];
  TH1F *histYieldToMB[4][3];
  TH1F *histYieldSist[4][3];
  TH1F *histYieldSistToMB[4][3];
  for (Int_t iPtTrigg = 0; iPtTrigg < 4; iPtTrigg++)
  {
    histYield[iPtTrigg][0] = new TH1F(Form("histYieldNear%d", iPtTrigg), "", 9, 0, 9);
    histYield[iPtTrigg][1] = new TH1F(Form("histYieldAway%d", iPtTrigg), "", 9, 0, 9);
    histYield[iPtTrigg][2] = new TH1F(Form("histYieldUE%d", iPtTrigg), "", 9, 0, 9);
    histYieldSist[iPtTrigg][0] = new TH1F(Form("histYieldSistNear%d", iPtTrigg), "", 9, 0, 9);
    histYieldSist[iPtTrigg][1] = new TH1F(Form("histYieldSistAway%d", iPtTrigg), "", 9, 0, 9);
    histYieldSist[iPtTrigg][2] = new TH1F(Form("histYieldSistUE%d", iPtTrigg), "", 9, 0, 9);
  }

  TF1 *boltzmann = new TF1("boltzmann", "x*[0]*TMath::Sqrt([2]*[2]+x*x)*TMath::Exp(-TMath::Sqrt([2]*[2]+x*x)/[1])", 0, 6);
  if (iPart == 0)
  {
    boltzmann->SetParameter(0, 0.04);
    boltzmann->SetParameter(1, 0.1);
  }
  if (iPart == 4)
  {
    boltzmann->SetParameter(0, 0.1);
    boltzmann->SetParameter(1, 1.5);
  }
  boltzmann->FixParameter(2, particleMass[iPart]);

  TF1 *fermiDir = new TF1("fermiDir", "x*[0]/(TMath::Exp(TMath::Sqrt([2]*[2]+x*x)/[1])+1)", 0, 6);
  fermiDir->SetParameter(0, 1);
  fermiDir->SetParameter(1, 160);
  fermiDir->FixParameter(2, particleMass[iPart]);

  TF1 *levy = new TF1("levy", "[3]*x/([0]*[1])*(([0]-1)*([0]-2))/([0]*[1]+[2]*([0]-2))*TMath::Power(1+(TMath::Sqrt([2]*[2]+x*x)-[2])/([0]*[1]),-[0])", 0, 6);
  levy->SetParameter(0, 7);
  levy->SetParameter(1, 0.8);
  levy->FixParameter(2, particleMass[iPart]);
  levy->SetParameter(3, 0.03);

  TF1 *mT = new TF1("mT", "x*[0]*TMath::Exp(-TMath::Sqrt([2]*[2]+x*x)/[1])", 0, 6);
  mT->SetParameter(0, 0.01);
  mT->SetParameter(1, 0.7);
  mT->SetParLimits(1, 0.1, 10);
  mT->FixParameter(2, particleMass[iPart]);

  for (Int_t iMult = 0; iMult < 9; iMult++)
  {
    fileSyst = new TFile("../../Uncertainty.root", "");
    if (!fileSyst)
    {
      cout << "File not found" << endl;
      return;
    }
    for (Int_t iPtTrigg = 0; iPtTrigg < 4; iPtTrigg++)
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
      can[iMult][iPtTrigg][iReg] = Plotter::CreateCanvas(Form("can%i%i%i", iMult, iPtTrigg, iReg));
      if (iReg == 2)
      {
        fProj[iMult][iPtTrigg][iReg]->Scale(50);
        fProjSyst[iMult][iPtTrigg][iReg]->Scale(50);
      }
      cout << fProj[iMult][iPtTrigg][iReg]->GetNbinsX() << " " << fProjRelSyst[iMult][iPtTrigg][iReg]->GetNbinsX() << endl;
      for (Int_t b = 1; b <= fProj[iMult][iPtTrigg][iReg]->GetNbinsX(); b++)
      {
        fProjSyst[iMult][iPtTrigg][iReg]->SetBinError(b, fProjRelSyst[iMult][iPtTrigg][iReg]->GetBinContent(b) * fProjSyst[iMult][iPtTrigg][iReg]->GetBinContent(b));
      }
      // fProj[iMult][iPtTrigg][iReg]->Fit(boltzmann,"mer");
      // fProj[iMult][iPtTrigg][iReg]->Fit(levy,"mer");
      // fProj[iMult][iPtTrigg][iReg]->Fit(fermiDir,"mer");
      // fProj[iMult][iPtTrigg][iReg]->Fit(mT,"mer");
      levy->SetLineColor(kBlue);
      boltzmann->SetLineColor(kMagenta);
      mT->SetLineColor(kGreen + 2);
      fProj[iMult][iPtTrigg][iReg]->GetYaxis()->SetMaxDigits(2);
      fProj[iMult][iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0, fProj[iMult][iPtTrigg][iReg]->GetMaximum() * 1.2);
      fProj[iMult][iPtTrigg][iReg]->GetXaxis()->SetRangeUser(0, 15);
      fProj[iMult][iPtTrigg][iReg]->DrawCopy();
      // fermiDir->DrawCopy("same");
      // levy->DrawCopy("same");
      // boltzmann->DrawCopy("same");
      // mT->DrawCopy("same");
      if (iMult == 0 && iReg == 0 && iPtTrigg == 0)
      {
        // leg->AddEntry(fermiDir, "Fermi-Dirac", "l");
        // leg->AddEntry(levy, "Levy-Tsallis", "l");
        // leg->AddEntry(boltzmann, "Boltzmann", "l");
        // leg->AddEntry(mT, "m_{T} scaling", "l");
      }
      // leg->Draw();

      pave[iMult][iPtTrigg][iReg] = new TPaveText();
      Plotter::SetPaveText(pave[iMult][iPtTrigg][iReg], 42, 0.05, 0, 0, 33, 0, 0.55, 0.97, 0.6, 0.95);
      pave[iMult][iPtTrigg][iReg]->AddText("ALICE, Work in Progress");
      pave[iMult][iPtTrigg][iReg]->AddText("pp, 13.6 TeV");
      pave[iMult][iPtTrigg][iReg]->AddText(Form("%s", multiplicityPave[iMult].Data()));
      pave[iMult][iPtTrigg][iReg]->AddText(Form("h-%s", finalNames[iPart].Data()));
      pave[iMult][iPtTrigg][iReg]->AddText(Form("%g < #font[52]{p}^{trigg}_{T} < %g GeV/#font[52]{c}", ptTriggBins[iPtTrigg], ptTriggBins[iPtTrigg + 1]));
      pave[iMult][iPtTrigg][iReg]->AddText("|#Delta#eta| < 1");
      pave[iMult][iPtTrigg][iReg]->AddText(paveRegions[iReg].Data());
      pave[iMult][iPtTrigg][iReg]->Draw("same");

      // yield=fProj[iMult][iPtTrigg][iReg]->Integral(0,8,"width");
      yield = 0.;
      erryield = 0.;
      Double_t erryieldSist = 0.;
      for (Int_t i = 1; i <= fProj[iMult][iPtTrigg][iReg]->GetNbinsX(); i++)
      {
        yield += fProj[iMult][iPtTrigg][iReg]->GetBinContent(i) * fProj[iMult][iPtTrigg][iReg]->GetBinWidth(i);
        erryield += pow(fProj[iMult][iPtTrigg][iReg]->GetBinError(i), 2) * pow(fProj[iMult][iPtTrigg][iReg]->GetBinWidth(i), 2);
        // assuming errors are fully correlated in pt
        erryieldSist += fProjSyst[iMult][iPtTrigg][iReg]->GetBinError(i) * fProj[iMult][iPtTrigg][iReg]->GetBinWidth(i);
      }
      erryield = sqrt(erryield);
      // cout << "Rel stat error of second pt interval: " << fProj[iMult][iPtTrigg][iReg]->GetBinError(2) / fProj[iMult][iPtTrigg][iReg]->GetBinContent(2) << endl;
      cout << "Rel stat error of yield: " << erryield / yield << endl;
      cout << "Rel sist error of yield: " << erryieldSist / yield << endl;
      cout << iMult << "___" << iPtTrigg << "_____" << yield << endl;
      // yield+=fermiDir->Integral(0,0.5);
      // yield+=fermiDir->Integral(6,100);

      histYield[iPtTrigg][iReg]->SetBinContent(9 - iMult, yield);
      histYield[iPtTrigg][iReg]->SetBinError(9 - iMult, erryield);
      histYieldSist[iPtTrigg][iReg]->SetBinContent(9 - iMult, yield);
      histYieldSist[iPtTrigg][iReg]->SetBinError(9 - iMult, erryieldSist);
      histYield[iPtTrigg][iReg]->GetXaxis()->SetBinLabel(iMult + 1, labels[iMult]);

      can[iMult][iPtTrigg][iReg]->SaveAs(Form("../../FitPt%s/Ptspectrum_%s_%s_ptTrigg%d.pdf", particleName[iPart].Data(), namesRegions[iReg].Data(), multiplicityNames[iMult].Data(), iPtTrigg));
      can[iMult][iPtTrigg][iReg]->Close();
      if (iMult == 0)
      {
        // canMult[iPtTrigg][iReg] = Plotter::CreateCanvas(Form("canMult%i%i", iPtTrigg, iReg));
        canMult[iPtTrigg][iReg] = new TCanvas(Form("canMult_pttrigg%i_ireg%i", iPtTrigg, iReg), Form("canMult_pttrigg%i_", iPtTrigg) + paveRegions[iReg], 1400, 800);
        canMult[iPtTrigg][iReg]->Divide(5, 2);
      }
      canMult[iPtTrigg][iReg]->cd(iMult + 1);
      canMult[iPtTrigg][iReg]->cd(iMult + 1)->SetLogy();
      fProj[iMult][iPtTrigg][iReg]->GetYaxis()->SetRangeUser(1e-8, fProj[iMult][iPtTrigg][iReg]->GetMaximum() * 5);
      fProj[iMult][iPtTrigg][iReg]->DrawCopy();
    }
  }
  TLegend *legPtTrigg = Plotter::CreateLegend(0.2, 0.3, 0.7, 0.9, 0.04);
  TCanvas *canYield = Plotter::CreateCanvas(Form("canY"));

  for (Int_t iPtTrigg = 0; iPtTrigg < 4; iPtTrigg++)
  {
    Plotter::SetHist(histYield[iPtTrigg][iReg], "", markers[iPtTrigg], colRegions[iReg][iPtTrigg], 1.);
    Plotter::SetHist(histYieldSist[iPtTrigg][iReg], "", markers[iPtTrigg], colRegions[iReg][iPtTrigg], 1.);
    // histYield[iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0.8 * histYield[0][iReg]->GetMinimum(), 1.6 * histYield[3][iReg]->GetMaximum());
    if (iReg == 0)
      histYield[iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0, 0.3);
    else if (iReg == 1)
      histYield[iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0, 0.3);
    else if (iReg == 2)
      histYield[iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0, 4);
    if (iPtTrigg == 0)
    {
      Plotter::SetHistAxes(histYield[iPtTrigg][iReg], "", "Y");
      histYield[iPtTrigg][iReg]->DrawCopy("ex0");
    }
    else
      histYield[iPtTrigg][iReg]->DrawCopy("same ex0");
    histYieldSist[iPtTrigg][iReg]->SetFillStyle(0);
    histYieldSist[iPtTrigg][iReg]->DrawCopy("same e2");
    legPtTrigg->AddEntry(histYield[iPtTrigg][iReg], Form("%g < #font[52]{p}^{trigg}_{T} < %g GeV/#font[52]{c}", ptTriggBins[iPtTrigg], ptTriggBins[iPtTrigg + 1]), "p");
  }

  TPaveText *paveYields = new TPaveText();
  Plotter::SetPaveText(paveYields, 42, 0.04, 0, 0, 33, 0, 0.7, 0.95, 0.6, 0.90);
  paveYields->AddText("ALICE, Work in Progress");
  paveYields->AddText("pp, 13.6 TeV");
  // paveYields->AddText(Form("h-%s, |#Delta#it{#eta}| < 1.1", finalNames[iPart].Data()));
  paveYields->AddText(Form("h-%s, |#it{#eta}_{trigg}| < 0.8, |#it{#eta}_{%s}| < 0.8", finalNames[iPart].Data(), finalNames[iPart].Data()));
  paveYields->AddText(Form("%s: |#Delta#it{#eta}| < 1.1 %s", paveRegions[iReg].Data(), PhiRegions[iReg].Data()));
  // paveYields->AddText(Form("|#it{#eta}_{trigg}| < 0.8, |#it{#eta}_{%s}| < 0.8", finalNames[iPart].Data()));
  paveYields->Draw("same");
  legPtTrigg->Draw();
  canYield->SaveAs(Form("../../YieldsVsMult_%s_%s.pdf", particleName[iPart].Data(), paveRegions[iReg].Data()));
  canYield->SaveAs(Form("../../YieldsVsMult_%s_%s.png", particleName[iPart].Data(), paveRegions[iReg].Data()));

  // RATIOS to MB
  TCanvas *canYieldToMB = Plotter::CreateCanvas(Form("canYToMB"));
  Int_t BinMB = 0;
  Float_t Err = 0;
  Float_t ErrSist = 0;
  for (Int_t iPtTrigg = 0; iPtTrigg < 4; iPtTrigg++)
  {
    histYieldToMB[iPtTrigg][0] = new TH1F(Form("histYieldNearToMB%d", iPtTrigg), "", 8, 0, 8);
    histYieldToMB[iPtTrigg][1] = new TH1F(Form("histYieldAwayToMB%d", iPtTrigg), "", 8, 0, 8);
    histYieldToMB[iPtTrigg][2] = new TH1F(Form("histYieldUEToMB%d", iPtTrigg), "", 8, 0, 8);
    histYieldSistToMB[iPtTrigg][0] = new TH1F(Form("histYieldSistNearToMB%d", iPtTrigg), "", 8, 0, 8);
    histYieldSistToMB[iPtTrigg][1] = new TH1F(Form("histYieldSistAwayToMB%d", iPtTrigg), "", 8, 0, 8);
    histYieldSistToMB[iPtTrigg][2] = new TH1F(Form("histYieldSistUEToMB%d", iPtTrigg), "", 8, 0, 8);
    for (Int_t b = 1; b <= histYieldToMB[iPtTrigg][0]->GetNbinsX(); b++)
    {
      histYieldToMB[iPtTrigg][iReg]->GetXaxis()->SetBinLabel(b, labels[b - 1]);
      BinMB = histYield[iPtTrigg][0]->GetNbinsX();
      histYieldToMB[iPtTrigg][iReg]->SetBinContent(b, histYield[iPtTrigg][iReg]->GetBinContent(b) / histYield[iPtTrigg][iReg]->GetBinContent(BinMB));
      histYieldSistToMB[iPtTrigg][iReg]->SetBinContent(b, histYieldToMB[iPtTrigg][iReg]->GetBinContent(b));
      Err = histYieldToMB[iPtTrigg][iReg]->GetBinContent(b) * sqrt(pow(histYield[iPtTrigg][iReg]->GetBinError(b) / histYield[iPtTrigg][iReg]->GetBinContent(b), 2) + pow(histYield[iPtTrigg][iReg]->GetBinError(BinMB) / histYield[iPtTrigg][iReg]->GetBinContent(BinMB), 2));
      ErrSist = histYieldToMB[iPtTrigg][iReg]->GetBinContent(b) * sqrt(pow(histYieldSist[iPtTrigg][iReg]->GetBinError(b) / histYieldSist[iPtTrigg][iReg]->GetBinContent(b), 2) + pow(histYieldSist[iPtTrigg][iReg]->GetBinError(BinMB) / histYieldSist[iPtTrigg][iReg]->GetBinContent(BinMB), 2));
      histYieldToMB[iPtTrigg][iReg]->SetBinError(b, Err);
      histYieldSistToMB[iPtTrigg][iReg]->SetBinError(b, ErrSist);
    }
    Plotter::SetHist(histYieldToMB[iPtTrigg][iReg], "", markers[iPtTrigg], colRegions[iReg][iPtTrigg], 1.);
    Plotter::SetHist(histYieldSistToMB[iPtTrigg][iReg], "", markers[iPtTrigg], colRegions[iReg][iPtTrigg], 1.);

    if (iPtTrigg == 0)
    {
      Plotter::SetHistAxes(histYieldToMB[iPtTrigg][iReg], "", "Y/Y_{MB}");
      histYieldToMB[iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0.7, 1.3);
      if (iReg == 2)
        histYieldToMB[iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0, 2);
      histYieldToMB[iPtTrigg][iReg]->DrawCopy("ex0");
    }
    else
      histYieldToMB[iPtTrigg][iReg]->DrawCopy("same ex0");
    histYieldSistToMB[iPtTrigg][iReg]->SetFillStyle(0);
    histYieldSistToMB[iPtTrigg][iReg]->DrawCopy("same e2");
  }

  // paveYields->Draw("same");
  legPtTrigg->Draw("same");
  canYieldToMB->SaveAs(Form("../../YieldsVsMultToMB_%s_%s.pdf", particleName[iPart].Data(), paveRegions[iReg].Data()));
  canYieldToMB->SaveAs(Form("../../YieldsVsMultToMB_%s_%s.png", particleName[iPart].Data(), paveRegions[iReg].Data()));

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

  for (Int_t iPtTrigg = 0; iPtTrigg < 4; iPtTrigg++)
  {
    canvasPtSpectra = new TCanvas(Form("canvasPtSpectra_pttrigg%i", iPtTrigg), Form("canvasPtSpectra_pttrigg%i", iPtTrigg), 700, 900);
    pad1 = new TPad("pad1", "pad1", 0, LLUpperPad, 1, 1); // xlow, ylow, xup, yup
    padL1 = new TPad("padL1", "padL1", 0, 0, 1, ULLowerPad);

    StylePad(pad1, 0.18, 0.01, 0.03, 0.);   // L, R, T, B
    StylePad(padL1, 0.18, 0.01, 0.02, 0.3); // L, R, T, B

    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendBorderSize(0);

    TLegend *legendAllMult = new TLegend(0.22, 0.03, 0.73, 0.28);
    legendAllMult->SetHeader("FT0M Multiplicity Percentile");
    legendAllMult->SetNColumns(3);
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
    LegendTitle->AddEntry("", Form("%s", paveRegions[iReg].Data()), "");

    LimSupSpectra = 99.99;
    if (iPtTrigg == 0)
    {
      if (iReg == 2)
      {
        LimInfSpectra = 0.4 * 1e-7;
        LimSupSpectra = 9999.99;
      }
      else{
        LimInfSpectra = 0.4 * 1e-6;
        LimSupSpectra = 9.99;
      }
    }
    else if (iPtTrigg == 1)
    {
      if (iReg == 2)
      {
        LimInfSpectra = 0.4 * 1e-5;
        LimSupSpectra = 9999.99;
      }
      else
        LimInfSpectra = 0.4 * 1e-5;
    }
    else if (iPtTrigg == 2)
    {
      if (iReg == 2)
      {
        LimInfSpectra = 0.4 * 1e-5;
        LimSupSpectra = 9999.99;
      }
      else
        LimInfSpectra = 0.4 * 1e-5;
    }
    else
    {
      if (iReg == 2)
      {
        LimInfSpectra = 0.4 * 1e-5;
        LimSupSpectra = 99999.99;
      }
      else
      {
        LimInfSpectra = 0.4 * 1e-4;
        LimSupSpectra = 999.99;
      }
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
    if (iReg == 2)
      ScaleFactorMB = pow(2, 10);
    else
      ScaleFactorMB = pow(2, 8);

    for (Int_t iMult = 0; iMult < 9; iMult++)
    {
      ScaleFactorFinal[iMult] = ScaleFactor[iMult];
      if (iMult == 0)
      {
        ColorMult[iMult] = ColorMB;
        MarkerMult[iMult] = MarkerMB;
        SizeMult[iMult] = SizeMB;
        ScaleFactorFinal[iMult] = ScaleFactorMB;
      }
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
      fProjScaled[iMult][iPtTrigg][iReg]->SetMarkerColor(ColorMult[iMult]);
      fProjScaled[iMult][iPtTrigg][iReg]->SetLineColor(ColorMult[iMult]);
      fProjScaled[iMult][iPtTrigg][iReg]->SetMarkerStyle(MarkerMult[iMult]);
      fProjScaled[iMult][iPtTrigg][iReg]->SetMarkerSize(SizeMult[iMult]);
      fProjSistScaled[iMult][iPtTrigg][iReg]->SetMarkerColor(ColorMult[iMult]);
      fProjSistScaled[iMult][iPtTrigg][iReg]->SetLineColor(ColorMult[iMult]);
      fProjSistScaled[iMult][iPtTrigg][iReg]->SetMarkerStyle(MarkerMult[iMult]);
      fProjSistScaled[iMult][iPtTrigg][iReg]->SetMarkerSize(SizeMult[iMult]);
      fProjScaled[iMult][iPtTrigg][iReg]->Draw("same ex0");
      fProjSistScaled[iMult][iPtTrigg][iReg]->SetFillStyle(0);
      fProjSistScaled[iMult][iPtTrigg][iReg]->Draw("same e2");
      sScaleFactorFinal[iMult] = Form(" (x2^{%i})", int(log2(ScaleFactorFinal[iMult])));
      if (ScaleFactorFinal[iMult] == 1)
        sScaleFactorFinal[iMult] = "";
      else if (ScaleFactorFinal[iMult] == 2)
        sScaleFactorFinal[iMult] = " (x2)";
      legendAllMult->AddEntry(fProjScaled[iMult][iPtTrigg][iReg], Form("%s", multiplicityPave[iMult].Data()) + sScaleFactorFinal[iMult] + " ", "pef");
    } // end loop on mult
    LegendTitle->Draw("");
    legendAllMult->Draw("");

    Int_t ChosenMult = 0; // MB
    TString TitleYSpectraRatio = Form("Ratio to %s", multiplicityPave[ChosenMult].Data());
    TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 10000, 0, 15.5);
    SetFont(hDummyRatio);
    StyleHistoYield(hDummyRatio, LimInfMultRatio, LimSupMultRatio, 1, 1, TitleXPt, TitleYSpectraRatio, "", 1, 1.15, YoffsetSpectraRatio);
    SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
    SetTickLength(hDummyRatio, tickXRatio, tickYRatio);
    hDummyRatio->GetXaxis()->SetRangeUser(0, UpRangePt[iPtTrigg]); // 15.5
    canvasPtSpectra->cd();
    padL1->Draw();
    padL1->cd();
    // gPad->SetLogy();
    hDummyRatio->Draw("same");

    for (Int_t iMult = 0; iMult < 9; iMult++)
    {
      fHistSpectrumStatMultRatio[iMult] = (TH1F *)fProj[iMult][iPtTrigg][iReg]->Clone(Form("fHistSpectrumStatRatio_%s_pttrigg%i", multiplicityNames[iMult].Data(), iPtTrigg));
      fHistSpectrumSistMultRatio[iMult] = (TH1F *)fProjSyst[iMult][iPtTrigg][iReg]->Clone(Form("fHistSpectrumSistRatio_%s_pttrigg%i", multiplicityNames[iMult].Data(), iPtTrigg));
      fHistSpectrumStatMultRatio[iMult]->Divide(fProj[ChosenMult][iPtTrigg][iReg]);
      fHistSpectrumSistMultRatio[iMult]->Divide(fProjSyst[ChosenMult][iPtTrigg][iReg]);
      ErrRatioCorr(fProj[iMult][iPtTrigg][iReg], fProj[ChosenMult][iPtTrigg][iReg], fHistSpectrumStatMultRatio[iMult], 0);
      ErrRatioCorr(fProjSyst[iMult][iPtTrigg][iReg], fProjSyst[ChosenMult][iPtTrigg][iReg], fHistSpectrumSistMultRatio[iMult], 0);
      fHistSpectrumStatMultRatio[iMult]->GetYaxis()->SetRangeUser(LimInfMultRatio, LimSupMultRatio);
      fHistSpectrumSistMultRatio[iMult]->GetYaxis()->SetRangeUser(LimInfMultRatio, LimSupMultRatio);
      fHistSpectrumStatMultRatio[iMult]->SetMarkerColor(ColorMult[iMult]);
      fHistSpectrumStatMultRatio[iMult]->SetLineColor(ColorMult[iMult]);
      fHistSpectrumStatMultRatio[iMult]->SetMarkerStyle(MarkerMult[iMult]);
      fHistSpectrumStatMultRatio[iMult]->SetMarkerSize(SizeMultRatio[iMult]);
      fHistSpectrumSistMultRatio[iMult]->SetMarkerColor(ColorMult[iMult]);
      fHistSpectrumSistMultRatio[iMult]->SetLineColor(ColorMult[iMult]);
      fHistSpectrumSistMultRatio[iMult]->SetMarkerStyle(MarkerMult[iMult]);
      fHistSpectrumSistMultRatio[iMult]->SetMarkerSize(SizeMultRatio[iMult]);
      if (iMult != ChosenMult)
      {
        fHistSpectrumStatMultRatio[iMult]->Draw("same ex0");
        fHistSpectrumSistMultRatio[iMult]->SetFillStyle(0);
        fHistSpectrumSistMultRatio[iMult]->Draw("same e2");
      }
    } // end loop on mult
    TString stringoutpdf = Form("../../Ptspectrum_%s_%s_ptTrigg%d", particleName[iPart].Data(), namesRegions[iReg].Data(), iPtTrigg);
    canvasPtSpectra->SaveAs(stringoutpdf + ".pdf");
    canvasPtSpectra->SaveAs(stringoutpdf + ".png");
    canvasPtSpectra->SaveAs(stringoutpdf + ".eps");
  } // end loop on pttrigg
}
