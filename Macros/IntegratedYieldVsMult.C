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

  TFile *file[nMultBins][nPtTriggBins];
  TFile *fileSyst;
  TH1F *fProj[nMultBins][nPtTriggBins][nRegions];
  TH1F *fProjRelSyst[nMultBins][nPtTriggBins][nRegions];
  TH1F *fProjSyst[nMultBins][nPtTriggBins][nRegions];
  TH1F *fProjRelSyst_apass4vsapass6[nMultBins][nPtTriggBins][nRegions];
  TH1F *fProjRelSyst_MCClosure[nMultBins][nPtTriggBins][nRegions];
  TH1F *fProjRelSyst_Final[nMultBins][nPtTriggBins][nRegions];
  TCanvas *can[nMultBins][nPtTriggBins][nRegions];
  TCanvas *canMult[nPtTriggBins][nRegions];
  TPaveText *pave[nMultBins][nPtTriggBins][nRegions];
  TLegend *leg = Plotter::CreateLegend(0.65, 0.95, 0.25, 0.5, 0.05);
  Double_t yield, erryield, erryieldSist, err;

  // const char *labels[11] = {"90-100%","80-90%","70-80%","60-70%","50-60%","40-50%","30-40%","20-30%","10-20%","0-10%","MB"};
  const char *labels[nMultBins] = {"70#minus100%", "50#minus70%", "40#minus50%", "30#minus40%", "20#minus30%", "10#minus20%", "1#minus10%", "0#minus1%", "0#minus100%"};
  TH1F *histYield[nPtTriggBins][nRegions];
  TH1F *histYieldToMB[nPtTriggBins][nRegions];
  TH1F *histYieldSist[nPtTriggBins][nRegions];
  TH1F *histYieldSistToMB[nPtTriggBins][nRegions];
  TH1F *histYieldToMBAllErrors[nPtTriggBins][nRegions];
  for (Int_t iPtTrigg = 0; iPtTrigg < nPtTriggBins; iPtTrigg++)
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

  TF1 *levy = new TF1("levy", "[nRegions]*x/([0]*[1])*(([0]-1)*([0]-2))/([0]*[1]+[2]*([0]-2))*TMath::Power(1+(TMath::Sqrt([2]*[2]+x*x)-[2])/([0]*[1]),-[0])", 0, 6);
  levy->SetParameter(0, 7);
  levy->SetParameter(1, 0.8);
  levy->FixParameter(2, particleMass[iPart]);
  levy->SetParameter(3, 0.03);

  TF1 *mT = new TF1("mT", "x*[0]*TMath::Exp(-TMath::Sqrt([2]*[2]+x*x)/[1])", 0, 6);
  mT->SetParameter(0, 0.01);
  mT->SetParameter(1, 0.7);
  mT->SetParLimits(1, 0.1, 10);
  mT->FixParameter(2, particleMass[iPart]);

  fileSyst = new TFile("../../Uncertainty.root", "");
  if (!fileSyst)
  {
    cout << "File not found" << endl;
    return;
  }

  for (Int_t iMult = 0; iMult < nMultBins; iMult++)
  {
    for (Int_t iPtTrigg = 0; iPtTrigg < nPtTriggBins; iPtTrigg++)
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
      fProjRelSyst[iMult][iPtTrigg][iReg]->SetName(Form("fhistRelSyst_%s_%s_%d", namesRegionsShort[iReg].Data(), multiplicityNames[iMult].Data(), iPtTrigg));
      // adding rel uncertainty associated with apass4 vs apass6
      fProjRelSyst_apass4vsapass6[iMult][iPtTrigg][iReg] = (TH1F *)fProjRelSyst[iMult][iPtTrigg][iReg]->Clone(Form("fhistRelSyst_apass4vsapass6_%s_%s_%d", namesRegionsShort[iReg].Data(), multiplicityNames[iMult].Data(), iPtTrigg));
      if (!fProjRelSyst_apass4vsapass6[iMult][iPtTrigg][iReg])
      {
        cout << "Histogram uncertainty pass4 vs pass6 not found" << endl;
        return;
      }
      for (Int_t b = 1; b <= fProjRelSyst_apass4vsapass6[iMult][iPtTrigg][iReg]->GetNbinsX(); b++)
      {
        if (iMult != 0)
          fProjRelSyst_apass4vsapass6[iMult][iPtTrigg][iReg]->SetBinContent(b, 0.02);
        else
          fProjRelSyst_apass4vsapass6[iMult][iPtTrigg][iReg]->SetBinContent(b, 0);
      }

      // adding rel uncertainty associated with closure test
      TFile *fileClosure = new TFile("ClosureUncertainty.root");
      if (!fileClosure)
      {
        cout << "File closure test uncertainties not found" << endl;
        return;
      }
      fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg] = (TH1F *)fileClosure->Get(Form("fhistRelSyst_MCClosure_%s_%s_%d", namesRegionsShort[iReg].Data(), multiplicityNames[6].Data(), iPtTrigg));
      if (!fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg])
      {
        cout << "Histogram closure test uncertainty not found" << endl;
        return;
      }

      // put all uncertainties together
      fProjRelSyst_Final[iMult][iPtTrigg][iReg] = (TH1F *)fProjRelSyst[iMult][iPtTrigg][iReg]->Clone(Form("fhistRelSyst_Final_%s_%s_%d", namesRegionsShort[iReg].Data(), multiplicityNames[iMult].Data(), iPtTrigg));
      for (Int_t b = 1; b <= fProjRelSyst_Final[iMult][iPtTrigg][iReg]->GetNbinsX(); b++)
      {
        if (iMult != 0)
          fProjRelSyst_Final[iMult][iPtTrigg][iReg]->SetBinContent(b, sqrt(pow(fProjRelSyst[iMult][iPtTrigg][iReg]->GetBinContent(b), 2) + pow(fProjRelSyst_apass4vsapass6[iMult][iPtTrigg][iReg]->GetBinContent(b), 2) + pow(fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->GetBinContent(b), 2)));
        else
          fProjRelSyst_Final[iMult][iPtTrigg][iReg]->SetBinContent(b, fProjRelSyst[iMult][iPtTrigg][iReg]->GetBinContent(b));
      }

      can[iMult][iPtTrigg][iReg] = Plotter::CreateCanvas(Form("can%i%i%i", iMult, iPtTrigg, iReg));
      if (iReg == 2)
      {
        fProj[iMult][iPtTrigg][iReg]->Scale(50);
        fProjSyst[iMult][iPtTrigg][iReg]->Scale(50);
      }

      for (Int_t b = 1; b <= fProj[iMult][iPtTrigg][iReg]->GetNbinsX(); b++)
      {
        fProjSyst[iMult][iPtTrigg][iReg]->SetBinError(b, fProjRelSyst_Final[iMult][iPtTrigg][iReg]->GetBinContent(b) * fProjSyst[iMult][iPtTrigg][iReg]->GetBinContent(b));
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
      pave[iMult][iPtTrigg][iReg]->AddText("ALICE, Preliminary");
      pave[iMult][iPtTrigg][iReg]->AddText("pp, 13.6 TeV");
      pave[iMult][iPtTrigg][iReg]->AddText(Form("%s", multiplicityPave[iMult].Data()));
      pave[iMult][iPtTrigg][iReg]->AddText(Form("h#minus%s correlation", finalNames[iPart].Data()));
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

      histYield[iPtTrigg][iReg]->SetBinContent(nMultBins - iMult, yield);
      histYield[iPtTrigg][iReg]->SetBinError(nMultBins - iMult, erryield);
      histYieldSist[iPtTrigg][iReg]->SetBinContent(nMultBins - iMult, yield);
      histYieldSist[iPtTrigg][iReg]->SetBinError(nMultBins - iMult, erryieldSist);
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

  // PLOT: uncertainties in 0-1%
  TLegend * legendUnc = Plotter::CreateLegend(0.15, 0.3, 0.65, 0.85, 0.04);
  for (Int_t iPtTrigg = 0; iPtTrigg < nPtTriggBins; iPtTrigg++)
  {
    for (Int_t iMult = 0; iMult < nMultBins; iMult++)
    {
      if (iMult != 1)
        continue;
      TCanvas *canSyst = new TCanvas(Form("canSyst%i%i", iMult, iPtTrigg), Form("canSyst%i%i", iMult, iPtTrigg), 1200, 800);
      fProjRelSyst[iMult][iPtTrigg][iReg]->SetLineColor(kBlue);
      fProjRelSyst_apass4vsapass6[iMult][iPtTrigg][iReg]->SetLineColor(kRed);
      fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->SetLineColor(kGreen + 2);
      fProjRelSyst_Final[iMult][iPtTrigg][iReg]->SetLineColor(kBlack);
      fProjRelSyst_Final[iMult][iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0, 0.5);
      fProjRelSyst_Final[iMult][iPtTrigg][iReg]->Draw();
      fProjRelSyst_apass4vsapass6[iMult][iPtTrigg][iReg]->Draw("same");
      fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->Draw("same hist");
      fProjRelSyst[iMult][iPtTrigg][iReg]->Draw("same");
      if (iPtTrigg == 0)
      {
        legendUnc->AddEntry(fProjRelSyst_Final[iMult][iPtTrigg][iReg], "Total", "l");
        legendUnc->AddEntry(fProjRelSyst_apass4vsapass6[iMult][iPtTrigg][iReg], "apass4 vs apass6", "l");
        legendUnc->AddEntry(fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg], "Closure test", "l");
        legendUnc->AddEntry(fProjRelSyst[iMult][iPtTrigg][iReg], "Other", "l");
      }
      legendUnc->Draw();
      canSyst->SaveAs(Form("../../Syst_%s_%s_ptTrigg%d.pdf", particleName[iPart].Data(), namesRegions[iReg].Data(), iPtTrigg));
      canSyst->SaveAs(Form("../../Syst_%s_%s_ptTrigg%d.png", particleName[iPart].Data(), namesRegions[iReg].Data(), iPtTrigg));
    }
  }

  TLegend *legPtTrigg = Plotter::CreateLegend(0.2, 0.3, 0.7, 0.9, 0.04);
  TCanvas *canYield = Plotter::CreateCanvas(Form("canY"));

  for (Int_t iPtTrigg = 0; iPtTrigg < nPtTriggBins; iPtTrigg++)
  {
    Plotter::SetHist(histYield[iPtTrigg][iReg], "", markers[iPtTrigg], colRegions[iReg][iPtTrigg], 1.);
    Plotter::SetHist(histYieldSist[iPtTrigg][iReg], "", markers[iPtTrigg], colRegions[iReg][iPtTrigg], 1.);
    // histYield[iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0.8 * histYield[0][iReg]->GetMinimum(), 1.6 * histYield[nRegions][iReg]->GetMaximum());
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
  paveYields->AddText("ALICE, Preliminary");
  paveYields->AddText("pp, 13.6 TeV");
  paveYields->AddText(Form("h#minus%s correlation, |#it{#eta}_{trigg}| < 0.8, |#it{#eta}_{%s}| < 0.8", finalNames[iPart].Data(), finalNames[iPart].Data()));
  paveYields->AddText(Form("%s: |#Delta#it{#eta}| < 1.1 %s", paveRegions[iReg].Data(), PhiRegions[iReg].Data()));
  paveYields->Draw("same");
  legPtTrigg->Draw();
  canYield->SaveAs(Form("../../YieldsVsMult_%s_%s.pdf", particleName[iPart].Data(), paveRegions[iReg].Data()));
  canYield->SaveAs(Form("../../YieldsVsMult_%s_%s.png", particleName[iPart].Data(), paveRegions[iReg].Data()));

  // RATIOS to MB
  Int_t ChosenMultYRatio = 0;
  TCanvas *canYieldToMB = Plotter::CreateCanvas(Form("canYToMB"));
  Int_t BinMB = 0;
  Float_t Err = 0;
  Float_t ErrSist = 0;
  TF1 *pol0[nPtTriggBins];
  for (Int_t iPtTrigg = 0; iPtTrigg < nPtTriggBins; iPtTrigg++)
  {
    // 70-100% not included!
    histYieldToMB[iPtTrigg][0] = new TH1F(Form("histYieldNearToMB%d", iPtTrigg), "", 7, 0, 7);
    histYieldToMB[iPtTrigg][1] = new TH1F(Form("histYieldAwayToMB%d", iPtTrigg), "", 7, 0, 7);
    histYieldToMB[iPtTrigg][2] = new TH1F(Form("histYieldUEToMB%d", iPtTrigg), "", 7, 0, 7);
    histYieldSistToMB[iPtTrigg][0] = new TH1F(Form("histYieldSistNearToMB%d", iPtTrigg), "", 7, 0, 7);
    histYieldSistToMB[iPtTrigg][1] = new TH1F(Form("histYieldSistAwayToMB%d", iPtTrigg), "", 7, 0, 7);
    histYieldSistToMB[iPtTrigg][2] = new TH1F(Form("histYieldSistUEToMB%d", iPtTrigg), "", 7, 0, 7);
    for (Int_t b = 1; b <= histYieldToMB[iPtTrigg][0]->GetNbinsX(); b++)
    {
      // histYieldToMB[iPtTrigg][iReg]->GetXaxis()->SetBinLabel(b, labels[b - 1]);
      histYieldToMB[iPtTrigg][iReg]->GetXaxis()->SetBinLabel(b, labels[b]); // skipping 70-100%
      if (ChosenMultYRatio == 0)
        BinMB = histYield[iPtTrigg][0]->GetNbinsX();
      else
      {
        return; // not implemented if we skip 70-100%
        BinMB = nMultBins - ChosenMultYRatio;
      }
      Int_t Shift = 1; // 0 to start from 70-100%
      histYieldToMB[iPtTrigg][iReg]->SetBinContent(b, histYield[iPtTrigg][iReg]->GetBinContent(b + Shift) / histYield[iPtTrigg][iReg]->GetBinContent(BinMB));
      histYieldSistToMB[iPtTrigg][iReg]->SetBinContent(b, histYieldToMB[iPtTrigg][iReg]->GetBinContent(b));
      Err = histYieldToMB[iPtTrigg][iReg]->GetBinContent(b) * sqrt(pow(histYield[iPtTrigg][iReg]->GetBinError(b + Shift) / histYield[iPtTrigg][iReg]->GetBinContent(b + Shift), 2) + pow(histYield[iPtTrigg][iReg]->GetBinError(BinMB) / histYield[iPtTrigg][iReg]->GetBinContent(BinMB), 2));
      ErrSist = histYieldToMB[iPtTrigg][iReg]->GetBinContent(b) * sqrt(pow(histYieldSist[iPtTrigg][iReg]->GetBinError(b + Shift) / histYieldSist[iPtTrigg][iReg]->GetBinContent(b + Shift), 2) + pow(histYieldSist[iPtTrigg][iReg]->GetBinError(BinMB) / histYieldSist[iPtTrigg][iReg]->GetBinContent(BinMB), 2));
      // Err = histYield[iPtTrigg][iReg]->GetBinError(b) / histYield[iPtTrigg][iReg]->GetBinContent(BinMB);
      // ErrSist = histYieldSist[iPtTrigg][iReg]->GetBinError(b) / histYieldSist[iPtTrigg][iReg]->GetBinContent(BinMB);
      histYieldToMB[iPtTrigg][iReg]->SetBinError(b, Err);
      histYieldSistToMB[iPtTrigg][iReg]->SetBinError(b, ErrSist);
    }
    Plotter::SetHist(histYieldToMB[iPtTrigg][iReg], "", markers[iPtTrigg], colRegions[iReg][iPtTrigg], 1.);
    Plotter::SetHist(histYieldSistToMB[iPtTrigg][iReg], "", markers[iPtTrigg], colRegions[iReg][iPtTrigg], 1.);

    if (iPtTrigg == 0)
    {
      Plotter::SetHistAxes(histYieldToMB[iPtTrigg][iReg], "", Form("Y/Y_%s", multiplicityPave[ChosenMultYRatio].Data()));
      histYieldToMB[iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0.7, 1.3);
      // histYieldToMB[iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0.7, 2);
      if (iReg == 2)
        histYieldToMB[iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0, 2);
      histYieldToMB[iPtTrigg][iReg]->DrawCopy("ex0");
    }
    else
      histYieldToMB[iPtTrigg][iReg]->DrawCopy("same ex0");
    histYieldSistToMB[iPtTrigg][iReg]->SetFillStyle(0);
    histYieldSistToMB[iPtTrigg][iReg]->DrawCopy("same e2");

    histYieldToMBAllErrors[iPtTrigg][iReg] = (TH1F *)histYieldToMB[iPtTrigg][iReg]->Clone(Form("histYieldToMBAllErrors_%i_%i", iPtTrigg, iReg));
    for (Int_t b = 1; b <= histYieldToMB[iPtTrigg][iReg]->GetNbinsX(); b++)
    {
      histYieldToMBAllErrors[iPtTrigg][iReg]->SetBinError(b, sqrt(pow(histYieldToMB[iPtTrigg][iReg]->GetBinError(b), 2) + pow(histYieldSistToMB[iPtTrigg][iReg]->GetBinError(b), 2)));
    }
    pol0[iPtTrigg] = new TF1(Form("plo0_%i", iPtTrigg), "pol0", 0, 8);
    pol0[iPtTrigg]->SetLineColor(colRegions[iReg][iPtTrigg]);
    cout << "\n\n************" << endl;
    cout << "Fit of " << iPtTrigg << " pt trigger bin" << endl;
    histYieldToMBAllErrors[iPtTrigg][iReg]->Fit(pol0[iPtTrigg], "R0");
    cout << "Chi /NDF " << pol0[iPtTrigg]->GetChisquare() << "/ " << pol0[iPtTrigg]->GetNDF() << " = " << pol0[iPtTrigg]->GetChisquare() / pol0[iPtTrigg]->GetNDF() << endl;
  }

  // paveYields->Draw("same");
  legPtTrigg->Draw("same");
  canYieldToMB->SaveAs(Form("../../YieldsVsMultToMB_%s_%s.pdf", particleName[iPart].Data(), paveRegions[iReg].Data()));
  canYieldToMB->SaveAs(Form("../../YieldsVsMultToMB_%s_%s.png", particleName[iPart].Data(), paveRegions[iReg].Data()));

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
    LegendTitle->AddEntry("", "#bf{ALICE Preliminary}", "");
    // option 1
    // LegendTitle->AddEntry("", "Run 3, pp #sqrt{#it{s}} = 13.6 TeV", "");
    // LegendTitle->AddEntry("", Form("h#minus%s correlation, %.0f < #it{p}_{T}^{trigg} < %.0f GeV/#it{c}", finalNames[iPart].Data(), ptTriggBins[iPtTrigg], ptTriggBins[iPtTrigg + 1]), "");
    // LegendTitle->AddEntry("", Form("|#it{#eta}^{trigg}| < 0.8, |#it{#eta}^{%s}| < 0.8", finalNames[iPart].Data()), "");
    // LegendTitle->AddEntry("", Form("%s: |#Delta#it{#eta}| < 1.1 %s", paveRegions[iReg].Data(), PhiRegions[iReg].Data()));
    // option 2
    LegendTitle->AddEntry("", Form("pp #sqrt{#it{s}} = 13.6 TeV, h#minus%s correlation", finalNames[iPart].Data()), "");
    LegendTitle->AddEntry("", Form("%.0f < #it{p}_{T}^{trigg} < %.0f GeV/#it{c}, |#it{#eta}^{trigg}| < 0.8, |#it{#eta}^{%s}| < 0.8", ptTriggBins[iPtTrigg], ptTriggBins[iPtTrigg + 1], finalNames[iPart].Data()), "");
    LegendTitle->AddEntry("", Form("%s: |#Delta#it{#eta}| < 1.1 %s", paveRegions[iReg].Data(), PhiRegions[iReg].Data()), "");

    // LimSupSpectra = 9999;
    LimSupSpectra = 30000;
    if (iPtTrigg == 0)
    {
      if (iReg == 2)
      {
        LimInfSpectra = 0.4 * 1e-7;
        LimSupSpectra = 99999.99;
      }
      else
        LimInfSpectra = 0.2 * 1e-5; // 0.4 * 1e-6 if not prob. density
    }
    else if (iPtTrigg == 1)
    {
      if (iReg == 2)
      {
        LimInfSpectra = 0.4 * 1e-6;
        LimSupSpectra = 999999;
      }
      else
        LimInfSpectra = 0.4 * 1e-4; // 0.4 * 1e-5 if not prob. density
    }
    else if (iPtTrigg == 2)
    {
      if (iReg == 2)
      {
        LimInfSpectra = 0.4 * 1e-6;
        LimSupSpectra = 999999;
      }
      else
        LimInfSpectra = 0.4 * 1e-4; // 0.4 * 1e-5 if not prob. density
    }
    else
    {
      if (iReg == 2)
      {
        LimInfSpectra = 0.4 * 1e-5;
        LimSupSpectra = 999999;
      }
      else
        LimInfSpectra = 0.4 * 1e-3; // 0.4 * 1e-4 if not prob. density
    }
    TH1F *hDummy = new TH1F("hDummy", "hDummy", 10000, 0, UpRangePt[iPtTrigg] + OffsetX[iPtTrigg]);
    for (Int_t i = 1; i <= hDummy->GetNbinsX(); i++)
      hDummy->SetBinContent(i, 1e-12);
    canvasPtSpectra->cd();
    SetFont(hDummy);
    StyleHistoYield(hDummy, LimInfSpectra, LimSupSpectra, 1, 1, TitleXPt, TitleYProbDensity, "", 1, 1.15, 1.6);
    SetHistoTextSize(hDummy, xTitle, xLabel, xOffset, xLabelOffset, yTitle, yLabel, yOffset, yLabelOffset);
    SetTickLength(hDummy, tickX, tickY);
    pad1->Draw();
    pad1->cd();
    gPad->SetLogy();
    gStyle->SetOptStat(0);
    hDummy->Draw("same");
    ScaleFactorMB = pow(2, 8); // 10 if not prob. density

    for (Int_t iMult = 0; iMult < nMultBins; iMult++)
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
      Float_t ProbDensityScaleFactor = histYield[iPtTrigg][iReg]->GetBinContent(nMultBins); // 0-100% integrated yields
      fProjScaled[iMult][iPtTrigg][iReg]->Scale(1. / ProbDensityScaleFactor);
      fProjSistScaled[iMult][iPtTrigg][iReg]->Scale(1. / ProbDensityScaleFactor);
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
      sScaleFactorFinal[iMult] = Form(" (x2^{%i})", int(log2(ScaleFactorFinal[iMult])));
      if (ScaleFactorFinal[iMult] == 1)
        sScaleFactorFinal[iMult] = "";
      else if (ScaleFactorFinal[iMult] == 2)
        sScaleFactorFinal[iMult] = " (x2)";
      if (iMult != (nMultBins - 1))
      {
        fProjScaled[iMult][iPtTrigg][iReg]->Draw("same ex0");
        fProjSistScaled[iMult][iPtTrigg][iReg]->SetFillStyle(0);
        fProjSistScaled[iMult][iPtTrigg][iReg]->Draw("same e2");
        legendAllMult->AddEntry(fProjScaled[iMult][iPtTrigg][iReg], Form("%s", multiplicityPave[iMult].Data()) + sScaleFactorFinal[iMult] + " ", "pef");
      }
    } // end loop on mult
    LegendTitle->Draw("");
    legendAllMult->Draw("");

    Int_t ChosenMult = 0; // MB
    TString TitleYSpectraRatio = Form("Ratio to %s", multiplicityPave[ChosenMult].Data());
    TH1F *hDummyRatio = new TH1F("hDummyRatio", "hDummyRatio", 10000, 0, UpRangePt[iPtTrigg] + OffsetX[iPtTrigg]);
    SetFont(hDummyRatio);
    StyleHistoYield(hDummyRatio, LimInfMultRatio, LimSupMultRatio, 1, 1, TitleXPt, TitleYSpectraRatio, "", 1, 1.15, YoffsetSpectraRatio);
    SetHistoTextSize(hDummyRatio, xTitleR, xLabelR, xOffsetR, xLabelOffsetR, yTitleR, yLabelR, yOffsetR, yLabelOffsetR);
    SetTickLength(hDummyRatio, tickXRatio, tickYRatio);
    canvasPtSpectra->cd();
    padL1->Draw();
    padL1->cd();
    gPad->SetLogy();
    hDummyRatio->Draw("same");

    for (Int_t iMult = 0; iMult < nMultBins; iMult++)
    {
      fHistSpectrumStatMultRatio[iMult] = (TH1F *)fProj[iMult][iPtTrigg][iReg]->Clone(Form("fHistSpectrumStatRatio_%s_pttrigg%i", multiplicityNames[iMult].Data(), iPtTrigg));
      fHistSpectrumSistMultRatio[iMult] = (TH1F *)fProjSyst[iMult][iPtTrigg][iReg]->Clone(Form("fHistSpectrumSistRatio_%s_pttrigg%i", multiplicityNames[iMult].Data(), iPtTrigg));
      fHistSpectrumStatMultRatio[iMult]->Divide(fProj[ChosenMult][iPtTrigg][iReg]);
      fHistSpectrumSistMultRatio[iMult]->Divide(fProjSyst[ChosenMult][iPtTrigg][iReg]);
      ErrRatioCorr(fProj[iMult][iPtTrigg][iReg], fProj[ChosenMult][iPtTrigg][iReg], fHistSpectrumStatMultRatio[iMult], 0);
      // ErrRatioCorr(fProjSyst[iMult][iPtTrigg][iReg], fProjSyst[ChosenMult][iPtTrigg][iReg], fHistSpectrumSistMultRatio[iMult], 0);
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
      if (iMult != ChosenMult && iMult != (nMultBins - 1)) // do not plot 0-100% and 70-100%
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

  TString stringoutroot = Form("../../YieldsPtIntegrated_%s_%s.root", particleName[iPart].Data(), namesRegions[iReg].Data());
  TFile *outputFile = new TFile(stringoutroot, "RECREATE");
  for (Int_t iPtTrigg = 0; iPtTrigg < nPtTriggBins; iPtTrigg++)
  {
    histYieldSist[iPtTrigg][iReg]->Write();
    histYield[iPtTrigg][iReg]->Write();
    histYieldSistToMB[iPtTrigg][iReg]->Write();
    histYieldToMB[iPtTrigg][iReg]->Write();
    for (Int_t iMult = 0; iMult < nMultBins; iMult++)
    {
      fProj[iMult][iPtTrigg][iReg]->Write();
      fProjSyst[iMult][iPtTrigg][iReg]->Write();
    }
  }
  outputFile->Close();
  cout << "\nI have created the output file: " << stringoutroot << endl;

  cout << "\n\n\e[35mWARNING: Uncertainty on closure test taken from mult 40-50\% for all multiplicity classes! Hardcoded in macro\n\n " << endl;
}
