#include "Plotter.h"
#include "TFile.h"
#include "Style.h"
#include "ErrRatioCorr.C"

void ClosureUncertainty()
{
  TFile *fileClosure;
  TFile *fileClosureMB;
  TH1F *fProj_MCClosure[nMultBins][nPtTriggBins][nRegions];
  TH1F *fProj_MCClosureMB[nPtTriggBins][nRegions];
  TH1F *fProjRelSyst_MCClosure_Old[nMultBins][nPtTriggBins][nRegions];
  TH1F *fProjRelSyst_MCClosureMB_Old[nPtTriggBins][nRegions];
  TH1F *fProjRelSyst_MCClosure[nMultBins][nPtTriggBins][nRegions];
  TLegend *legendMult;

  for (Int_t iPtTrigg = 0; iPtTrigg < nPtTriggBins; iPtTrigg++)
  {
    TCanvas *can = new TCanvas(Form("canMultComp_pttrigg%i", iPtTrigg), Form("canMultComp_pttrigg%i", iPtTrigg), 1500, 500);

    can->Divide(3, 1);
    for (Int_t iReg = 0; iReg < nRegions; iReg++)
    {
      legendMult = Plotter::CreateLegend(0.2, 0.3, 0.7, 0.9, 0.04);
      fileClosureMB = new TFile(Form("../../ClosureUncer_Trigg%d.root", iPtTrigg));
      if (!fileClosureMB)
      {
        cout << "File closure test MB not found" << endl;
        return;
      }
      fProj_MCClosureMB[iPtTrigg][iReg] = (TH1F *)fileClosureMB->Get(Form("fHistRatio%s", namesRegionsShort[iReg].Data()));
      fProj_MCClosureMB[iPtTrigg][iReg]->SetName(Form("fhistRelSyst_MCClosureMB_%s_%d", namesRegionsShort[iReg].Data(), iPtTrigg));
      if (!fProj_MCClosureMB[iPtTrigg][iReg])
      {
        cout << "Histogram closure test MB uncertainty not found" << endl;
        return;
      }
      // let's get reco/gen = f from (reco - gen)/gen
      for (Int_t b = 1; b <= fProj_MCClosureMB[iPtTrigg][iReg]->GetNbinsX(); b++)
      {
        fProj_MCClosureMB[iPtTrigg][iReg]->SetBinContent(b, fProj_MCClosureMB[iPtTrigg][iReg]->GetBinContent(b) + 1);
      }
      fProj_MCClosureMB[iPtTrigg][iReg]->Smooth();

      // old way of computing uncertanty for MB
      fProjRelSyst_MCClosureMB_Old[iPtTrigg][iReg] = (TH1F *)fProj_MCClosureMB[iPtTrigg][iReg]->Clone(Form("fhistRelSyst_MCClosureMB_Old%s_%d", namesRegionsShort[iReg].Data(), iPtTrigg));
      for (Int_t b = 1; b <= fProj_MCClosureMB[iPtTrigg][iReg]->GetNbinsX(); b++)
      {
        fProjRelSyst_MCClosureMB_Old[iPtTrigg][iReg]->SetBinContent(b, TMath::Abs(1 - fProj_MCClosureMB[iPtTrigg][iReg]->GetBinContent(b)));
      }
      fProjRelSyst_MCClosureMB_Old[iPtTrigg][iReg]->Smooth();
      // end old way

      for (Int_t iMult = 0; iMult < nMultBins; iMult++)
      {
        if (iMult != 2 && iMult != 6)
          continue; // keep only 1-10%, 40-50% and 0-100%
        fileClosure = new TFile(Form("../../ClosureUncer_Trigg%d_%s.root", iPtTrigg, multiplicityNamesShort[iMult].Data()));
        if (!fileClosure)
        {
          cout << "File closure test not found" << endl;
          return;
        }
        fProj_MCClosure[iMult][iPtTrigg][iReg] = (TH1F *)fileClosure->Get(Form("fHistRatio%s", namesRegionsShort[iReg].Data()));
        fProj_MCClosure[iMult][iPtTrigg][iReg]->SetName(Form("fhistRelSyst_MCClosure_%s_%s_%d", namesRegionsShort[iReg].Data(), multiplicityNames[iMult].Data(), iPtTrigg));
        if (!fProj_MCClosure[iMult][iPtTrigg][iReg])
        {
          cout << "Histogram closure test uncertainty not found" << endl;
          return;
        }

        // let's get reco/gen = f from (reco - gen)/gen
        for (Int_t b = 1; b <= fProj_MCClosure[iMult][iPtTrigg][iReg]->GetNbinsX(); b++)
        {
          fProj_MCClosure[iMult][iPtTrigg][iReg]->SetBinContent(b, fProj_MCClosure[iMult][iPtTrigg][iReg]->GetBinContent(b) + 1);
        }
        fProj_MCClosure[iMult][iPtTrigg][iReg]->Smooth();

        // new way of getting uncertianty
        fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg] = (TH1F *)fProj_MCClosure[iMult][iPtTrigg][iReg]->Clone(Form("fhistRelSyst_MCClosure_%s_%s_%d", namesRegionsShort[iReg].Data(), multiplicityNames[iMult].Data(), iPtTrigg));
        fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->Divide(fProj_MCClosureMB[iPtTrigg][iReg]);
        for (Int_t b = 1; b <= fProj_MCClosure[iMult][iPtTrigg][iReg]->GetNbinsX(); b++)
        {
          fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->SetBinContent(b, TMath::Abs(1 - fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->GetBinContent(b)));
        }
        fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->Smooth();
        // end new way

        // old way of computing uncertainty
        fProjRelSyst_MCClosure_Old[iMult][iPtTrigg][iReg] = (TH1F *)fProj_MCClosure[iMult][iPtTrigg][iReg]->Clone(Form("fhistRelSyst_MCClosure_Old%s_%s_%d", namesRegionsShort[iReg].Data(), multiplicityNames[iMult].Data(), iPtTrigg));
        for (Int_t b = 1; b <= fProj_MCClosure[iMult][iPtTrigg][iReg]->GetNbinsX(); b++)
        {
          fProjRelSyst_MCClosure_Old[iMult][iPtTrigg][iReg]->SetBinContent(b, TMath::Abs(1 - fProj_MCClosure[iMult][iPtTrigg][iReg]->GetBinContent(b)));
        }
        fProjRelSyst_MCClosure_Old[iMult][iPtTrigg][iReg]->Smooth();
        // end old way

        // comparison between closure tests in different multiplicity classes
        fProj_MCClosure[iMult][iPtTrigg][iReg]->SetLineColor(ColorMult[iMult]);
        fProj_MCClosureMB[iPtTrigg][iReg]->SetLineColor(kBlack);
        fProj_MCClosure[iMult][iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0.5, 1.5);
        Plotter::SetHist(fProj_MCClosure[iMult][iPtTrigg][iReg], "", 20, ColorMult[iMult], 1.);
        Plotter::SetHist(fProj_MCClosureMB[iPtTrigg][iReg], "", 20, kBlack, 1.);
        fProj_MCClosure[iMult][iPtTrigg][iReg]->GetYaxis()->SetTitle("rec/gen");
        fProj_MCClosure[iMult][iPtTrigg][iReg]->GetXaxis()->SetTitleSize(0.05);
        fProj_MCClosure[iMult][iPtTrigg][iReg]->GetYaxis()->SetTitleSize(0.05);
        fProj_MCClosure[iMult][iPtTrigg][iReg]->GetXaxis()->SetLabelSize(0.05);
        fProj_MCClosure[iMult][iPtTrigg][iReg]->GetYaxis()->SetLabelSize(0.05);
        fProj_MCClosureMB[iPtTrigg][iReg]->GetYaxis()->SetTitle("rec/gen");
        fProj_MCClosureMB[iPtTrigg][iReg]->GetXaxis()->SetTitleSize(0.05);
        fProj_MCClosureMB[iPtTrigg][iReg]->GetYaxis()->SetTitleSize(0.05);
        fProj_MCClosureMB[iPtTrigg][iReg]->GetXaxis()->SetLabelSize(0.05);
        fProj_MCClosureMB[iPtTrigg][iReg]->GetYaxis()->SetLabelSize(0.05);

        can->cd(iReg + 1);
        fProj_MCClosure[iMult][iPtTrigg][iReg]->Draw("same");
        fProj_MCClosureMB[iPtTrigg][iReg]->Draw("same");

        legendMult->AddEntry(fProj_MCClosure[iMult][iPtTrigg][iReg], Form("%s %s", namesRegionsShort[iReg].Data(), multiplicityPave[iMult].Data()), "l");
        if (iMult == 2)
          legendMult->AddEntry(fProj_MCClosureMB[iPtTrigg][iReg], Form("%s MB", namesRegionsShort[iReg].Data()), "l");
      }
      legendMult->Draw();
    }
  }

  for (Int_t iPtTrigg = 0; iPtTrigg < nPtTriggBins; iPtTrigg++)
  {
    TCanvas *canUnc = new TCanvas(Form("canUncer_pttrigg%i", iPtTrigg), Form("canUncer_pttrigg%i", iPtTrigg), 1500, 500);
    canUnc->Divide(3, 1);
    for (Int_t iReg = 0; iReg < nRegions; iReg++)
    {
      for (Int_t iMult = 0; iMult < nMultBins; iMult++)
      {
        if (iMult != 2 && iMult != 6)
          continue; // keep only 1-10%, 40-50% and 0-100%
        canUnc->cd(iReg + 1);
        fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->GetYaxis()->SetTitle("Rel. uncertainty");
        fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->GetXaxis()->SetTitleSize(0.05);
        fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->GetYaxis()->SetTitleSize(0.05);
        fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->GetXaxis()->SetLabelSize(0.05);
        fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->GetYaxis()->SetLabelSize(0.05);
        fProjRelSyst_MCClosureMB_Old[iPtTrigg][iReg]->GetYaxis()->SetTitle("Rel. uncertainty");
        fProjRelSyst_MCClosureMB_Old[iPtTrigg][iReg]->GetXaxis()->SetTitleSize(0.05);
        fProjRelSyst_MCClosureMB_Old[iPtTrigg][iReg]->GetYaxis()->SetTitleSize(0.05);
        fProjRelSyst_MCClosureMB_Old[iPtTrigg][iReg]->GetXaxis()->SetLabelSize(0.05);
        fProjRelSyst_MCClosureMB_Old[iPtTrigg][iReg]->GetYaxis()->SetLabelSize(0.05);
        fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->SetLineColor(ColorMult[iMult]);
        fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->GetYaxis()->SetRangeUser(0, 0.5);
        fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->Draw("hist same");
        fProjRelSyst_MCClosureMB_Old[iPtTrigg][iReg]->Draw("same hist");
      }
      legendMult->Draw();
    }
  }

  TFile *fileOut = new TFile("ClosureUncertainty.root", "recreate");
  for (Int_t iPtTrigg = 0; iPtTrigg < nPtTriggBins; iPtTrigg++)
  {
    for (Int_t iReg = 0; iReg < nRegions; iReg++)
    {
      for (Int_t iMult = 0; iMult < nMultBins; iMult++)
      {
        if (iMult != 2 && iMult != 6)
          continue; // keep only 1-10%, 40-50%
        fProjRelSyst_MCClosure[iMult][iPtTrigg][iReg]->Write();
      }
    }
  }
  fileOut->Close();
  cout << "File ClosureUncertainty.root created" << endl;
}
