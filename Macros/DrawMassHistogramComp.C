#include <Riostream.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TAttFill.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>

void DrawMassHistogramComp(){

  TString histinput1[7] = {"hstrangecorrelationfilter/h3dMassK0Short","hstrangecorrelationfilter/h3dMassLambda","hstrangecorrelationfilter/h3dMassAntiLambda",
                            "hstrangecorrelationfilter/h3dMassXiMinus","hstrangecorrelationfilter/h3dMassXiPlus","hstrangecorrelationfilter/h3dMassOmegaMinus","hstrangecorrelationfilter/h3dMassOmegaPlus"};
  TString histinput2[7] = {"correlate-strangeness/h2dMassK0Short","correlate-strangeness/h2dMassLambda","correlate-strangeness/h2dMassAntiLambda",
                            "correlate-strangeness/h2dMassXiMinus","correlate-strangeness/h2dMassXiPlus","correlate-strangeness/h2dMassOmegaMinus","correlate-strangeness/h2dMassOmegaPlus"};

  TString histoutput2[7] = {"histMassK0ShortOld_","histMassLambdaOld_","histMassAntiLambdaOld_","histMassXiMinusOld_","histMassXiPlusOld_","histMassOmegaMinusOld_","histMassOmegaPlusOld_"};
  TString histoutput1[7] = {"histMassK0Short_","histMassLambda_","histMassAntiLambda_","histMassXiMinus_","histMassXiPlus_","histMassOmegaMinus_","histMassOmegaPlus_"};

  double parameter0[7]={3e+04,-7e+05,1e+05,-3e+04,-3e+04,1e+03,1e+03};
  double parameter1[7]={0,};
  double parameter2[7]={0,};
  double parameter3[7]={0.498,1.115,1.115,1.321,1.321,1.672,1.672};
  double parameter4[7]={0.01,0.003,0.003,0.004,0.004,0.01,0.01};
  double parameterwidth0[7]={0.0038897,0.00120573,0.00120573,0.00160825,0.00160825,0.0018886,0.0018886};
  double parameterwidth1[7]={0.00057291,0.000232949,0.000232949,0.000286217,0.000286217,0.000324539,0.000324539};
  double parameterwidth2[7]={0.0,0.00306607,0.00306607,0.00392055,0.00392055,0.00606385,0.00606385};
  double parameterwidth3[7]={0.0,2.02,2.02,1.56117,1.56117,1.77199,1.77199};
  double parametermean0[7]={0.495837,1.11548,1.11548,1.32202,1.32202,1.67177,1.67177};
  double parametermean1[7]={0.000202455,0.0000331134,1.33227e-15,1.33227e-15,0.00006717,0.0000586899,0.0000586899};
  double parametermean2[7]={0.00,-0.00254089,-0.00254089,-0.00203195,-0.00203195,-0.00135193,-0.00135193};
  double parametermean3[7]={0.00,0.96783,0.96783,0.255695,0.255695,2.04058,2.04058};
  TFile *file1 = new TFile("AnalysisResults.root", "READ");
  TFile *file2 = new TFile("AnalysisResultsOld.root", "READ");
  //TDirectoryFile *my_dir = (TDirectoryFile*)file->Get("hstrangecorrelationfilter");
  TFile *file_Tuple = new TFile ("draw.root","recreate");

      for (int i=0; i<7; i++){
        
       TCanvas *c1 = new TCanvas("c1", "", 800,600);

        TH3D *inputhist1 = (TH3D*) file1->Get(histinput1[i]);
        TH2D *inputhist2 = (TH2D*) file2->Get(histinput2[i]);
        c1->SetTicks(1,1);
        c1->SetLeftMargin(0.163);
        c1->SetBottomMargin(0.16);
        c1->SetRightMargin(0.147);
        c1->SetTopMargin(0.05);
        c1->SetFrameFillStyle(0);
        c1->SetFillStyle(0);
        gStyle->SetOptStat(0);


        for (int j=0; j<10; j++){

            TString histnum1=histoutput1[i]+(j+1);
            TString histnum2=histoutput2[i]+(j+1);
            
            TH1D *histout1 = new TH1D(histnum1, histnum1, 100,5000,5400);
            TH1D *histout2 = new TH1D(histnum2, histnum2, 100,5000,5400);
            //slice the 2D hists to 1D hist.
            int binleft = inputhist1->GetXaxis()->FindBin(j+1e-6);
            int binright = inputhist1->GetXaxis()->FindBin(j+1-1e-6);
            inputhist1->ProjectionY(histnum1,binleft,binright,1,inputhist1->GetZaxis()->GetNbins());
            inputhist2->ProjectionY(histnum2,binleft,binright);
            //histout1->GetYaxis()->SetTitle("Entry");
            double numEvents = histout1->GetEntries();
            if (numEvents < 500) continue;
            histout1->SetTitle(histnum1);
            histout1->SetLineColor(kBlack);
            histout1->SetLineWidth(2);
            double integral1 = histout1->Integral();
            double targetValue = 1.0;
            histout1->Scale(targetValue / integral1);
            histout1->Draw("hist");
            double integral2 = histout2->Integral();
            histout2->Scale(targetValue / integral2);
           // histout2->Scale(0.1);
            histout2->SetTitle(histnum1);
            histout2->SetLineColor(kRed);
            histout2->SetLineWidth(2);
            histout2->Draw("hist same");

            TLegend *legend = new TLegend(0.2, 0.7, 0.4, 0.9);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);

            // 将曲线添加到图例中，并设置标签
            legend->AddEntry(histout1, "New Train", "l");
            legend->AddEntry(histout2, "Previous Train", "l");
            legend->SetTextSize(0.04); 
            legend->Draw("same");

            /*
            TH1D *hBg = new TH1D("hBg", "", 10000, histout1->GetXaxis()->GetXmin(), histout1->GetXaxis()->GetXmax());
  
            for(Int_t ii=1; ii<hBg->GetNbinsX()+1; ii++){
                double x = hBg->GetBinCenter(ii);
    
                double mean = parametermean0[i]+parametermean1[i]*(j+0.5)+parametermean2[i]*TMath::Exp(-parametermean3[i]*(j+0.5));
                double width = parameterwidth0[i]+parameterwidth1[i]*(j+0.5)+parameterwidth2[i]*TMath::Exp(-parameterwidth3[i]*(j+0.5));
            bool isThisBgRegion = TMath::Abs(x-mean)>3*width && TMath::Abs(x-mean)<6*width;
    
            if( isThisBgRegion )
              hBg->SetBinContent(ii, histout1->GetBinContent( histout1->FindBin(x)));
            }

            hBg->SetFillColorAlpha(2,0.5);
            histout1->Draw("same");
            hBg->Draw("same") ;*/
            c1->Write(histnum1);
            TString outputfilename = TString::Format("./NewPlots/%s.pdf", histnum1.Data());
            c1->SaveAs(outputfilename);
        }
      }
     file_Tuple->Close();
 }
