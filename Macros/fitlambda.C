#include <Riostream.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TAttFill.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>

      Double_t BackgroundFunc(const Double_t* x, const Double_t* par)
    {
    Double_t xVal = x[0];
    if (xVal > 1.10 && xVal < 1.135) {
        TF1::RejectPoint();
        return 0;
    }
    Double_t background = par[0] + par[1] * xVal + par[2]*xVal*xVal;
  //  Double_t background = par[0];
    return background;
    }

void fitlambda(){

    string dir = "/Users/cuikai/Downloads/";
    string name = "AnalysisResults.root";
    string filename = dir + name;
    
    TString histinput[7] = {"h2dMassK0Short","h2dMassLambda","h2dMassAntiLambda","h2dMassXiMinus","h2dMassXiPlus","h2dMassOmegaMinus","h2dMassOmegaPlus"};
    TString histoutput[7] = {"histMassK0Short_","histMassLambda_","histMassAntiLambda_","histMassXiMinus_","histMassXiPlus_","histMassOmegaMinus_","histMassOmegaPlus_"};
    TString meanName[7] = {"K0Short_mean","Lambda_mean","AntiLambda_mean","XiMinus_mean","XiPlus_mean","OmegaMinus_mean","OmegaPlus_mean"};
    TString widthName[7] = {"K0Short_width","Lambda_width","AntiLambda_width","XiMinus_width","XiPlus_width","OmegaMinus_width","OmegaPlus_width"};
    double parameter0[7]={3e+04,-7e+05,1e+05,-3e+04,-3e+04,1e+03,1e+03};
    double parameter1[7]={0,};
    double parameter3[7]={0,};
    double parameter4[7]={0.498,1.115,1.115,1.321,1.321,1.672,1.672};
    double parameter5[7]={0.01,0.003,0.003,0.004,0.004,0.01,0.01};
    double parameterwidth0[7]={0.00354,0.00127,0.00127,0.0016,0.0019,0.00188,0.00185};
    double parameterwidth1[7]={0.000609,0.000172,0.000172,0.000286,0.00022,0.00032,0.00033};
    double parameterwidth2[7]={0.0,0.00261,0.00261,0.00392,0.0037,0.0060,0.0055};
    double parameterwidth3[7]={0.0,2.02,2.02,1.56,1.63,1.77,1.51};
    double parametermean0[7]={0.495,1.114,1.114,1.32095,1.32095,1.67197,1.67177};
    double parametermean1[7]={0.000250,0.000314,0.000314,0.00006717,0.00006717,0.00,0.0000586867};
    double parametermean2[7]={0.00,0.140,0.140,-0.00169506,-0.00169506,-0.01296,-0.0135192};
    double parametermean3[7]={0.00,11.9,11.9,1.28465,1.28465,1.89985,2.04057};
    double bin1[7]={0.4,1.07,1.07,1.25,1.25,1.6,1.6};
    double bin2[7]={0.45,1.10,1.101,1.30,1.30,1.651,1.651};
    double bin3[7]={0.535,1.135,1.135,1.34,1.34,1.69,1.69};
    double bin4[7]={0.6,1.25,1.25,1.42,1.77,1.77};
    double peakPosition[7][20]={0,};
    double peakWidth[7][20]={0,};
    double uncertaintyMean[7][20]={0,};
    double uncertaintyWidth[7][20]={0,};
    double par0,par1,par2=0;
    //get the inputhists from the .root file.
    TFile *f = new TFile( filename.c_str() );
    TDirectoryFile *my_dir = (TDirectoryFile*)f->Get("correlate-strangeness");

    //set some fit cunction: liner to background and liner + gaus to signal
    TF1 *fullfit = new TF1( "fullfit", "[0]+[1]*x+[2]*x*x+gaus(3)", 1.10, 1.135);
    TF1 *bgfit = new TF1( "bgfit", "[0]+[1]*x", 1.07, 1.25);
    
    //create a .root file to save the outputhists.
    TFile *file_Tuple = new TFile ("tuple.root","recreate");
    for (int i=0; i<7; i++){

        Double_t range[2]={bin2[i],bin3[i]};
        TH2F *inputhist = (TH2F*)my_dir->Get( histinput[i] );
        TCanvas *c = new TCanvas();
              TF1 *bgleft = new TF1( "bgleft", BackgroundFunc, 1.07, 1.25, 3);
              bgleft->SetParameter(0, -1.65e7);
              bgleft->SetParameter(1, 8.93e7);
              bgleft->SetParameter(2, -1.02e8);
  

        for (int j=0; j<10; j++){
            TString histnum=histoutput[i]+(j+1);
            
            TH1D *histout = new TH1D(histnum, histnum, 100,5000,5400);
            //slice the 2D hists to 1D hist.
            int binleft = inputhist->GetXaxis()->FindBin(j+1e-6);
            int binright = inputhist->GetXaxis()->FindBin(j+1-1e-6);
            inputhist->ProjectionY(histnum,binleft,binright);
            histout->GetYaxis()->SetTitle("Entry");
            double numEvents = histout->GetEntries();
            if (numEvents < 500) continue;
            //do fit
            histout->Fit("bgleft","","",bin1[i],bin4[i]);
            //set some config parameters.
              par0 = bgleft->GetParameter(0);
              par1 = bgleft->GetParameter(1);
              par2 = bgleft->GetParameter(2);
              bgleft->SetLineColor(kBlue);
              bgleft->Draw("same");
              fullfit->FixParameter(0,par0);
              fullfit->FixParameter(1,par1);
              fullfit->FixParameter(2,par2);
              fullfit->SetParameter(3,parameter3[i]);
              fullfit->SetParameter(4,parameter4[i]);
              fullfit->SetParameter(5,parameter5[i]);
              
            histout->Fit("fullfit","","",bin2[i],bin3[i]);
            fullfit->Draw("same");

            //save the parameters we get from the fit before.
            peakPosition[i][j]=fullfit->GetParameter(4);
            uncertaintyMean[i][j]=fullfit->GetParError(4);
            peakWidth[i][j]= abs(fullfit->GetParameter(5));
            uncertaintyWidth[i][j]=fullfit->GetParError(5);
            TLatex* latex = new TLatex;
            latex->SetNDC();
            latex->SetTextAlign(22);
            latex->SetTextFont(43);
            latex->SetTextSize(20);
            char buffer[50];
            sprintf(buffer, "Mean = %.4f", peakPosition[i][j]);
            latex->DrawLatex(0.3, 0.8, buffer);

            TLatex* latex2 = new TLatex;
            latex2->SetNDC();
            latex2->SetTextAlign(22);
            latex2->SetTextFont(43);
            latex2->SetTextSize(20);
            char buffer2[50];
            sprintf(buffer2, "Width = %.4f", peakWidth[i][j]);
            latex2->DrawLatex(0.3, 0.7, buffer2);           
         
            //cout<<peakPosition[i][j]<<"  "<<peakWidth[i][j]<<endl;
            //histout->Fit("bgfit","","",0.54,0.6);
            //cout<<bg1->GetParameter(0)<<endl;

            //save the hist to the outputfile.
            //histout->Write();
            c->Write(histnum);
            TString outputfilename = TString::Format("./FitPlots/%s.pdf", histnum.Data());
            c->SaveAs(outputfilename);
        }
    }
    //set two fit function.
    TF1 *meanfitK0 = new TF1( "meanfitK0", "[0]+[1]*x", 0, 10);
    TF1 *widthfitK0 = new TF1( "widthfitK0", "[0]+[1]*x", 0, 10);
    TF1 *widthfit = new TF1( "widthfit", "[0]+[1]*x+[2]*TMath::Exp(-[3]*x)", 0, 10);
    TF1 *meanfit = new TF1( "meanfit", "[0]+[1]*x+[2]*TMath::Exp(-[3]*x)", 0, 10);
    TF1 *meanfitnew = new TF1( "meanfitnew", "[0]+[1]*x", 1, 10);
    double x[20],y[20],z[20],errors[20]={0};
    
   for(int k =0; k<7; k++){

        for(int l=0; l<10; l++){
            x[l]=peakPosition[k][l];
            y[l]=peakWidth[k][l];
            errors[l]=uncertaintyWidth[k][l];
            z[l]=l+1;
        }
        
        TCanvas *c = new TCanvas("c","",900,530);
        TH1D * gr2 = new TH1D("gr2","",10,0,10);
        widthfit->SetParameters(parameterwidth0[k], parameterwidth1[k], parameterwidth2[k], parameterwidth3[k]);
        
        for (int num2 = 0; num2 < 10; num2++)
        {
            gr2->SetBinContent(num2+1, y[num2]);
            gr2->SetBinError(num2+1, errors[num2]);
        }
        //gr2->SetMarkerColor(4);
        gr2->SetMarkerStyle(20);
        gr2->SetTitle(widthName[k]);

        //gr2->SetLineColor(0);
        //gr2->RemovePoint(0);
        /*int numPoints = gr2->GetN();
        int tempIndex = 0;
        for (int i = 0; i < numPoints; i++) {
                double x, y;
                gr2->GetPoint(i, x, y);
                if (y == 0) {
                    gr2->RemovePoint(i);
                    i--;  // 由于删除了一个数据点，需要减小索引值
                    numPoints--;
                }
        }*/
        gr2->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        //gr2->GetXaxis()->SetLabelSize(0.04);
        //gr2->GetXaxis()->SetLabelOffset(0.01);
        gr2->GetYaxis()->SetTitle("Width");
        if (k == 0) gr2->Fit("widthfitK0","","",1,10);
        else {
            gr2->Fit("widthfit","","",0,10);
            //gr2->Fit("widthfit1","","",0,10);
        }
        gr2->Draw();
        
        TLatex* latex5 = new TLatex;
        latex5->SetNDC();
        latex5->SetTextAlign(22);
        latex5->SetTextFont(43);
        latex5->SetTextSize(14);
        char buffer5[50];
        sprintf(buffer5, "Angular coefficient = %.6f +/- %.6f", widthfit->GetParameter(1),widthfit->GetParError(1));
        latex5->DrawLatex(0.4, 0.8, buffer5);
        

        TLatex* latex6 = new TLatex;
        latex6->SetNDC();
        latex6->SetTextAlign(22);
        latex6->SetTextFont(43);
        latex6->SetTextSize(14);
        char buffer6[50];
        sprintf(buffer6, "Linear coefficient = %.6f +/- %.6f", widthfit->GetParameter(0),widthfit->GetParError(0));
        latex6->DrawLatex(0.4, 0.7, buffer6);
        c->Write(widthName[k]);

        TString outputfilenamewidth = TString::Format("./FitPlots/%s.pdf", widthName[k].Data());
        c->SaveAs(outputfilenamewidth);  
    }
   for(int k =0; k<7; k++){

        for(int l=0; l<10; l++){
            x[l]=peakPosition[k][l];
            y[l]=peakWidth[k][l];
            errors[l]=uncertaintyMean[k][l];
            z[l]=l+1;
        }
        
        TCanvas *c = new TCanvas("c","",900,530);
       // TGraphErrors *gr = new TGraphErrors(10,z,y);
        TH1D * gr = new TH1D("Mean","",10,0,10);
        meanfit->SetParameters(parametermean0[k], parametermean1[k], parametermean2[k], parametermean3[k]);
        meanfit->SetParLimits(1,0,1);
        for (int num = 0; num < 10; num++)
        {
            gr->SetBinContent(num+1, x[num]);
            gr->SetBinError(num+1, errors[num]);
        }
        //gr2->SetMarkerColor(4);
        gr->SetMarkerStyle(20);
        gr->SetTitle(meanName[k]);
        //gr2->SetLineColor(0);
        //gr->RemovePoint(0);

        gr->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
        gr->GetYaxis()->SetTitle("Mean(GeV/c^{2})");
        if (k == 3 | k == 4)
            {
              gr->SetMinimum(1.318);
              gr->SetMaximum(1.324);
            }
        /*if (k == 5 | k == 6)
            {
              gr->SetMinimum(1.670);
              gr->SetMaximum(1.674);
            }*/
        if(k == 0)  gr->Fit("meanfitK0","","",1,10);
        else {
            gr->Fit("meanfit","","",1,10);
            meanfit->SetLineColor(kBlue);
            meanfit->Draw("same");
            //gr->Fit("meanfitnew","","",1,10);
            //meanfitnew->SetLineColor(kRed);
            //meanfitnew->Draw("same");

        }
        //gr->Draw();
        //TLegend* legend = new TLegend(0.1, 0.7, 0.3, 0.9);
        //legend->AddEntry(meanfit, "Label1", "");
        //legend->AddEntry(graph2, "Label2", "option2");
        // 添加更多的图例项...

        //legend->Draw("same");
        TLatex* latex3 = new TLatex;
        latex3->SetNDC();
        latex3->SetTextAlign(22);
        latex3->SetTextFont(43);
        latex3->SetTextSize(14);
        char buffer3[50];
       // sprintf(buffer3, "Angular coefficient = %.6f +/- %.6f", meanfit->GetParameter(1),meanfit->GetParError(1));
        sprintf(buffer3, "Blue Line : [0]+[1]*x+[2]*TMath::Exp(-[3]*x)");
        latex3->DrawLatex(0.4, 0.8, buffer3);
        

        TLatex* latex4 = new TLatex;
        latex4->SetNDC();
        latex4->SetTextAlign(22);
        latex4->SetTextFont(43);
        latex4->SetTextSize(14);
        char buffer4[50];
        //sprintf(buffer4, "Linear coefficient = %.6f +/- %.6f", meanfit->GetParameter(0),meanfit->GetParError(0));
        sprintf(buffer4, "Red Line : [0]+[1]*x");
        latex4->DrawLatex(0.4, 0.7, buffer4);
        c->Write(meanName[k]);
        TString outputfilenamemean = TString::Format("./FitPlots/%s.pdf", meanName[k].Data());
        c->SaveAs(outputfilenamemean);  
    }
  file_Tuple->Close();
  
  }
