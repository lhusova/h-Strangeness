void prepareFlowCoeff(){

  TFile * flow = new TFile("../data/Flow/v2_pp_HM.root");
  //https://www.hepdata.net/record/ins1784041?version=1
  TFile * spectraPion = new TFile("../data/Flow/spectra_Pion_published.root");
  TFile * spectraKaon = new TFile("../data/Flow/spectra_Kaon_published.root");
  TFile * spectraProton = new TFile("../data/Flow/spectra_Proton_published.root");

  TH1F * pions = (TH1F *) spectraPion->Get("Table 1/Hist1D_y1");
  pions->SetName("pions");
  TH1F * pionsStats = (TH1F *) spectraPion->Get("Table 1/Hist1D_y1_e1");
  TH1F * pionsSyst = (TH1F *) spectraPion->Get("Table 1/Hist1D_y1_e2");

  TH1F * kaons = (TH1F *) spectraKaon->Get("Table 3/Hist1D_y1");
  kaons->SetName("kaons");
  TH1F * kaonsStats = (TH1F *) spectraKaon->Get("Table 3/Hist1D_y1_e1");
  TH1F * kaonsSyst = (TH1F *) spectraKaon->Get("Table 3/Hist1D_y1_e2");

  TH1F * protons = (TH1F *) spectraProton->Get("Table 5/Hist1D_y1");
  protons->SetName("protons");
  TH1F * protonsStats = (TH1F *) spectraProton->Get("Table 5/Hist1D_y1_e1");
  TH1F * protonsSyst = (TH1F *) spectraProton->Get("Table 5/Hist1D_y1_e2");

  TGraphErrors * v2pion = (TGraphErrors *) flow->Get("fFinalGraphPion");
  TGraphErrors * v2pionSyst = (TGraphErrors *) flow->Get("fFinalGraphPionSyst");
  //
  Double_t * pionv2 = v2pion->GetY();
  Double_t * pionv2Stat = v2pion->GetEY();
  Double_t * pionv2Syst = v2pionSyst->GetEY();
  //
  TGraphErrors * v2kaon = (TGraphErrors *) flow->Get("fFinalGraphKaon");
  TGraphErrors * v2kaonSyst = (TGraphErrors *) flow->Get("fFinalGraphKaonSyst");
  //
  Double_t * kaonv2 = v2kaon->GetY();
  Double_t * kaonv2Stat = v2kaon->GetEY();
  Double_t * kaonv2Syst = v2kaonSyst->GetEY();
  //
  TGraphErrors * v2proton = (TGraphErrors *) flow->Get("fFinalGraphProton");
  TGraphErrors * v2protonSyst = (TGraphErrors *) flow->Get("fFinalGraphProtonSyst");
  //
  Double_t * protonv2 = v2proton->GetY();
  Double_t * protonv2Stat = v2proton->GetEY();
  Double_t * protonv2Syst = v2protonSyst->GetEY();

  TGraphErrors * v2K0s = (TGraphErrors *) flow->Get("fFinalGraphK0");
  TGraphErrors * v2K0Syst = (TGraphErrors *) flow->Get("fFinalGraphK0Syst");
  //
  Double_t * k0sv2 = v2K0s->GetY();
  Double_t * k0sv2Stat = v2K0s->GetEY();
  Double_t * k0sv2Syst = v2K0Syst->GetEY();

  TGraphErrors * v2Lambda = (TGraphErrors *) flow->Get("fFinalGraphLam");
  TGraphErrors * v2LambdaSyst = (TGraphErrors *) flow->Get("fFinalGraphLamSyst");
  //
  Double_t * lambdav2 = v2Lambda->GetY();
  Double_t * lambdav2Stat = v2Lambda->GetEY();
  Double_t * lambdav2Syst = v2LambdaSyst->GetEY();

  Double_t binsV2Pion[]=      {0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0, 2.25,2.5,3.0,3.5,4.0,5.0,6.0,10.0};
  Double_t binsV2KaonProton[]={0.2,0.4,0.6,0.8,1.0,1.25,1.5,1.75,2.0,      2.5,3.0,    4.0,    6.0,10.0};

  for (size_t i = 0; i < pions->GetXaxis()->GetNbins(); i++) {
    pions->SetBinError(i+1,TMath::Sqrt(TMath::Power(pionsStats->GetBinContent(i+1),2)+TMath::Power(pionsSyst->GetBinContent(i+1),2)));
    kaons->SetBinError(i+1,TMath::Sqrt(TMath::Power(kaonsStats->GetBinContent(i+1),2)+TMath::Power(kaonsSyst->GetBinContent(i+1),2)));
    protons->SetBinError(i+1,TMath::Sqrt(TMath::Power(protonsStats->GetBinContent(i+1),2)+TMath::Power(protonsSyst->GetBinContent(i+1),2)));
  }
  // pions->DrawCopy();
  TGraphErrors * v2_charged  = (TGraphErrors*) v2kaon->Clone();
  v2_charged->SetName("v2_charged");
  Double_t v2ch;
  Double_t v2ch_err;
  Double_t bin;

  for (size_t i = 0; i < 13; i++) {
    if(i==0) {
      v2ch = pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+1]))*pionv2[i]+kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))*kaonv2[i];
      v2ch = v2ch/(pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+1]))+kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])));
      pions->GetXaxis()->SetRangeUser(binsV2Pion[i],binsV2Pion[i+1]);
      kaons->GetXaxis()->SetRangeUser(binsV2KaonProton[i],binsV2KaonProton[i+1]);
      bin =(pions->GetMean()*pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+1]))+kaons->GetMean()*kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])))/(kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))+pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+1])));
      cout << bin << "   " << pions->GetMean() << endl;
      pions->GetXaxis()->SetRangeUser(0.1,20);
      kaons->GetXaxis()->SetRangeUser(0.2,20);

      // v2ch_err = TMath::Sqrt( TMath::Power(pionv2Stat[i],2)+TMath::Power(pionv2Syst[i],2)*TMath::Power(pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+1])),2)+
      //                         TMath::Power(kaonv2Stat[i],2)+TMath::Power(kaonv2Syst[i],2)*TMath::Power(kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])),2) )/
      //                         (pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+1]))+kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])));
    }
    else if(i==8){
      v2ch = pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+2]))*v2pion->Eval(2.25)+
             kaons->Integral(kaons->FindBin(binsV2Pion[i]),kaons->FindBin(binsV2Pion[i+1]))*kaonv2[i]+
             protons->Integral(protons->FindBin(binsV2Pion[i]),protons->FindBin(binsV2Pion[i+1]))*protonv2[i];
      v2ch = v2ch/(pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+2]))+kaons->Integral(kaons->FindBin(binsV2Pion[i]),kaons->FindBin(binsV2Pion[i+1]))+protons->Integral(protons->FindBin(binsV2Pion[i]),protons->FindBin(binsV2Pion[i+1])));

      pions->GetXaxis()->SetRangeUser(binsV2Pion[i],binsV2Pion[i+2]);
      kaons->GetXaxis()->SetRangeUser(binsV2KaonProton[i],binsV2KaonProton[i+1]);
      protons->GetXaxis()->SetRangeUser(binsV2KaonProton[i],binsV2KaonProton[i+1]);
      bin =(pions->GetMean()*pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+2]))+
            kaons->GetMean()*kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))+
            protons->GetMean()*protons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])))/
            (kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))+pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+2]))+protons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])));
      cout << bin << "   " << pions->GetMean() << endl;
      pions->GetXaxis()->SetRangeUser(0.1,20);
      kaons->GetXaxis()->SetRangeUser(0.2,20);
      protons->GetXaxis()->SetRangeUser(0.3,20);
    }
    else if(i==9){
      v2ch = pions->Integral(pions->FindBin(binsV2Pion[i+1]),pions->FindBin(binsV2Pion[i+2]))*pionv2[i+1]+
             kaons->Integral(kaons->FindBin(binsV2Pion[i]),kaons->FindBin(binsV2Pion[i+1]))*kaonv2[i]+
             protons->Integral(protons->FindBin(binsV2Pion[i]),protons->FindBin(binsV2Pion[i+1]))*protonv2[i];
      v2ch = v2ch/(pions->Integral(pions->FindBin(binsV2Pion[i+1]),pions->FindBin(binsV2Pion[i+2]))+kaons->Integral(kaons->FindBin(binsV2Pion[i]),kaons->FindBin(binsV2Pion[i+1]))+protons->Integral(protons->FindBin(binsV2Pion[i]),protons->FindBin(binsV2Pion[i+1])));

      pions->GetXaxis()->SetRangeUser(binsV2Pion[i+1],binsV2Pion[i+2]);
      kaons->GetXaxis()->SetRangeUser(binsV2KaonProton[i],binsV2KaonProton[i+1]);
      protons->GetXaxis()->SetRangeUser(binsV2KaonProton[i],binsV2KaonProton[i+1]);
      bin =(pions->GetMean()*pions->Integral(pions->FindBin(binsV2Pion[i+1]),pions->FindBin(binsV2Pion[i+2]))+
            kaons->GetMean()*kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))+
            protons->GetMean()*protons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])))/
            (kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))+pions->Integral(pions->FindBin(binsV2Pion[i+1]),pions->FindBin(binsV2Pion[i+2]))+protons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])));
      cout << bin << "   " << pions->GetMean() << endl;
      pions->GetXaxis()->SetRangeUser(0.1,20);
      kaons->GetXaxis()->SetRangeUser(0.2,20);
      protons->GetXaxis()->SetRangeUser(0.3,20);
    }
    else if(i==10){
      v2ch = pions->Integral(pions->FindBin(binsV2Pion[i+1]),pions->FindBin(binsV2Pion[i+3]))*v2pion->Eval(3.5)+
             kaons->Integral(kaons->FindBin(binsV2Pion[i]),kaons->FindBin(binsV2Pion[i+1]))*kaonv2[i]+
             protons->Integral(protons->FindBin(binsV2Pion[i]),protons->FindBin(binsV2Pion[i+1]))*protonv2[i];
      v2ch = v2ch/(pions->Integral(pions->FindBin(binsV2Pion[i+1]),pions->FindBin(binsV2Pion[i+3]))+kaons->Integral(kaons->FindBin(binsV2Pion[i]),kaons->FindBin(binsV2Pion[i+1]))+protons->Integral(protons->FindBin(binsV2Pion[i]),protons->FindBin(binsV2Pion[i+1])));

      pions->GetXaxis()->SetRangeUser(binsV2Pion[i+1],binsV2Pion[i+3]);
      kaons->GetXaxis()->SetRangeUser(binsV2KaonProton[i],binsV2KaonProton[i+1]);
      protons->GetXaxis()->SetRangeUser(binsV2KaonProton[i],binsV2KaonProton[i+1]);
      bin =(pions->GetMean()*pions->Integral(pions->FindBin(binsV2Pion[i+1]),pions->FindBin(binsV2Pion[i+3]))+
            kaons->GetMean()*kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))+
            protons->GetMean()*protons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])))/
            (kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))+pions->Integral(pions->FindBin(binsV2Pion[i+1]),pions->FindBin(binsV2Pion[i+3]))+protons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])));
      cout << bin << "   " << pions->GetMean() << endl;
      pions->GetXaxis()->SetRangeUser(0.1,20);
      kaons->GetXaxis()->SetRangeUser(0.2,20);
      protons->GetXaxis()->SetRangeUser(0.3,20);

    }
    else if(i==11){
      v2ch = pions->Integral(pions->FindBin(binsV2Pion[i+2]),pions->FindBin(binsV2Pion[i+4]))*v2pion->Eval(5)+
             kaons->Integral(kaons->FindBin(binsV2Pion[i]),kaons->FindBin(binsV2Pion[i+1]))*kaonv2[i]+
             protons->Integral(protons->FindBin(binsV2Pion[i]),protons->FindBin(binsV2Pion[i+1]))*protonv2[i];
      v2ch = v2ch/(pions->Integral(pions->FindBin(binsV2Pion[i+2]),pions->FindBin(binsV2Pion[i+4]))+kaons->Integral(kaons->FindBin(binsV2Pion[i]),kaons->FindBin(binsV2Pion[i+1]))+protons->Integral(protons->FindBin(binsV2Pion[i]),protons->FindBin(binsV2Pion[i+1])));

      pions->GetXaxis()->SetRangeUser(binsV2Pion[i+2],binsV2Pion[i+4]);
      kaons->GetXaxis()->SetRangeUser(binsV2KaonProton[i],binsV2KaonProton[i+1]);
      protons->GetXaxis()->SetRangeUser(binsV2KaonProton[i],binsV2KaonProton[i+1]);
      bin =(pions->GetMean()*pions->Integral(pions->FindBin(binsV2Pion[i+2]),pions->FindBin(binsV2Pion[i+4]))+
            kaons->GetMean()*kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))+
            protons->GetMean()*protons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])))/
            (kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))+pions->Integral(pions->FindBin(binsV2Pion[i+2]),pions->FindBin(binsV2Pion[i+4]))+protons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])));
      cout << bin << "   " << pions->GetMean() << endl;
      pions->GetXaxis()->SetRangeUser(0.1,20);
      kaons->GetXaxis()->SetRangeUser(0.2,20);
      protons->GetXaxis()->SetRangeUser(0.3,20);

    }
    else if(i==12){
      v2ch = pions->Integral(pions->FindBin(binsV2Pion[i+3]),pions->FindBin(binsV2Pion[i+4]))*pionv2[i+3]+
             kaons->Integral(kaons->FindBin(binsV2Pion[i]),kaons->FindBin(binsV2Pion[i+1]))*kaonv2[i]+
             protons->Integral(protons->FindBin(binsV2Pion[i]),protons->FindBin(binsV2Pion[i+1]))*protonv2[i];
      v2ch = v2ch/(pions->Integral(pions->FindBin(binsV2Pion[i+3]),pions->FindBin(binsV2Pion[i+4]))+kaons->Integral(kaons->FindBin(binsV2Pion[i]),kaons->FindBin(binsV2Pion[i+1]))+protons->Integral(protons->FindBin(binsV2Pion[i]),protons->FindBin(binsV2Pion[i+1])));

      pions->GetXaxis()->SetRangeUser(binsV2Pion[i+3],binsV2Pion[i+4]);
      kaons->GetXaxis()->SetRangeUser(binsV2KaonProton[i],binsV2KaonProton[i+1]);
      protons->GetXaxis()->SetRangeUser(binsV2KaonProton[i],binsV2KaonProton[i+1]);
      bin =(pions->GetMean()*pions->Integral(pions->FindBin(binsV2Pion[i+3]),pions->FindBin(binsV2Pion[i+4]))+
            kaons->GetMean()*kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))+
            protons->GetMean()*protons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])))/
            (kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))+pions->Integral(pions->FindBin(binsV2Pion[i+3]),pions->FindBin(binsV2Pion[i+4]))+protons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])));
      cout << bin << "   " << pions->GetMean() << endl;
      pions->GetXaxis()->SetRangeUser(0.1,20);
      kaons->GetXaxis()->SetRangeUser(0.2,20);
      protons->GetXaxis()->SetRangeUser(0.3,20);

    }
    else {
      v2ch = pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+1]))*pionv2[i]+
             kaons->Integral(kaons->FindBin(binsV2Pion[i]),kaons->FindBin(binsV2Pion[i+1]))*kaonv2[i]+
             protons->Integral(protons->FindBin(binsV2Pion[i]),protons->FindBin(binsV2Pion[i+1]))*protonv2[i];
      v2ch = v2ch/(pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+1]))+kaons->Integral(kaons->FindBin(binsV2Pion[i]),kaons->FindBin(binsV2Pion[i+1]))+protons->Integral(protons->FindBin(binsV2Pion[i]),protons->FindBin(binsV2Pion[i+1])));
      pions->GetXaxis()->SetRangeUser(binsV2Pion[i],binsV2Pion[i+1]);
      kaons->GetXaxis()->SetRangeUser(binsV2KaonProton[i],binsV2KaonProton[i+1]);
      protons->GetXaxis()->SetRangeUser(binsV2KaonProton[i],binsV2KaonProton[i+1]);
      bin =(pions->GetMean()*pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+1]))+
            kaons->GetMean()*kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))+
            protons->GetMean()*protons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])))/
            (kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1]))+pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+1]))+protons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])));
      cout << bin << "   " << pions->GetMean() << endl;
      pions->GetXaxis()->SetRangeUser(0.1,20);
      kaons->GetXaxis()->SetRangeUser(0.2,20);
      protons->GetXaxis()->SetRangeUser(0.3,20);

      // v2ch_err = TMath::Sqrt( TMath::Power(pionv2Stat[i],2)+TMath::Power(pionv2Syst[i],2)*TMath::Power(pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+1])),2)+
      //                         TMath::Power(kaonv2Stat[i],2)+TMath::Power(kaonv2Syst[i],2)*TMath::Power(kaons->Integral(kaons->FindBin(binsV2KaonProton[i]),kaons->FindBin(binsV2KaonProton[i+1])),2) +
      //                         TMath::Power(protonv2Stat[i],2)+TMath::Power(protonv2Syst[i],2)*TMath::Power(protons->Integral(protons->FindBin(binsV2KaonProton[i]),protons->FindBin(binsV2KaonProton[i+1])),2))/
      //                         (pions->Integral(pions->FindBin(binsV2Pion[i]),pions->FindBin(binsV2Pion[i+1]))+kaons->Integral(kaons->FindBin(binsV2Pion[i]),kaons->FindBin(binsV2Pion[i+1]))+protons->Integral(protons->FindBin(binsV2Pion[i]),protons->FindBin(binsV2Pion[i+1])));
    }
    // bin = (binsV2KaonProton[i+1]+binsV2KaonProton[i])/2;
    v2_charged->SetPoint(i,bin,v2ch);
    v2_charged->SetPointError(i,0.01,TMath::Sqrt(TMath::Power(protonv2Stat[i],2)+TMath::Power(protonv2Syst[i],2)));
  }
  v2_charged->Draw("ap");


  TGraphErrors * v2_Lambda_Combined = (TGraphErrors*) v2Lambda->Clone();
  v2_Lambda_Combined->SetName("v2_Lambda_Combined");

  for (size_t i = 0; i < v2_Lambda_Combined->GetN(); i++) {
    v2_Lambda_Combined->SetPointError(i,0.01,TMath::Sqrt(TMath::Power(lambdav2Stat[i],2)+TMath::Power(lambdav2Syst[i],2)));
  }

  TGraphErrors * v2_K0_Combined = (TGraphErrors*) v2K0s->Clone();
  v2_K0_Combined->SetName("v2_K0_Combined");

  for (size_t i = 0; i < v2_K0_Combined->GetN(); i++) {
    v2_K0_Combined->SetPointError(i,0.01,TMath::Sqrt(TMath::Power(k0sv2Stat[i],2)+TMath::Power(k0sv2Syst[i],2)));
  }

  TFile * newFile = TFile::Open("../data/Flow/prapared_flow.root","RECREATE");
  v2_charged->Write();
  v2_Lambda_Combined->Write();
  v2_K0_Combined->Write();

}
