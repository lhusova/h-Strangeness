std::vector<double> GetAxisEdges(TH1C *hist){
  std::vector<double> returnArray;
  //cout<<"bins: "<<hist->GetNbinsX()<<endl;
  for(Int_t ib=1; ib<hist->GetNbinsX()+2; ib++){
    //cout<<"edges: "<<hist->GetBinLowEdge(ib)<<endl;
    returnArray.emplace_back(hist->GetBinLowEdge(ib));
  }
  return returnArray;
}

THnF* GetTHnF(TFile *file, TString newHistoName = "newHistogram", TString histoName = "correlate-strangeness/sameEvent/Signal/K0Short"){
  // returns THnF with expanded axes
  // step 1: grab the necessary 1D axes from the reference TH1Cs
  TH1C *hDPhi = (TH1C*) file->Get("correlate-strangeness/axes/hDeltaPhiAxis");
  TH1C *hDEta = (TH1C*) file->Get("correlate-strangeness/axes/hDeltaEtaAxis");
  TH1C *hPtAssoc = (TH1C*) file->Get("correlate-strangeness/axes/hPtAssocAxis");
  TH1C *hVtxZ = (TH1C*) file->Get("correlate-strangeness/axes/hVertexZAxis");
  TH1C *hMult = (TH1C*) file->Get("correlate-strangeness/axes/hMultAxis");
  TH1C *hPtTrigg = (TH1C*) file->Get("correlate-strangeness/axes/hPtTriggerAxis");

  std::vector<std::vector<double>> expandedAxes;
  std::vector<double> axisDPhi = GetAxisEdges(hDPhi);
  std::vector<double> axisDEta = GetAxisEdges(hDEta);
  std::vector<double> axisPtAssoc = GetAxisEdges(hPtAssoc);
  std::vector<double> axisVtxZ = GetAxisEdges(hVtxZ);
  std::vector<double> axisMult = GetAxisEdges(hMult);
  std::vector<double> axisPtTrigg = GetAxisEdges(hPtTrigg);

  expandedAxes.emplace_back(axisDPhi);
  expandedAxes.emplace_back(axisDEta);
  expandedAxes.emplace_back(axisPtAssoc);
  expandedAxes.emplace_back(axisPtTrigg);
  expandedAxes.emplace_back(axisVtxZ);
  expandedAxes.emplace_back(axisMult);

  Int_t nbins[6];
  nbins[0] = hDPhi->GetNbinsX();
  nbins[1] = hDEta->GetNbinsX();
  nbins[2] = hPtAssoc->GetNbinsX();
  nbins[3] = hPtTrigg->GetNbinsX();
  nbins[4] = hVtxZ->GetNbinsX();
  nbins[5] = hMult->GetNbinsX();

  //  The constructor we want to use
  //  THn(const char *name, const char *title, Int_t dim, const Int_t *nbins,
  //      const std::vector<std::vector<double>> &xbins);
  THnF *hReturnHisto = new THnF(newHistoName.Data(), "", 6, nbins, expandedAxes);

  // get the information to populate
  cout<<"Loading and expanding histogram: "<<histoName.Data()<<endl;
  THnF *hNd = (THnF*) file->Get(histoName.Data());

  Int_t coordinate[6];
  Int_t coordinateNew[6];
  double content;

  for(Int_t i1=0; i1<nbins[0]; i1++){
    for(Int_t i2=0; i2<nbins[1]; i2++){
      for(Int_t i3=0; i3<nbins[2]; i3++){
        for(Int_t i4=0; i4<nbins[3]; i4++){
          for(Int_t i5=0; i5<nbins[4]; i5++){
            for(Int_t i6=0; i6<nbins[5]; i6++){
              coordinate[0] = i1; coordinate[1] = i2; coordinate[2] = i3; coordinate[3] = i4; coordinate[4] = i5; coordinate[5] = i6;
              coordinateNew[0] = i1+1; coordinateNew[1] = i2+1; coordinateNew[2] = i3+1; coordinateNew[3] = i4+1; coordinateNew[4] = i5+1; coordinateNew[5] = i6+1;
              content = hNd->GetBinContent(coordinate);
              hReturnHisto->SetBinContent(coordinateNew, content);
            }
          }
        }
      }
    }
  }
  return hReturnHisto;
}
