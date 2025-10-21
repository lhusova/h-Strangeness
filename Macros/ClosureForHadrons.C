#include "Plotter.h"

void ClosureForHadrons(){

    TFile * fFile = new TFile("../data/AnalysisResults_hhClosure_multEff.root");

    TH3F * fHistRecHadrons = (TH3F *)fFile->Get("h-strange-correlation/hAsssocTrackEtaVsPtVsPhi");
    TH3F * fHistGenHadrons = (TH3F *)fFile->Get("h-strange-correlation__MCgen/ClosureTest/hHadron");
    fHistRecHadrons->GetXaxis()->SetRangeUser(0.251,15);
    fHistGenHadrons->GetXaxis()->SetRangeUser(0.251,15);
    fHistGenHadrons->Sumw2();

    TH1F * fHistRecoPt = (TH1F *) fHistRecHadrons->ProjectionX();
    fHistRecoPt->SetName("fHistRecoPt");
    TH1F * fHistGenPt = (TH1F *)fHistGenHadrons->ProjectionX();
    fHistGenPt -> SetName("fHistGenPt");

    SetHist(fHistGenPt, "", 20, kRed+1, 1.2);
    SetHist(fHistRecoPt, "", 71, kBlue+1, 1.2);

    TCanvas * can = CreateCanvas("can",700,750,true);
    TPad * padRatio = new TPad("padRatio","",0.001,0.001,0.999,0.3);
    padRatio->SetMargin(0.12,0.02,0.25,0.01);
    padRatio->Draw();
    can->GetPadSave()->cd();
    can->GetPadSave()->SetLogy();
    SetHistAxes(fHistGenPt, "p_{T} (GeV/c)", "#");
    fHistGenPt->DrawCopy();
    fHistRecoPt->DrawCopy("same");

    padRatio->cd();
    TH1F * fHistRatio = (TH1F *)fHistRecoPt->Clone();
    fHistRatio->SetName("fHistRatio");
    fHistRatio->Divide(fHistGenPt);
    SetHistAxesSmallPad(fHistRatio, "p_{T} (GeV/c)", "rec/gen")  ;
    fHistRatio->DrawCopy();
    DrawConstant(kBlack,1,0.2,15);

    TH2F* histTriggersRec = (TH2F*) fFile->Get("h-strange-correlation/sameEvent/TriggerParticlesHadron");
    histTriggersRec->GetXaxis()->SetRangeUser(2,15);
    TH3F * histTriggersGen = (TH3F*) fFile->Get("h-strange-correlation__MCgen/ClosureTest/hTrigger");
    histTriggersGen->Sumw2();
    histTriggersGen->GetXaxis()->SetRangeUser(2,15);

    TH1F * fHistRecoPtTrigg = (TH1F *) histTriggersRec->ProjectionX();
    fHistRecoPtTrigg->SetName("fHistRecoPtTrigg");
    TH1F * fHistGenPtTrigg = (TH1F *)histTriggersGen->ProjectionX();
    fHistGenPtTrigg -> SetName("fHistGenPtTrigg");

    SetHist(fHistGenPtTrigg, "", 20, kRed+1, 1.2);
    SetHist(fHistRecoPtTrigg, "", 71, kBlue+1, 1.2);

    TCanvas * canTrigg = CreateCanvas("canTrigg",700,750,true);
    TPad * padRatioTrigg = new TPad("padRatioTrigg","",0.001,0.001,0.999,0.3);
    padRatioTrigg->SetMargin(0.12,0.02,0.25,0.01);
    padRatioTrigg->Draw();
    canTrigg->GetPadSave()->cd();
    canTrigg->GetPadSave()->SetLogy();
    SetHistAxes(fHistGenPtTrigg, "p_{T} (GeV/c)", "#");
    fHistGenPtTrigg->DrawCopy();
    fHistRecoPtTrigg->DrawCopy("same");

    padRatioTrigg->cd();
    TH1F * fHistRatioTrigg = (TH1F *)fHistRecoPtTrigg->Clone();
    fHistRatioTrigg->SetName("fHistRatio");
    fHistRatioTrigg->Divide(fHistGenPtTrigg);
    SetHistAxesSmallPad(fHistRatioTrigg, "p_{T} (GeV/c)", "rec/gen")  ;
    fHistRatioTrigg->DrawCopy();
    DrawConstant(kBlack,1,2,15);
}