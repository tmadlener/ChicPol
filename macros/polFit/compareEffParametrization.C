#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "TStyle.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TH2.h"
#include "TH1.h"
#include "TFrame.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TF1.h"

void compareEffParametrization(){

    gStyle->SetOptStat(kFALSE);
    
    // 2012 binning
    double ptBins[] = {2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 3.75, 4., 4.25, 4.5, 4.75, 5., 5.5, 6., 6.5, 7., 7.5, 8., 8.5, 9., 9.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5, 30., 35., 40., 50., 70.};
    double etaBins[] = {0, 0.2, 0.3, 0.6, 0.8, 1., 1.2, 1.4, 1.6};
    int nPt = sizeof(ptBins)/sizeof(ptBins[0])-1;
    int nEta = sizeof(etaBins)/sizeof(etaBins[0])-1;

    TFile *ftruth = TFile::Open("EffFiles/ParametrizedFactDataEff_2016_01_05_Central.root");
    TFile *fsigmoid = TFile::Open("EffFiles/singleMuonEff_noTracking_L3ptg2_final.root");

    for(int ieta = 0; ieta < nEta; ieta++){

        std::stringstream graph, sigmoid;
        graph << "gEff_DATA_PT_AETA" << ieta;
        sigmoid << "fitTotEff_DATA_pt_etaBin" << ieta;
        TGraph *gdata = (TGraph*)fsigmoid->Get(graph.str().c_str());
        TGraph *gtruth = (TGraph*)ftruth->Get(graph.str().c_str());
        TF1* func = (TF1*)fsigmoid->Get(sigmoid.str().c_str());

        TCanvas *c = new TCanvas();
        TH1F *h = c->DrawFrame(2.,0.,30.,1.);
        h->SetXTitle("p_{T} [GeV]");
        h->SetYTitle("#varepsilon");
        h->Draw();
        stringstream etaRange;
        etaRange << etaBins[ieta] << " < |#eta| < " << etaBins[ieta+1];
        TLegend *l = new TLegend(0.7, 0.15, 0.9, 0.35, etaRange.str().c_str());
        l->SetFillStyle(0);
        l->SetBorderSize(0);
        gdata->SetMarkerStyle(20);
        gdata->Draw("P");
        //gtruth->SetMarkerStyle(1);
        //gtruth->SetMarkerColor(kBlue);
        gtruth->SetLineColor(kBlue);
        gtruth->SetLineWidth(2);
        gtruth->SetLineStyle(1);
        gtruth->Draw("L");
        func->SetLineColor(kRed);
        func->Draw("SAME");

        TLine *line = new TLine();
        line->SetLineStyle(3);
        line->SetLineWidth(2);
        if(ieta < 6)
            line->DrawLine(4.5,0.,4.5,1.);
        else if(ieta == 6)
            line->DrawLine(3.5,0.,3.5,1.);
        else
            line->DrawLine(3.,0.,3.,1.);

        l->AddEntry(gdata, "data", "LP");
        l->AddEntry(gtruth, "truth param.", "L");
        l->AddEntry(func, "sigmoid", "L");
        l->Draw();

        std::stringstream name;
        name << "eff/parametrization_comparison_abseta" << ieta << ".pdf";
        c->SaveAs(name.str().c_str());

    } // ieta

}// void
