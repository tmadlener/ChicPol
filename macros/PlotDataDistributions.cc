/*
 * PlotDataDistributions.cc
 *
 *  Created on: Feb 14, 2014
 *      Author: valentinknuenz
 */

#include "rootIncludes.inc"
#include "commonVar.h"
#include "RooUtils.h"
#include "THStack.h"
#include "TCutG.h"
#include "TFile.h"
#include "TH1.h"
#include "RooAbsData.h"

using namespace RooFit;

void plotRegionCompDist(RooWorkspace *ws, TTree* tree, int rapBin, int ptBin, int nState,
                        int nbins, double PlotMin, double PlotMax, char varChar[200], bool NormalizeRegions,
                        char xTitle[200], double enlargeYby, double setMinYTimesMaxY, int LegendPosition,
                        char savename[200], bool saveSRhistos);
void plot2Ddist(RooWorkspace *ws, TTree* tree, int rapBin, int ptBin, int nState,
		int nbinsx, int nbinsy, double PlotxMin, double PlotxMax, double PlotyMin, double PlotyMax, char varxChar[200], char varyChar[200], bool Normalize,
		char xTitle[200], char yTitle[200], char savename[200], int region);

//==============================================

void PlotDataDistributions(const std::string &infilename, char infilenameData[200], int rapBin, int ptBin, int nState, int PlottingDataDists){
  RooWorkspace* ws = getFromTFile<RooWorkspace>(infilename, "ws_masslifetime");
  TFile* infileData = new TFile(infilenameData,"R");
  TTree* tree = (TTree*)infileData->Get("selectedData");

  int nbins,nbinsx, nbinsy, LegendPosition;
  double PlotMin, PlotMax, enlargeYby, setMinYTimesMaxY, PlotxMin, PlotxMax, PlotyMin, PlotyMax;
  bool NormalizeRegions, Normalize;
  char savename[200], varChar[200], varxChar[200], varyChar[200], xTitle[200], yTitle[200];
  bool saveSRhistos;
  int region;

  bool plot1D=false;
  bool plot2D=false;


  switch (PlottingDataDists) {
  case 0:
    plot1D=true;
    plot2D=true;
    break;
  case 1:
    plot1D=true;
    break;
  case 2:
    plot2D=true;
    break;
  default:
    std::cerr << "I do not know what do do with this value of Plotting" << std::endl;
  }




  if(plot1D){

    std::cout << ">>>>Plotting chic pT, compare all regions" << std::endl;

    nbins=100; LegendPosition=2;
    PlotMin=0; PlotMax=100;
    //PlotMin=5; PlotMax=30;
    NormalizeRegions=false;
    saveSRhistos=true;
    enlargeYby=1.2; setMinYTimesMaxY=5e-4;
    sprintf(savename,"ChicPt");
    sprintf(varChar,"chic->Pt()");
    sprintf(xTitle,"%s",ws->var("chicPt")->GetTitle());

    plotRegionCompDist(ws, tree, rapBin, ptBin, nState,
                       nbins, PlotMin, PlotMax, varChar, NormalizeRegions,
                       xTitle, enlargeYby, setMinYTimesMaxY, LegendPosition,
                       savename, saveSRhistos);

    std::cout << ">>>>Plotting chic rap, compare all regions" << std::endl;

    nbins=100; LegendPosition=2;
    PlotMin=-2.; PlotMax=2.;
    NormalizeRegions=false;
    saveSRhistos=true;
    enlargeYby=1.2; setMinYTimesMaxY=5e-4;
    sprintf(savename,"ChicRap");
    sprintf(varChar,"chic->Rapidity()");
    sprintf(xTitle,"%s",ws->var("chicRap")->GetTitle());

    plotRegionCompDist(ws, tree, rapBin, ptBin, nState,
                       nbins, PlotMin, PlotMax, varChar, NormalizeRegions,
                       xTitle, enlargeYby, setMinYTimesMaxY, LegendPosition,
                       savename, saveSRhistos);

    std::cout << ">>>>Plotting gamma pT, compare all regions" << std::endl;

    nbins=100; LegendPosition=2;
    PlotMin=0; PlotMax=10;
    NormalizeRegions=false;
    saveSRhistos=true;
    enlargeYby=1.2; setMinYTimesMaxY=5e-4;
    sprintf(savename,"GammaPt");
    sprintf(varChar,"photon->Pt()");
    sprintf(xTitle,"p^{#gamma}_{T}");

    plotRegionCompDist(ws, tree, rapBin, ptBin, nState,
                       nbins, PlotMin, PlotMax, varChar, NormalizeRegions,
                       xTitle, enlargeYby, setMinYTimesMaxY, LegendPosition,
                       savename, saveSRhistos);

    std::cout << ">>>>Plotting gamma rap, compare all regions" << std::endl;

    nbins=100; LegendPosition=2;
    PlotMin=-2.; PlotMax=2.;
    NormalizeRegions=false;
    saveSRhistos=true;
    enlargeYby=1.2; setMinYTimesMaxY=5e-4;
    sprintf(savename,"GammaRap");
    sprintf(varChar,"photon->Rapidity()");
    sprintf(xTitle,"#eta^{#gamma}");

    plotRegionCompDist(ws, tree, rapBin, ptBin, nState,
                       nbins, PlotMin, PlotMax, varChar, NormalizeRegions,
                       xTitle, enlargeYby, setMinYTimesMaxY, LegendPosition,
                       savename, saveSRhistos);

    std::cout << ">>>>Plotting jpsi pt, compare all regions" << std::endl;

    nbins=100; LegendPosition=2;
    PlotMin=9.; PlotMax=60; // NOTE: tmadlener: this are just some "educated guess" min and max values!
    NormalizeRegions=false;
    saveSRhistos=true;
    enlargeYby=1.2; setMinYTimesMaxY=5e-4;
    sprintf(savename,"JpsiPt");
    sprintf(varChar,"jpsi->Pt()");
    sprintf(xTitle,"p^{#psi}_{T}");

    plotRegionCompDist(ws, tree, rapBin, ptBin, nState,
                       nbins, PlotMin, PlotMax, varChar, NormalizeRegions,
                       xTitle, enlargeYby, setMinYTimesMaxY, LegendPosition,
                       savename, saveSRhistos);


    std::cout << ">>>>Plotting jpsi rap, compare all regions" << std::endl;

    nbins=100; LegendPosition=2;
    PlotMin=-2.; PlotMax=2.;
    NormalizeRegions=false;
    saveSRhistos=true;
    enlargeYby=1.2; setMinYTimesMaxY=5e-4;
    sprintf(savename,"JpsiRap");
    sprintf(varChar,"jpsi->Rapidity()");
    sprintf(xTitle,"#eta^{#psi}");

    plotRegionCompDist(ws, tree, rapBin, ptBin, nState,
                       nbins, PlotMin, PlotMax, varChar, NormalizeRegions,
                       xTitle, enlargeYby, setMinYTimesMaxY, LegendPosition,
                       savename, saveSRhistos);


    std::cout << ">>>>Plotting ctau errors, compare all regions" << std::endl;

    nbins=25; LegendPosition=2;
    PlotMin=0; PlotMax=0.06;
    NormalizeRegions=true;
    saveSRhistos=false;
    enlargeYby=1.2; setMinYTimesMaxY=5e-4;
    sprintf(savename,"CtauError_MassRegionComp");
    sprintf(varChar,"JpsictErr");
    sprintf(xTitle,"%s",ws->var("JpsictErr")->GetTitle());

    plotRegionCompDist(ws, tree, rapBin, ptBin, nState,
                       nbins, PlotMin, PlotMax, varChar, NormalizeRegions,
                       xTitle, enlargeYby, setMinYTimesMaxY, LegendPosition,
                       savename, saveSRhistos);

    std::cout << ">>>>Plotting ctau, compare all regions" << std::endl;

    nbins=50; LegendPosition=2;
    PlotMin=-0.3; PlotMax=1.;
    NormalizeRegions=true;
    saveSRhistos=false;
    enlargeYby=1.2; setMinYTimesMaxY=5e-4;
    sprintf(savename,"Ctau_MassRegionComp");
    sprintf(varChar,"Jpsict");
    sprintf(xTitle,"%s",ws->var("Jpsict")->GetTitle());

    plotRegionCompDist(ws, tree, rapBin, ptBin, nState,
                       nbins, PlotMin, PlotMax, varChar, NormalizeRegions,
                       xTitle, enlargeYby, setMinYTimesMaxY, LegendPosition,
                       savename, saveSRhistos);

    std::cout << ">>>>Plotting chic mass, compare all regions" << std::endl;

    nbins=50; LegendPosition=2;
    PlotMin=onia::massChiSBMin; PlotMax=onia::massChiSBMax;
    NormalizeRegions=false;
    saveSRhistos=false;
    enlargeYby=1.2; setMinYTimesMaxY=5e-4;
    sprintf(savename,"ChicMass_MassRegionComp");
    sprintf(varChar,"chic_rf->M()");
    sprintf(xTitle,"%s",ws->var("chicMass")->GetTitle());

    plotRegionCompDist(ws, tree, rapBin, ptBin, nState,
                       nbins, PlotMin, PlotMax, varChar, NormalizeRegions,
                       xTitle, enlargeYby, setMinYTimesMaxY, LegendPosition,
                       savename, saveSRhistos);

    std::cout << ">>>>Plotting probKVF, compare all regions" << std::endl;

    nbins=100; LegendPosition=2;
    PlotMin=0; PlotMax=1.;
    NormalizeRegions=true;
    saveSRhistos=false;
    enlargeYby=1.2; setMinYTimesMaxY=5e-4;
    sprintf(savename,"ProbKVF_MassRegionComp");
    sprintf(varChar,"probKVF");
    sprintf(xTitle,"Kinematical vertex fit probability");

    plotRegionCompDist(ws, tree, rapBin, ptBin, nState,
                       nbins, PlotMin, PlotMax, varChar, NormalizeRegions,
                       xTitle, enlargeYby, setMinYTimesMaxY, LegendPosition,
                       savename, saveSRhistos);


    std::cout << ">>>>Plotting probKVF zoom, compare all regions" << std::endl;

    nbins=100; LegendPosition=2;
    PlotMin=0; PlotMax=0.1;
    NormalizeRegions=true;
    saveSRhistos=false;
    enlargeYby=1.2; setMinYTimesMaxY=5e-4;
    sprintf(savename,"ProbKVF_zoom_MassRegionComp");
    sprintf(varChar,"probKVF");
    sprintf(xTitle,"Kinematical vertex fit probability");

    plotRegionCompDist(ws, tree, rapBin, ptBin, nState,
                       nbins, PlotMin, PlotMax, varChar, NormalizeRegions,
                       xTitle, enlargeYby, setMinYTimesMaxY, LegendPosition,
                       savename, saveSRhistos);

  }

  if(plot2D){

    std::cout << ">>>>Plotting 2D chicPt vs chicRap" << std::endl;

    nbinsx=50; nbinsy=50;
    PlotxMin=-1.3; PlotxMax=1.3;
    PlotyMin=8; PlotyMax=50;
    Normalize=false;
    sprintf(varxChar,"chicRap");
    sprintf(varyChar,"chicPt");
    sprintf(xTitle,"y^{#chi}");
    sprintf(yTitle,"p^{#chi}_{T} [GeV]");

    sprintf(savename,"chicRap_vs_chicPt");
    region=0;

    plot2Ddist(ws, tree, rapBin, ptBin, nState,
               nbinsx, nbinsy, PlotxMin, PlotxMax, PlotyMin, PlotyMax, varxChar, varyChar, Normalize,
               xTitle, yTitle, savename, region);


    nbinsx=50; nbinsy=50;
    PlotxMin=-1.3; PlotxMax=1.3;
    PlotyMin=8; PlotyMax=50;
    Normalize=false;
    sprintf(varxChar,"JpsiRap");
    sprintf(varyChar,"JpsiPt");
    sprintf(xTitle,"y^{#psi}");
    sprintf(yTitle,"p^{#psi}_{T} [GeV]");

    sprintf(savename,"JpsiRap_vs_JpsiPt");
    region=0;

    plot2Ddist(ws, tree, rapBin, ptBin, nState,
               nbinsx, nbinsy, PlotxMin, PlotxMax, PlotyMin, PlotyMax, varxChar, varyChar, Normalize,
               xTitle, yTitle, savename, region);


  }

  delete ws;
}





void plot2Ddist(RooWorkspace *ws, TTree* tree, int rapBin, int ptBin, int nState,
		int nbinsx, int nbinsy, double PlotxMin, double PlotxMax, double PlotyMin, double PlotyMax, char varxChar[200], char varyChar[200], bool Normalize,
		char xTitle[200], char yTitle[200], char savename[200], int region){
  TGaxis::SetMaxDigits(3);




  char varHistChar[200];
  char selChar[200];
  char selRegionChar[200];
  char regionLegend[200];

  std::stringstream binName;
  binName << "data_rap" << rapBin << "_pt" << ptBin<< "_SR";
  RooDataSet *data = (RooDataSet*)ws->data(binName.str().c_str());

  //calculate mass ranges
  double sig1MaxMass = ws->var("var_sig1MaxMass")->getVal();
  double sig1MinMass = ws->var("var_sig1MinMass")->getVal();
  double sig2MaxMass = ws->var("var_sig2MaxMass")->getVal();
  double sig2MinMass = ws->var("var_sig2MinMass")->getVal();
  double lsbMaxMass = ws->var("var_lsbMaxMass")->getVal();
  double lsbMinMass = ws->var("var_lsbMinMass")->getVal();
  double rsbMaxMass = ws->var("var_rsbMaxMass")->getVal();
  double rsbMinMass = ws->var("var_rsbMinMass")->getVal();

  double regionMin[onia::nRegions+1], regionMax[onia::nRegions+1];

  regionMin[0]=lsbMinMass;
  regionMax[0]=rsbMaxMass;
  regionMin[1]=lsbMinMass;
  regionMax[1]=lsbMaxMass;
  regionMin[2]=sig1MinMass;
  regionMax[2]=sig1MaxMass;
  regionMin[3]=sig2MinMass;
  regionMax[3]=sig2MaxMass;
  regionMin[4]=rsbMinMass;
  regionMax[4]=rsbMaxMass;

  sprintf(varHistChar,"%s:%s>>hPlot",varyChar,varxChar);
  sprintf(selRegionChar,"chicMass>%f && chicMass<%f",regionMin[region],regionMax[region]);
  RooDataSet *dataRegion = (RooDataSet*)data->reduce(selRegionChar);

  //TH2F* hPlot = new TH2F("hPlot","hPlot", nbinsx, PlotxMin, PlotxMax, nbinsy, PlotyMin, PlotyMax);

  cout<<"PlotxMin "<<PlotxMin<<endl;
  cout<<"PlotxMax "<<PlotxMax<<endl;
  cout<<"PlotyMin "<<PlotyMin<<endl;
  cout<<"PlotyMax "<<PlotyMax<<endl;

  cout<<"dataRegion->numEntries() "<<dataRegion->numEntries()<<endl;

  TH2F* hPlot = (TH2F*)dataRegion->createHistogram("hPlot",*ws->var(varxChar),Binning(nbinsx, PlotxMin, PlotxMax), YVar(*ws->var(varyChar),Binning(nbinsy, PlotyMin, PlotyMax)));

  if(Normalize) hPlot->Scale(1./hPlot->Integral());
  hPlot->SetTitle("");
  hPlot->SetStats(0);

  hPlot->GetYaxis()->SetTitle(yTitle);
  hPlot->GetYaxis()->SetTitleOffset(1.5);
  hPlot->GetXaxis()->SetTitle(xTitle);
  hPlot->GetXaxis()->SetTitleOffset(1.1);
  hPlot->GetZaxis()->SetTitle("");
  //hPlot->GetZaxis()->SetTitleOffset(1.4);



  double left=0.55, top=0.9, textSize=0.03;
  TLatex *latex=new TLatex();
  latex->SetTextFont(42);
  latex->SetNDC(kTRUE);
  latex->SetTextSize(textSize);
  double step=textSize*1.6;

  gStyle->SetPadBottomMargin(0.08); //0.12
  gStyle->SetPadLeftMargin(0.12); //0.12
  gStyle->SetPadRightMargin(0.15); //0.05
  gStyle->SetPadTopMargin(0.05); //0.05

  TCanvas *c1;


  char ptchar[200];
  char rapchar[200];
  char regionchar[200];

  if(region==0) sprintf(regionchar,"J/#psi SR");
  if(region==1) sprintf(regionchar,"J/#psi SR, #chi LSB");
  if(region==2) sprintf(regionchar,"J/#psi SR, #chi SR1");
  if(region==3) sprintf(regionchar,"J/#psi SR, #chi SR2");
  if(region==4) sprintf(regionchar,"J/#psi SR, #chi RSB");

  if(ptBin==0)
    sprintf(ptchar,"%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTBins[rapBin]]);
  else
    sprintf(ptchar,"%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]);

  if(rapBin==0) sprintf(rapchar,"%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin],onia::KinParticleChar,onia::rapForPTRange[2]);
  else if(rapBin==1) sprintf(rapchar,"|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[rapBin]);
  else sprintf(rapchar,"%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin-1],onia::KinParticleChar,onia::rapForPTRange[rapBin]);


  top=0.9675; left=0.275;




  for(int LinLog=0; LinLog<2; LinLog++){

    c1=new TCanvas("c1","",800,700);

    hPlot->Draw("colz");
    latex->DrawLatex(left,top,Form("%s, %s, %s",regionchar, ptchar, rapchar));

    if(LinLog==0) c1->SetLogz(false);
    if(LinLog==1) c1->SetLogz(true);

    std::stringstream savePlot;
    if(LinLog==0) savePlot << "Figures/PlotDataDistributions/2D_"<<savename<<"_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
    if(LinLog==1) savePlot << "Figures/PlotDataDistributions/2D_"<<savename<<"_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
    c1->SaveAs(savePlot.str().c_str());


  }

  delete c1;

  return;
}







void plotRegionCompDist(
                        RooWorkspace *ws, TTree* tree, int rapBin, int ptBin, int nState,
                        int nbins, double PlotMin, double PlotMax, char varChar[200], bool NormalizeRegions,
                        char xTitle[200], double enlargeYby, double setMinYTimesMaxY, int LegendPosition,
                        char savename[200], bool saveSRhistos
                        ){
  TGaxis::SetMaxDigits(3);

  char yTitle[200];
  if(NormalizeRegions) sprintf(yTitle,"arb. units");
  else sprintf(yTitle,"Frequency");

  double binWidth=(PlotMax-PlotMin)/double(nbins);

  TH1F* hRegionComp[onia::nRegions+1];
  for(int iRegion=0;iRegion<onia::nRegions+1;iRegion++){
    hRegionComp[iRegion] = new TH1F(Form("hRegionComp_region%d",iRegion),Form("hRegionComp_region%d",iRegion),nbins,PlotMin,PlotMax);
  }

  //calculate mass ranges
  double sig1MaxMass = ws->var("var_sig1MaxMass")->getVal();
  double sig1MinMass = ws->var("var_sig1MinMass")->getVal();
  double sig2MaxMass = ws->var("var_sig2MaxMass")->getVal();
  double sig2MinMass = ws->var("var_sig2MinMass")->getVal();
  double lsbMaxMass = ws->var("var_lsbMaxMass")->getVal();
  double lsbMinMass = ws->var("var_lsbMinMass")->getVal();
  double rsbMaxMass = ws->var("var_rsbMaxMass")->getVal();
  double rsbMinMass = ws->var("var_rsbMinMass")->getVal();

  double regionMin[onia::nRegions+1], regionMax[onia::nRegions+1];

  regionMin[0]=lsbMinMass;
  regionMax[0]=rsbMaxMass;
  regionMin[1]=lsbMinMass;
  regionMax[1]=lsbMaxMass;
  regionMin[2]=sig1MinMass;
  regionMax[2]=sig1MaxMass;
  regionMin[3]=sig2MinMass;
  regionMax[3]=sig2MaxMass;
  regionMin[4]=rsbMinMass;
  regionMax[4]=rsbMaxMass;

  char varHistChar[2000];
  char selChar[2000];
  char selRegionChar[2000];
  char regionLegend[2000];

  Double_t ptMin;
  Double_t ptMax;
  Double_t yMin;
  Double_t yMax;

  if(rapBin==0){
    yMin = onia::rapForPTRange[0];
    yMax = onia::rapForPTRange[onia::kNbRapForPTBins];
  } else{
    yMin = onia::rapForPTRange[rapBin-1];
    yMax = onia::rapForPTRange[rapBin];
  }

  if(ptBin==0){
    ptMin = onia::pTRange[rapBin][0];
    ptMax = onia::pTRange[rapBin][onia::kNbPTBins[0]];
  } else{
    ptMin = onia::pTRange[rapBin][ptBin-1];
    ptMax = onia::pTRange[rapBin][ptBin];
  }

  //jpsi signal region
  std::stringstream cutSR;
  cutSR << "TMath::Abs(jpsi->M()-" << ws->var("CBmass_p0_jpsi")->getVal() << ") < " << onia::nSigMass << "*(" << ws->var("CBsigma_p0_jpsi")->getVal() << "+"<< ws->var("CBsigma_p1_jpsi")->getVal() << "*TMath::Abs(jpsi->Rapidity()) + "<< ws->var("CBsigma_p2_jpsi")->getVal() << "*TMath::Abs(jpsi->Rapidity())*TMath::Abs(jpsi->Rapidity()))";
  //cout<<cutSR.str().c_str()<<endl;

  sprintf(selChar,"TMath::Abs(chic->Rapidity())>%f && TMath::Abs(chic->Rapidity())<%f && chic->Pt()>%f && chic->Pt()<%f && %s",yMin,yMax,ptMin,ptMax,cutSR.str().c_str());

  //cout<<selChar<<endl;

  for(int iRegion=0;iRegion<onia::nRegions+1;iRegion++){
    sprintf(varHistChar,"%s>>hRegionComp_region%d",varChar,iRegion);
    sprintf(selRegionChar,"%s && chic->M()>%f && chic->M()<%f",selChar,regionMin[iRegion],regionMax[iRegion]);
    //cout<<selRegionChar<<endl;
    tree->Draw(varHistChar,selRegionChar);
    if(NormalizeRegions) hRegionComp[iRegion]->Scale(1./hRegionComp[iRegion]->Integral());
    hRegionComp[iRegion]->SetTitle("");
    hRegionComp[iRegion]->SetStats(0);
    hRegionComp[iRegion]->SetLineColor(onia::colorRegions[iRegion]);
  }


  double minY = 0.;
  double maxY = 0.;
  maxY=hRegionComp[0]->GetMaximum()*enlargeYby;
  minY=setMinYTimesMaxY*maxY;

  hRegionComp[0]->GetYaxis()->SetTitle(yTitle);
  hRegionComp[0]->GetYaxis()->SetTitleOffset(1.6);
  hRegionComp[0]->GetXaxis()->SetTitle(xTitle);
  hRegionComp[0]->GetXaxis()->SetTitleOffset(1.3);
  hRegionComp[0]->GetYaxis()->SetRangeUser(minY,maxY);

  TLegend* PlotLegend;
  if(LegendPosition==1) PlotLegend = new TLegend(0.15,0.7,0.25,0.925);
  if(LegendPosition==2) PlotLegend = new TLegend(0.8,0.7,0.9,0.925);
  if(LegendPosition==3) PlotLegend = new TLegend(0.8,0.125,0.9,0.35);
  if(LegendPosition==4) PlotLegend = new TLegend(0.15,0.125,0.25,0.35);

  PlotLegend->SetFillColor(kWhite);
  PlotLegend->SetFillStyle(0);
  PlotLegend->SetTextFont(42);
  PlotLegend->SetTextSize(0.035);
  PlotLegend->SetBorderSize(0.);

  for(int iRegion=0;iRegion<onia::nRegions+1;iRegion++){
    if(iRegion==0) sprintf(regionLegend,"All");
    if(iRegion==1) sprintf(regionLegend,"LSB");
    if(iRegion==2) sprintf(regionLegend,"SR1");
    if(iRegion==3) sprintf(regionLegend,"SR2");
    if(iRegion==4) sprintf(regionLegend,"RSB");
    PlotLegend->AddEntry(hRegionComp[iRegion],regionLegend,"l");
  }

  gStyle->SetPadBottomMargin(0.1); //0.12
  gStyle->SetPadLeftMargin(0.12); //0.12
  gStyle->SetPadRightMargin(0.075); //0.05
  gStyle->SetPadTopMargin(0.05); //0.05

  TCanvas *c1;


  for(int LinLog=0; LinLog<2; LinLog++){

    c1=new TCanvas("c1","",1000,900);

    char DrawChar[100];
    sprintf(DrawChar,"");

    hRegionComp[0]->Draw(DrawChar);

    sprintf(DrawChar,"same");
    for(int iRegion=1;iRegion<onia::nRegions+1;iRegion++){
      hRegionComp[iRegion]->Draw(DrawChar);
      if(saveSRhistos&&iRegion>1&&iRegion<4){
        std::stringstream saveHist;
        saveHist << "Figures/PlotDataDistributions/TH1_SR"<<iRegion-1<<"_"<<savename<<"_rap" << rapBin << "_pt" << ptBin << ".root";
        hRegionComp[iRegion]->SaveAs(saveHist.str().c_str());
      }
    }

    PlotLegend->Draw("same");

    if(LinLog==0) c1->SetLogy(false);
    if(LinLog==1) c1->SetLogy(true);

    std::stringstream savePlot;
    if(LinLog==0) savePlot << "Figures/PlotDataDistributions/"<<savename<<"_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
    if(LinLog==1) savePlot << "Figures/PlotDataDistributions/"<<savename<<"_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
    c1->SaveAs(savePlot.str().c_str());


  }

  delete c1;

  return;
}
