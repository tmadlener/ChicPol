/*
 * PlotJpsiMassLifetime.cc
 *
 *  Created on: Apr 23, 2014
 *      Author: valentinknuenz
 */



#include "rootIncludes.inc"
#include "commonVar.h"
#include "RooUtils.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooFormulaVar.h"

using namespace RooFit;

double plotMass(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom);
double plotMassRap(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom);
double plotLifeSig(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom);
double plotLifeBg(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom, int SB);
double plotJpsiPedagogical(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom, int LinLog);
void latexFloatingLifetimePars(RooWorkspace *ws, TLatex* latex, int region);
//==============================================

void PlotMassLifetime(const std::string &infilename, int rapBin, int ptBin, int nState, int PlottingJpsi){

  TFile *infile = new TFile(infilename.c_str(), "UPDATE");
  if(!infile){
    std::cout << "Error: failed to open file with dataset" << std::endl;
  }
  RooWorkspace *ws=(RooWorkspace *)infile->Get("ws_masslifetime");
  if(!ws){
    std::cout << "Error: failed to open workspace " << std::endl;
  }

  nState=4;

  bool SpeedPlotting=true;
  int nSpeedPlotting=250;
  bool plotZoom=false;

  //bool correctResolutionForPlotting=true;
  //double resCorrFactor=1.0525;

  double jpsi_chi2ndfMass=0.;
  double jpsi_chi2ndfLifeSR=0.;
  double jpsi_chi2ndfLifeLSB=0.;
  double jpsi_chi2ndfLifeRSB=0.;
  double jpsi_chi2ndfbuffer=0.;

  switch (PlottingJpsi) {
  case 1:
    std::cout << ">>>>PlottingJpsi mass" << std::endl;
    jpsi_chi2ndfMass = plotMass(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
    std::cout << ">>>>PlottingJpsi lifetime sidebands" << std::endl;
    jpsi_chi2ndfLifeLSB = plotLifeBg(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false, 1);
    jpsi_chi2ndfbuffer = plotLifeBg(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true, 1);
    jpsi_chi2ndfLifeRSB = plotLifeBg(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false, 2);
    jpsi_chi2ndfbuffer = plotLifeBg(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true, 2);
    std::cout << ">>>>PlottingJpsi lifetime signal region" << std::endl;
    jpsi_chi2ndfLifeSR = plotLifeSig(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
    jpsi_chi2ndfbuffer = plotLifeSig(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
    std::cout << ">>>>PlottingJpsi mass-rap" << std::endl;
    jpsi_chi2ndfbuffer = plotMassRap(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
    std::cout << ">>>>PlottingJpsi pedagogical lin" << std::endl;
    jpsi_chi2ndfbuffer = plotJpsiPedagogical(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true, 1);
    std::cout << ">>>>PlottingJpsi pedagogical log" << std::endl;
    jpsi_chi2ndfbuffer = plotJpsiPedagogical(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true, 0);
    break;
  case 2:
    std::cout << ">>>>PlottingJpsi mass" << std::endl;
    jpsi_chi2ndfMass = plotMass(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
    break;
  case 3:
    std::cout << ">>>>PlottingJpsi lifetime sidebands" << std::endl;
    jpsi_chi2ndfLifeLSB = plotLifeBg(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false, 1);
    jpsi_chi2ndfbuffer = plotLifeBg(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true, 1);
    jpsi_chi2ndfLifeRSB = plotLifeBg(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false, 2);
    jpsi_chi2ndfbuffer = plotLifeBg(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true, 2);
    break;
  case 4:
    std::cout << ">>>>PlottingJpsi lifetime signal region" << std::endl;
    jpsi_chi2ndfLifeSR = plotLifeSig(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
    jpsi_chi2ndfbuffer = plotLifeSig(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
    break;
  case 5:
    std::cout << ">>>>PlottingJpsi mass-rap" << std::endl;
    jpsi_chi2ndfbuffer = plotMassRap(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
    break;
  case 6:
    std::cout << ">>>>PlottingJpsi lifetime sidebands" << std::endl;
    jpsi_chi2ndfLifeLSB = plotLifeBg(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false, 1);
    jpsi_chi2ndfbuffer = plotLifeBg(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true, 1);
    jpsi_chi2ndfLifeRSB = plotLifeBg(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false, 2);
    jpsi_chi2ndfbuffer = plotLifeBg(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true, 2);
    std::cout << ">>>>PlottingJpsi lifetime signal region" << std::endl;
    jpsi_chi2ndfLifeSR = plotLifeSig(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, false);
    jpsi_chi2ndfbuffer = plotLifeSig(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true);
    break;
  case 7:
    std::cout << ">>>>PlottingJpsi pedagogical lin" << std::endl;
    jpsi_chi2ndfbuffer = plotJpsiPedagogical(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true, 1);
    std::cout << ">>>>PlottingJpsi pedagogical log" << std::endl;
    jpsi_chi2ndfbuffer = plotJpsiPedagogical(ws, rapBin, ptBin, nState, SpeedPlotting, nSpeedPlotting, true, 0);
    break;
  default:
    std::cerr << "I do not know what do do with this value of PlottingJpsi" << std::endl;
  }
  infile->Close();

  cout<<"writing chi2s to ws"<<endl;

  TFile *infileOut = new TFile(infilename.c_str(), "UPDATE");
  if(!infile){
    std::cout << "Error: failed to open file with dataset" << std::endl;
  }
  RooWorkspace *wsOut=(RooWorkspace *)infileOut->Get("ws_masslifetime");
  if(!ws){
    std::cout << "Error: failed to open workspace " << std::endl;
  }

  if(PlottingJpsi==1 || PlottingJpsi==2) {RooRealVar var_jpsi_chi2ndf_mass("var_jpsi_chi2ndf_mass","var_jpsi_chi2ndf_mass",jpsi_chi2ndfMass); var_jpsi_chi2ndf_mass.setVal(jpsi_chi2ndfMass); if(!wsOut->var("var_jpsi_chi2ndf_mass")) wsOut->import(var_jpsi_chi2ndf_mass); else wsOut->var("var_jpsi_chi2ndf_mass")->setVal(jpsi_chi2ndfMass);}
  if(PlottingJpsi==1 || PlottingJpsi==4 || PlottingJpsi==6) {RooRealVar var_chi2ndf_Lifetime_SR("var_chi2ndf_Lifetime_SR","var_chi2ndf_Lifetime_SR",jpsi_chi2ndfLifeSR); var_chi2ndf_Lifetime_SR.setVal(jpsi_chi2ndfLifeSR); if(!wsOut->var("var_chi2ndf_Lifetime_SR")) wsOut->import(var_chi2ndf_Lifetime_SR); else wsOut->var("var_chi2ndf_Lifetime_SR")->setVal(jpsi_chi2ndfLifeSR);}
  if(PlottingJpsi==1 || PlottingJpsi==3 || PlottingJpsi==6) {RooRealVar var_chi2ndf_Lifetime_LSB("var_chi2ndf_Lifetime_LSB","var_chi2ndf_Lifetime_LSB",jpsi_chi2ndfLifeLSB); var_chi2ndf_Lifetime_LSB.setVal(jpsi_chi2ndfLifeLSB); if(!wsOut->var("var_chi2ndf_Lifetime_LSB")) wsOut->import(var_chi2ndf_Lifetime_LSB); else wsOut->var("var_chi2ndf_Lifetime_LSB")->setVal(jpsi_chi2ndfLifeLSB);}
  if(PlottingJpsi==1 || PlottingJpsi==3 || PlottingJpsi==6) {RooRealVar var_chi2ndf_Lifetime_RSB("var_chi2ndf_Lifetime_RSB","var_chi2ndf_Lifetime_RSB",jpsi_chi2ndfLifeRSB);	var_chi2ndf_Lifetime_RSB.setVal(jpsi_chi2ndfLifeRSB); if(!wsOut->var("var_chi2ndf_Lifetime_RSB")) wsOut->import(var_chi2ndf_Lifetime_RSB); else wsOut->var("var_chi2ndf_Lifetime_RSB")->setVal(jpsi_chi2ndfLifeRSB);}



  std::cout<< "write wsOut" <<std::endl;
  wsOut->Write();
  std::cout<< "print wsOut" <<std::endl;
  wsOut->Print("v");

  infileOut->Close();


  delete ws;
}

//==============================================
double plotMass(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom){
  int  nbins=90; //0.005 bin size
  TGaxis::SetMaxDigits(3);

  //nbins=32;

  gSystem->mkdir("Fit/jpsiFit",kTRUE);

  double binWidth=(onia::massMax-onia::massMin)/double(nbins)*1000;

  RooRealVar *JpsiMass = ws->var("JpsiMass");
  assert( 0 != JpsiMass );
  RooRealVar *JpsiRap = ws->var("JpsiRap");
  assert( 0 != JpsiRap );

  RooPlot *massFrame = JpsiMass->frame(Bins(nbins));
  assert ( 0 != massFrame );
  massFrame->SetName(Form("mass_plot_rap%d_pt%d",rapBin,ptBin));
  massFrame->SetTitle("");
  massFrame->GetYaxis()->SetTitle(Form("Events / %1.0f MeV",binWidth));
  massFrame->GetYaxis()->SetTitleOffset(1.3);

  RooPlot *massFramePull = JpsiMass->frame(Bins(nbins));
  assert ( 0 != massFramePull );
  massFramePull->SetName(Form("pullmass_plot_rap%d_pt%d",rapBin,ptBin));
  massFramePull->SetTitle("");
  massFramePull->GetYaxis()->SetTitle("pull");
  massFramePull->GetXaxis()->SetTitleSize(0.08);
  massFramePull->GetYaxis()->SetTitleSize(0.08);
  massFramePull->GetXaxis()->SetLabelSize(0.08);
  massFramePull->GetYaxis()->SetLabelSize(0.08);
  massFramePull->GetYaxis()->SetTitleOffset(0.4);
  massFramePull->GetYaxis()->SetRangeUser(-5.99,5.99);

  RooAbsData *data= ws->data(Form("jpsi_data_rap%d_pt%d",rapBin,ptBin));
  assert ( 0 != data );
  RooFitResult* fitRlt = dynamic_cast<RooFitResult*>(ws->obj(Form("m_fitresult_rap%d_pt%d",rapBin,ptBin)));
  assert ( 0 != fitRlt);

  double Sigma = -1.0, SigmaErr = -1.0;

  if(ws->var("CBsigma_jpsi")!=NULL){
    Sigma=ws->var("CBsigma_jpsi")->getVal();
    SigmaErr=ws->var("CBsigma_jpsi")->getError();
  }
  else{
    Sigma=ws->function("CBsigma_jpsi")->getVal();
    SigmaErr=0;
  }

  double Meanp0 = -1.0, Meanp0Err = -1.0;
  getVarFromWorkspace(ws, "CBmass_p0_jpsi", Meanp0, Meanp0Err);
  double Meanp1 = -1.0, Meanp1Err = -1.0;
  getVarFromWorkspace(ws, "CBmass_p1_jpsi", Meanp1, Meanp1Err);
  double Meanp2 = -1.0, Meanp2Err = -1.0;
  getVarFromWorkspace(ws, "CBmass_p2_jpsi", Meanp2, Meanp2Err);

  double Mean=Meanp0;

  double Sigmap0 = -1.0, Sigmap0Err = -1.0;
  getVarFromWorkspace(ws, "CBsigma_p0_jpsi", Sigmap0, Sigmap0Err);
  double Sigmap1 = -1.0, Sigmap1Err = -1.0;
  getVarFromWorkspace(ws, "CBsigma_p1_jpsi", Sigmap1, Sigmap1Err);
  double Sigmap2 = -1.0, Sigmap2Err = -1.0;
  getVarFromWorkspace(ws, "CBsigma_p2_jpsi", Sigmap2, Sigmap2Err);

  double Alphap0 = -1.0, Alphap0Err = -1.0;
  getVarFromWorkspace(ws, "CBalpha_p0_jpsi", Alphap0, Alphap0Err);
  double Alphap1 = -1.0, Alphap1Err = -1.0;
  getVarFromWorkspace(ws, "CBalpha_p1_jpsi", Alphap1, Alphap1Err);

  double cbN = -1.0, cbNErr = -1.0;
  getVarFromWorkspace(ws, "CBn_jpsi", cbN, cbNErr);
  double lambda = -1.0, lambdaErr = -1.0;
  getVarFromWorkspace(ws, "bkgLambda_jpsi", lambda, lambdaErr);
  double fracCB1 = -1.0, fracCB1Err = -1.0;
  getVarFromWorkspace(ws, "fracCB1_jpsi", fracCB1, fracCB1Err);
  double fracBkg = -1.0, fracBkgErr = -1.0;
  getVarFromWorkspace(ws, "fracBkg_jpsi", fracBkg, fracBkgErr);


  std::stringstream masssnapshotname;
  masssnapshotname << "m_snapshot_rap" << rapBin << "_pt" << ptBin;
  ws->loadSnapshot(masssnapshotname.str().c_str());

  RooAbsPdf *massPdf = ws->pdf("massModel_jpsi");
  assert ( 0 != massPdf );
  RooAbsPdf *bkgMassShape = ws->pdf("bkgMassShape_jpsi");
  assert ( 0 != bkgMassShape );
  RooAbsPdf *sigMassShape_jpsi = ws->pdf("sigMassShape_jpsi");
  assert ( 0 != sigMassShape_jpsi );
  RooAbsPdf *gaussMassShape_jpsi = ws->pdf("gaussMassShape_jpsi");
  assert ( 0 != gaussMassShape_jpsi );
  RooAbsPdf *FullmassPdf = ws->pdf("FullmassPdf");
  assert ( 0 != FullmassPdf );

  int nEntries = data->numEntries();

  int JpsiFitColor=633;

  cout<<"plot data"<<endl;
  data->plotOn(massFrame,MarkerSize(onia::markerSize_ML), Name("myHist"));
  cout<<"plot FullmassPdf"<<endl;
  FullmassPdf->plotOn(massFrame,
                      Normalization(nEntries,2),
                      LineWidth(2),
                      LineColor(JpsiFitColor), Name("myCurve")
                      );

  //------get chi2------------SHOULD be DONE after PLOTTING------
  int parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.
  int nBins_Mass=massFrame->GetNbinsX();
  double chi2Pre_Mass=massFrame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
  int ndof_Mass=nBins_Mass-parsFit;  //num of degree of freedom
  double chi2_Mass=chi2Pre_Mass*ndof_Mass;

  RooHist* hpull_mass = massFrame->pullHist("myHist","myCurve",kTRUE);

  hpull_mass->SetMarkerSize(0.8);
  for(int i=0;i<hpull_mass->GetN();i++){
    hpull_mass->SetPointEYlow(i,0.);
    hpull_mass->SetPointEYhigh(i,0.);
  }
  massFramePull->addPlotable(hpull_mass,"P");

  cout<<"plot sigMassShape_jpsi"<<endl;
  sigMassShape_jpsi->plotOn(massFrame,
                            Normalization(nEntries*(1.-fracBkg),2),
                            LineStyle(2),
                            LineColor(onia::ColorSumJpsiSignal),
                            LineWidth(2),
                            ProjWData(*JpsiRap, *data));
  cout<<"plot bkgMassShape"<<endl;
  bkgMassShape->plotOn(massFrame,
                       Normalization(nEntries*fracBkg,2),
                       LineStyle(2),
                       LineColor(onia::ColorMuMuBG),
                       LineWidth(2),
                       ProjWData(*JpsiRap, *data));

  //cout<<"plot data"<<endl;
  //data->plotOn(massFrame,MarkerSize(onia::markerSize_ML));


  double minY = 0.;

  double maxY = 0.;
  if(nState == 4) maxY = massFrame->GetMaximum()*0.3;
  if(nState == 5) maxY = massFrame->GetMaximum()*0.4;
  //double lineWidth = 2.0;
  //TLine *lineSBLow = new TLine(sbLowMass, minY, sbLowMass, maxY);
  //TLine *lineSBHigh = new TLine(sbHighMass, minY, sbHighMass, maxY);
  //TLine *lineSigLow = new TLine(sigMinMass, minY, sigMinMass, maxY);
  //TLine *lineSigHigh = new TLine(sigMaxMass, minY, sigMaxMass, maxY);
  //lineSBLow->SetLineWidth(lineWidth);lineSBHigh->SetLineWidth(lineWidth);
  //lineSigLow->SetLineWidth(lineWidth);lineSigHigh->SetLineWidth(lineWidth);
  //lineSBLow->SetLineColor(kBlue);lineSBHigh->SetLineColor(kBlue);
  //lineSigLow->SetLineColor(kRed);lineSigHigh->SetLineColor(kRed);
  //lineSBLow->SetLineStyle(7);lineSBHigh->SetLineStyle(7);
  //lineSigLow->SetLineStyle(5);lineSigHigh->SetLineStyle(5);

  TH1* legendBlue = data->createHistogram("legendBlue",*JpsiMass,Binning(50)) ; legendBlue->SetLineColor(JpsiFitColor) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
  TH1* legendBlueDash = data->createHistogram("legendBlueDash",*JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(kBlue) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
  TH1* legendRed = data->createHistogram("legendRed",*JpsiMass,Binning(50)) ; legendRed->SetLineColor(onia::ColorNPJpsi) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
  TH1* legendBlack = data->createHistogram("legendBlack",*JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
  TH1* legendGreen = data->createHistogram("legendGreen",*JpsiMass,Binning(50)) ; legendGreen->SetLineColor(onia::ColorSumJpsiSignal) ; legendGreen->SetLineStyle(2) ; legendGreen->SetLineWidth(2) ;
  TH1* legendGreenDash = data->createHistogram("legendGreenDash",*JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(onia::ColorSumJpsiSignal) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
  TH1* legendPink = data->createHistogram("legendPink",*JpsiMass,Binning(50)) ; legendPink->SetLineColor(onia::ColorMuMuBG) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

  TLegend* MassLegend=new TLegend(0.78,0.35,0.905,0.5);
  MassLegend->SetFillColor(kWhite);
  MassLegend->SetTextFont(42);
  MassLegend->SetTextSize(0.035);
  MassLegend->SetBorderSize(0.);
  MassLegend->AddEntry(legendBlue,"sum","l");
  MassLegend->AddEntry(legendGreen,"J/#psi","l");
  MassLegend->AddEntry(legendPink,"#mu#mu BG","l");

  double left=0.7, top=0.9, textSize=0.03;
  TLatex *latex=new TLatex();

  gStyle->SetPadBottomMargin(0.08); //0.12
  gStyle->SetPadLeftMargin(0.09); //0.12
  gStyle->SetPadRightMargin(0.035); //0.05
  gStyle->SetPadTopMargin(0.05); //0.05

  TCanvas *c1;
  TPad *pad1;
  TPad *pad2;

  for(int LinLog=0; LinLog<2; LinLog++){
    cout<<"LinLog "<<LinLog<<endl;

    left=0.65, top=0.9, textSize=0.03;
    latex->SetTextFont(42);
    latex->SetNDC(kTRUE);
    latex->SetTextSize(textSize);
    double step=textSize*1.6;

    c1=new TCanvas("c1","",1000,900);

    c1->cd();
    pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
    pad1->SetGridy();
    pad1->SetBottomMargin(0.2);
    pad1->Draw();
    c1->cd();
    pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
    pad2->Draw();

    if(LinLog==1) massFrame->SetMinimum(massFrame->GetMaximum()*4e-3);
    if(LinLog==1) massFrame->SetMaximum(massFrame->GetMaximum()*3.);

    pad2->cd(0);
    if(LinLog==0) pad2->SetLogy(false);
    if(LinLog==1) pad2->SetLogy(true);

    massFrame->Draw(); MassLegend->Draw();
    //lineSBLow->Draw("same"); lineSBHigh->Draw("same"); lineSigLow->Draw("same"); lineSigHigh->Draw("same");
    top=0.885; textSize=0.030; latex->SetTextSize(textSize);
    left=0.74;

    latex->SetTextColor(kRed);
    latex->DrawLatex(left,top,"J/#psi");
    latex->SetTextColor(kBlack);

    top-=step;
    textSize=0.03; latex->SetTextSize(textSize);


    if(rapBin==0) latex->DrawLatex(left,top,Form("%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin],onia::KinParticleChar,onia::rapForPTRange[onia::kNbRapForPTBins]));
    else if(rapBin==1) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[rapBin]));
    else latex->DrawLatex(left,top,Form("%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin-1],onia::KinParticleChar,onia::rapForPTRange[rapBin]));
    top-=step;
    if(ptBin==0)
      latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTBins[rapBin]]));
    else
      latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]));

    //left=0.675;
    top-=step;
    top-=step;
    if(ws->var("CBn_jpsi")->isConstant()){
      latex->DrawLatex(left,top,Form("n^{CB}_{#psi}  =  %.1f",cbN));
      top-=step;
    }

    latex->DrawLatex(left,top,Form("effective #sigma =  %.3f MeV",ws->var("var_massres")->getVal()*1000));
    top-=step;
    latex->DrawLatex(left,top,Form("#frac{B}{B+S} (#pm3#sigma)  =  %.3f",ws->var("var_frac_jpsi_BGInSR")->getVal()));


    left=0.15; top=0.885; textSize=0.03;
    latex->SetTextSize(textSize);
    latex->DrawLatex(left,top,Form("#chi^{2}/ndf = %.1f / %d", chi2_Mass, ndof_Mass));
    top-=step;
    if(!ws->var("CBmass_p0_jpsi")->isConstant()){
      latex->DrawLatex(left,top,Form("#mu_{#psi}   =  %.3f #pm %.3f MeV",Meanp0*1000, Meanp0Err*1000));
      top-=step;
    }
    if(!ws->var("CBmass_p1_jpsi")->isConstant()){
      latex->DrawLatex(left,top,Form("#mu_{p1}   =  %.3f #pm %.3f MeV",Meanp1*1000, Meanp1Err*1000));
      top-=step;
    }
    if(!ws->var("CBmass_p2_jpsi")->isConstant()){
      latex->DrawLatex(left,top,Form("#mu_{p2}   =  %.3f #pm %.3f MeV",Meanp2*1000, Meanp2Err*1000));
      top-=step;
    }

    if(!ws->var("CBsigma_p0_jpsi")->isConstant()){
      latex->DrawLatex(left,top,Form("#sigma^{p0}_{#psi}  =  %.3f #pm %.3f MeV",Sigmap0*1000, Sigmap0Err*1000));
      top-=step;
    }
    if(!ws->var("CBsigma_p1_jpsi")->isConstant()){
      latex->DrawLatex(left,top,Form("#sigma^{p1}_{#psi}  =  %.3f #pm %.3f MeV",Sigmap1*1000, Sigmap1Err*1000));
      top-=step;
    }
    if(!ws->var("CBsigma_p2_jpsi")->isConstant()){
      latex->DrawLatex(left,top,Form("#sigma^{p2}_{#psi}  =  %.3f #pm %.3f MeV",Sigmap2*1000, Sigmap2Err*1000));
      top-=step;
    }

    if(!ws->var("CBalpha_p0_jpsi")->isConstant()){
      latex->DrawLatex(left,top,Form("#alpha^{CB}_{#psi}  =  %.3f #pm %.3f",Alphap0, Alphap0Err));
      top-=step;
    }
    if(!ws->var("CBn_jpsi")->isConstant()){
      latex->DrawLatex(left,top,Form("n^{CB}_{#psi}  =  %.3f #pm %.3f",cbN, cbNErr));
      top-=step;
    }
    if(!ws->var("CBalpha_p1_jpsi")->isConstant()){
      latex->DrawLatex(left,top,Form("#alpha_{p1}  =  %.3f #pm %.3f",Alphap1, Alphap1Err));
      top-=step;
    }
    if(!ws->var("bkgLambda_jpsi")->isConstant()){
      latex->DrawLatex(left,top,Form("#lambda_{BG} =  %.3f #pm %.3f GeV^{-1}",-1.*lambda, lambdaErr)); // change to has positive value
      top-=step;
    }
    if(!ws->var("fracBkg_jpsi")->isConstant()){
      latex->DrawLatex(left,top,Form("f^{#psi}_{#mu#muBG}  =  %.3f #pm %.3f",fracBkg, fracBkgErr));
      top-=step;
    }


    pad1->cd(0); pad1->SetLogy(0);
    massFramePull->Draw();

    std::stringstream saveMass;
    if(LinLog==0) saveMass << "Fit/jpsiFit/mass_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
    if(LinLog==1) saveMass << "Fit/jpsiFit/mass_log_rap" << rapBin << "_pt" << ptBin << ".pdf";

    c1->SaveAs(saveMass.str().c_str());

  }

  delete c1;
  delete legendBlue;
  delete legendBlueDash;
  delete legendRed;
  delete legendBlack;
  delete legendGreen;
  delete legendGreenDash;
  delete legendPink;

  return chi2_Mass/ndof_Mass;
}



//==============================================
double plotMassRap(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom){
  int  nbins=90; //0.005 bin size
  TGaxis::SetMaxDigits(3);
  gSystem->mkdir("Fit/jpsiFit",kTRUE);

  const int nRapSteps=4;
  double RapStepBorders[nRapSteps+1]={0.,0.3,0.6,0.9,1.2};
  const int colorRap[nRapSteps+2]={632, 632, 632, 600, 418, 616};


  nbins=32;

  double binWidth=(onia::massMax-onia::massMin)/double(nbins)*1000;


  RooRealVar *JpsiMass = ws->var("JpsiMass");
  assert( 0 != JpsiMass );
  RooRealVar *JpsiRap = ws->var("JpsiRap");
  assert( 0 != JpsiRap );

  RooPlot *massFrame = JpsiMass->frame(Bins(nbins));
  assert ( 0 != massFrame );
  massFrame->SetName(Form("mass_plot_rap%d_pt%d",rapBin,ptBin));
  massFrame->SetTitle("");
  massFrame->GetYaxis()->SetTitle(Form("Events / %1.0f MeV",binWidth));
  massFrame->GetYaxis()->SetTitleOffset(1.3);


  RooAbsData *data= ws->data(Form("jpsi_data_rap%d_pt%d",rapBin,ptBin));
  assert ( 0 != data );
  RooFitResult* fitRlt = dynamic_cast<RooFitResult*>(ws->obj(Form("m_fitresult_rap%d_pt%d",rapBin,ptBin)));
  assert ( 0 != fitRlt);


  double plotRapMax=1.2;
  TH1* hist_rap_vs_mass = (TH1*)data->createHistogram("hist",*JpsiRap,Binning(30,-plotRapMax, plotRapMax),YVar(*JpsiMass,Binning(90,onia::massMin, onia::massMax)));

  char XTitle[200];
  char YTitle[200];
  sprintf(YTitle,"M^{J/#psi} [GeV]");
  sprintf(XTitle,"y^{J/#psi}");
  double yOffset=1.5;
  hist_rap_vs_mass->GetYaxis()->SetTitleOffset(yOffset);
  hist_rap_vs_mass->SetStats(0);
  hist_rap_vs_mass->SetTitle(0);
  hist_rap_vs_mass->GetYaxis()->SetTitle(YTitle);
  hist_rap_vs_mass->GetXaxis()->SetTitle(XTitle);

  hist_rap_vs_mass->GetYaxis()->SetRangeUser(onia::massMin, onia::massMax);
  hist_rap_vs_mass->GetXaxis()->SetRangeUser(-plotRapMax, plotRapMax);


  TF1 *RegionLine[4];

  char Formula[200];
  sprintf(Formula,"(%f-%f*(%f+%f*abs(x)+%f*x*x))",ws->var("CBmass_p0_jpsi")->getVal(), onia::nSigMass, ws->var("CBsigma_p0_jpsi")->getVal(), ws->var("CBsigma_p1_jpsi")->getVal(), ws->var("CBsigma_p2_jpsi")->getVal());
  RegionLine[1]=new TF1("RegionLine1",Formula,-plotRapMax, plotRapMax);
  RegionLine[1]->SetLineColor(kRed);
  RegionLine[1]->SetLineWidth( 3 );
  RegionLine[1]->SetLineStyle( 2 );
  sprintf(Formula,"(%f+%f*(%f+%f*abs(x)+%f*x*x))",ws->var("CBmass_p0_jpsi")->getVal(), onia::nSigMass, ws->var("CBsigma_p0_jpsi")->getVal(), ws->var("CBsigma_p1_jpsi")->getVal(), ws->var("CBsigma_p2_jpsi")->getVal());
  RegionLine[2]=new TF1("RegionLine2",Formula,-plotRapMax, plotRapMax);
  RegionLine[2]->SetLineColor(kRed);
  RegionLine[2]->SetLineWidth( 3 );
  RegionLine[2]->SetLineStyle( 2 );
  sprintf(Formula,"(%f-%f*(%f+%f*abs(x)+%f*x*x))",ws->var("CBmass_p0_jpsi")->getVal(), onia::nSigBkgLow, ws->var("CBsigma_p0_jpsi")->getVal(), ws->var("CBsigma_p1_jpsi")->getVal(), ws->var("CBsigma_p2_jpsi")->getVal());
  RegionLine[3]=new TF1("RegionLine3",Formula,-plotRapMax, plotRapMax);
  RegionLine[3]->SetLineColor(kBlack);
  RegionLine[3]->SetLineWidth( 3 );
  RegionLine[3]->SetLineStyle( 2 );
  sprintf(Formula,"(%f+%f*(%f+%f*abs(x)+%f*x*x))",ws->var("CBmass_p0_jpsi")->getVal(), onia::nSigBkgHigh, ws->var("CBsigma_p0_jpsi")->getVal(), ws->var("CBsigma_p1_jpsi")->getVal(), ws->var("CBsigma_p2_jpsi")->getVal());
  RegionLine[4]=new TF1("RegionLine4",Formula,-plotRapMax, plotRapMax);
  RegionLine[4]->SetLineColor(kBlack);
  RegionLine[4]->SetLineWidth( 3 );
  RegionLine[4]->SetLineStyle( 2 );

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

  TCanvas *c2=new TCanvas("c2","",800,700);

  hist_rap_vs_mass->Draw("colz");
  RegionLine[1]->Draw("same");
  RegionLine[2]->Draw("same");
  RegionLine[3]->Draw("same");
  RegionLine[4]->Draw("same");


  char ptchar[200];
  char rapchar[200];

  if(ptBin==0)
    sprintf(ptchar,"%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTBins[rapBin]]);
  else
    sprintf(ptchar,"%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]);

  if(rapBin==0) sprintf(rapchar,"%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin],onia::KinParticleChar,onia::rapForPTRange[onia::kNbRapForPTBins]);
  else if(rapBin==1) sprintf(rapchar,"|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[rapBin]);
  else sprintf(rapchar,"%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin-1],onia::KinParticleChar,onia::rapForPTRange[rapBin]);


  top=0.965; left=0.35;

  latex->DrawLatex(left,top,Form("J/#psi, %s, %s",ptchar, rapchar));

  std::stringstream saveMassRap;
  saveMassRap << "Fit/jpsiFit/rap_vs_mass_rap" << rapBin << "_pt" << ptBin << ".pdf";
  c2->SaveAs(saveMassRap.str().c_str());
  std::stringstream saveMassRapLog;
  saveMassRapLog << "Fit/jpsiFit/rap_vs_mass_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
  c2->SetLogz(true);
  c2->SaveAs(saveMassRapLog.str().c_str());


  RooAbsPdf *massPdf = ws->pdf("massModel_jpsi");
  assert ( 0 != massPdf );
  RooAbsPdf *bkgMassShape = ws->pdf("bkgMassShape_jpsi");
  assert ( 0 != bkgMassShape );
  RooAbsPdf *sigMassShape_jpsi = ws->pdf("sigMassShape_jpsi");
  assert ( 0 != sigMassShape_jpsi );

  RooRealVar new_fBG("new_fBG","new_fBG",ws->var("fracBkg_jpsi")->getVal()/1.36);
  RooAddPdf FullmassPdf= RooAddPdf("FullmassPdf","FullmassPdf",RooArgList(*bkgMassShape,*sigMassShape_jpsi),RooArgList(new_fBG));

  TLegend* MassLegend=new TLegend(0.7,0.55,0.88,0.75);
  MassLegend->SetFillColor(kWhite);
  MassLegend->SetTextFont(42);
  MassLegend->SetTextSize(0.035);
  MassLegend->SetBorderSize(0.);

  char legendchar[200];


  RooAbsData* dataRapStep[nRapSteps];
  TH1* legendPhantom[nRapSteps];
  for(int iRapStep=0;iRapStep<nRapSteps;iRapStep++){
    cout<<"iRapStep "<<iRapStep<<endl;
    cout<<"color "<<colorRap[iRapStep+2]<<endl;

    std::stringstream cutRap;
    cutRap << "TMath::Abs(JpsiRap) > " << RapStepBorders[iRapStep] << " && TMath::Abs(JpsiRap) < " << RapStepBorders[iRapStep+1];
    dataRapStep[iRapStep] = data->reduce(Cut(cutRap.str().c_str()));
    JpsiRap->setRange(Form("RapRegion%d_pos",iRapStep),RapStepBorders[iRapStep],RapStepBorders[iRapStep+1]);
    JpsiRap->setRange(Form("RapRegion%d_neg",iRapStep),-RapStepBorders[iRapStep+1],-RapStepBorders[iRapStep]);

    dataRapStep[iRapStep]->plotOn(massFrame,MarkerSize(0.8), MarkerColor(colorRap[iRapStep+2]));
    int nEntries = dataRapStep[iRapStep]->numEntries();

    RooDataSet *dataJpsiRap = (RooDataSet*)dataRapStep[iRapStep]->reduce(SelectVars(RooArgSet(*JpsiRap)),Name("dataJpsiRap"));

    int nevt=1000000;
    //int nevt=1000;
    cout<<"generating datasets"<<endl;
    RooDataSet *SignalPseudoData = sigMassShape_jpsi->generate(*JpsiMass,ProtoData(*dataJpsiRap),NumEvents(nevt));
    cout<<"finished generating datasets"<<endl;

    int nbinsHists=100;
    TH1F* sigMassShape_jpsi_asHist = new TH1F("sigMassShape_jpsi_asHist","sigMassShape_jpsi_asHist", nbinsHists, onia::massMin, onia::massMax);
    SignalPseudoData->fillHistogram(sigMassShape_jpsi_asHist,RooArgList(*JpsiMass));
    RooDataHist* signalMassShape_asRooDataHist = new RooDataHist("signalMassShape_asRooDataHist","signalMassShape_asRooDataHist", RooArgList(*JpsiMass), sigMassShape_jpsi_asHist);
    RooHistPdf* signalMassShape_asHistPdf = new RooHistPdf("signalMassShape_asHistPdf","signalMassShape_asHistPdf", RooArgSet(*JpsiMass), *signalMassShape_asRooDataHist, 3);

    RooRealVar new_fBG("new_fBG","new_fBG",ws->var("fracBkg_jpsi")->getVal());
    RooAddPdf FullmassPdf= RooAddPdf("FullmassPdf","FullmassPdf",RooArgList(*bkgMassShape,*signalMassShape_asHistPdf),RooArgList(new_fBG));

    FullmassPdf.plotOn(massFrame,
                       LineWidth(2), LineColor(colorRap[iRapStep+2]),
                       Normalization(nEntries,2)
                       );


    legendPhantom[iRapStep] = dataRapStep[iRapStep]->createHistogram(Form("legendRap%d",iRapStep),*JpsiMass,Binning(50)) ; legendPhantom[iRapStep]->SetLineColor(colorRap[iRapStep+2]) ; legendPhantom[iRapStep]->SetLineStyle(kSolid) ; legendPhantom[iRapStep]->SetLineWidth(2) ;
    sprintf(legendchar,"%1.1f < |y^{J/#psi}| <  %1.1f", RapStepBorders[iRapStep],RapStepBorders[iRapStep+1]);
    MassLegend->AddEntry(legendPhantom[iRapStep],legendchar,"l");


  }



  gStyle->SetPadBottomMargin(0.08); //0.12
  gStyle->SetPadLeftMargin(0.12); //0.12
  gStyle->SetPadRightMargin(0.02); //0.05
  gStyle->SetPadTopMargin(0.05); //0.05


  TCanvas *c1=new TCanvas("c1","",800,700);

  massFrame->Draw(); MassLegend->Draw();
  left=0.7; top=0.90; textSize=0.030; latex->SetTextSize(textSize);
  //if(nState == 4)
  //	latex->DrawLatex(left,top,"J/#psi");
  //if(nState == 5)
  //	latex->DrawLatex(left,top,"#psi(2S)");
  //top-=step;
  textSize=0.03; latex->SetTextSize(textSize);


  if(ptBin==0)
    latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTBins[rapBin]]));
  else
    latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]));


  std::stringstream saveMass;
  saveMass << "Fit/jpsiFit/mass_rapDep_rap" << rapBin << "_pt" << ptBin << ".pdf";
  c1->SaveAs(saveMass.str().c_str());



  delete c1;
  return 0.;
}



//==============================================
double plotLifeSig(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom){
  int nbins=onia::LifetimePlotBins;
  TGaxis::SetMaxDigits(3);

  double PlotMin, PlotMax;
  if(zoom){
    PlotMin=onia::ctPlotMinZoom;
    PlotMax=onia::ctPlotMaxZoom;
  }
  else{
    PlotMin=onia::ctPlotMin;
    PlotMax=onia::ctPlotMax;
  }

  double binWidth=(PlotMax-PlotMin)/double(nbins)*1000;

  //ws->var("jpsi_promptMean")->setMin(-10.);
  //ws->var("jpsi_promptMean")->setVal(-1.);

  RooRealVar JpsiMass(*ws->var("JpsiMass"));
  RooRealVar Jpsict(*ws->var("Jpsict"));
  RooRealVar *JpsictErr = ws->var("JpsictErr");
  RooFormulaVar *JpsictErrDeformed = (RooFormulaVar*)ws->function("JpsictErrDeformed");

  RooPlot *ctauFrame=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(PlotMin, PlotMax));
  ctauFrame->SetName(Form("ctausig_plot_rap%d_pt%d",rapBin,ptBin));
  ctauFrame->GetYaxis()->SetTitle(Form("Events / %1.0f micron",binWidth));
  ctauFrame->SetTitle("");
  ctauFrame->GetYaxis()->SetTitleOffset(1.3);

  RooPlot *ctauFramePull=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(PlotMin, PlotMax));
  ctauFramePull->SetName(Form("pullctausig_plot_rap%d_pt%d",rapBin,ptBin));
  ctauFramePull->SetTitle("");
  ctauFramePull->GetYaxis()->SetTitle("pull");
  ctauFramePull->GetXaxis()->SetTitleSize(0.08);
  ctauFramePull->GetYaxis()->SetTitleSize(0.08);
  ctauFramePull->GetXaxis()->SetLabelSize(0.08);
  ctauFramePull->GetYaxis()->SetLabelSize(0.08);
  ctauFramePull->GetYaxis()->SetTitleOffset(0.4);
  ctauFramePull->GetYaxis()->SetRangeUser(-5.99,5.99);


  RooDataSet *dataSR=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));
  RooAbsData* dataSRProj;
  if(SpeedPlotting){
    //dataSR1Proj = dataSR1->reduce(EventRange(0,nSpeedPlotting));
    dataSRProj = dataSR->reduce(EventRange(0,nSpeedPlotting));
  }
  else{
    //dataSR1Proj = dataSR1->reduce(EventRange(0,1e10));
    dataSRProj = dataSR->reduce(EventRange(0,1e10));
  }
  cout<<"number of events in dataSR = "<<dataSR->numEntries()<<endl;
  cout<<"number of events in dataSRProj = "<<dataSRProj->numEntries()<<endl;

  double nentries=dataSR->numEntries();

  double minY = 0.;
  double maxY = 0.;
  double enlargeYby=onia::enlargeYby_ML;


  RooFitResult* fitRlt = (RooFitResult*)ws->obj(Form("jpsi_l_fitresult_rap%d_pt%d",rapBin,ptBin));
  if(!fitRlt) { cout<<">>=======Error: no Fit result object in workspace=========="<<endl; return 0.; }
  fitRlt->Print();
  cout<<"dataSR->numEntries: "<<dataSR->numEntries()<<endl;

  bool correctResolutionForPlotting=false;
  double resCorrFactor=1.0525;
  if(ptBin==0) resCorrFactor=1.0825;
  if(ptBin==1) resCorrFactor=1.0725;
  if(correctResolutionForPlotting){
    ws->var("jpsi_ctResolution")->setVal(ws->var("jpsi_ctResolution")->getVal()*resCorrFactor);
    ws->var("jpsi_ctResolution2")->setVal(ws->var("jpsi_ctResolution2")->getVal()*resCorrFactor);
  }

  RooRealVar *promptMean_=(RooRealVar*)ws->var("jpsi_promptMean");
  RooRealVar *ctResolution_=(RooRealVar*)ws->var("jpsi_ctResolution");
  RooRealVar *ctResolution2_=(RooRealVar*)ws->var("jpsi_ctResolution2");
  RooRealVar *fracGauss2_=(RooRealVar*)ws->var("jpsi_fracGauss2");
  //RooRealVar *fracGauss3_=(RooRealVar*)ws->var("jpsi_fracGauss3");

  RooRealVar *nonPromptTau_=(RooRealVar*)ws->var("jpsi_nonPromptTau");
  RooRealVar *fBkg_ = (RooRealVar*)ws->var("jpsi_fBkg");
  RooRealVar *fPrompt_ = (RooRealVar*)ws->var("jpsi_fPrompt");

  double promptMean = promptMean_->getVal();
  double promptMeanErr = promptMean_->getError();
  double promptCtRe = ctResolution_->getVal();
  double promptCtReErr = ctResolution_->getError();
  double promptCtRe2 = ctResolution2_->getVal();
  double promptCtRe2Err = ctResolution2_->getError();
  double fracGauss2 = fracGauss2_->getVal();
  double fracGauss2Err = fracGauss2_->getError();
  //double fracGauss3 = fracGauss3_->getVal();
  //double fracGauss3Err = fracGauss3_->getError();
  double nonPromptTau = nonPromptTau_->getVal();
  double nonPromptTauErr = nonPromptTau_->getError();

  double fBkg = fBkg_->getVal();
  double fBkgErr = fBkg_->getError();
  double fPrompt = fPrompt_->getVal();
  double fPromptErr = fPrompt_->getError();
  double fNonPrompt =  1.-fBkg-fPrompt;

  cout<<"promptMean: "<<promptMean<<endl;
  cout<<"promptCtRe: "<<promptCtRe<<endl;
  cout<<"nonPromptTau: "<<nonPromptTau<<endl;
  cout<<"fBkg: "<<fBkg<<endl;
  cout<<"fBkgErr: "<<fBkgErr<<endl;
  cout<<"fPrompt: "<<fPrompt<<endl;
  cout<<"fNonPrompt: "<<fNonPrompt<<endl;

  RooAddPdf *ModelLife = (RooAddPdf*)ws->pdf("jpsi_fulllifetimeSR");

  RooAddModel *Prompt = (RooAddModel*)ws->pdf("jpsi_TotalPromptLifetime");
  RooAbsPdf *nonPromptSSD = (RooAbsPdf*)ws->pdf("jpsi_nonPromptSSD");
  RooAbsPdf *backgroundlifetime = (RooAbsPdf*)ws->pdf("jpsi_backgroundlifetime");

  //RooAbsPdf *promptLifetime = (RooAbsPdf*)ws->pdf("jpsi_promptLifetime");
  //RooAbsPdf *promptLifetime2 = (RooAbsPdf*)ws->pdf("jpsi_promptLifetime2");
  //RooAbsPdf *promptLifetime3 = (RooAbsPdf*)ws->pdf("jpsi_promptLifetime3");










  int parsFit;
  //ploting signal region
  dataSR->plotOn(ctauFrame,MarkerSize(onia::markerSize_ML), Name("myHist"));

  maxY = ctauFrame->GetMaximum()*enlargeYby;
  minY = 5e-1;
  ctauFrame->GetYaxis()->SetRangeUser(minY,maxY);

  int JpsiFitColor=633;

  ModelLife->plotOn(ctauFrame,
                    ProjWData(*JpsictErr, *dataSRProj),
                    LineWidth(2),
                    LineColor(JpsiFitColor),
                    NumCPU(1),
                    Name("myCurve"));



  //------get chi2------------SHOULD DONE after PLOTTING------
  parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.

  int nBins_LSig=ctauFrame->GetNbinsX();
  double chi2Pre_LSig=ctauFrame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
  int ndof_LSig=nBins_LSig-parsFit;  //num of degree of freedom
  double chi2_LSig=chi2Pre_LSig*ndof_LSig;

  //RooHist* hpull_ctauSig = ctauFrame->pullHist() ;
  RooHist* hpull_ctauSig = ctauFrame->pullHist("myHist","myCurve",kTRUE);

  /*
    manual chi2 calculation

    RooHist* hresid_ctauSig = ctauFrame->residHist("myHist","myCurve",true,kTRUE);

    double chi2manual=0;
    for(int iX=0;iX<hresid_ctauSig->GetN()+1;iX++){
    double buffx, buffy;
    hresid_ctauSig->GetPoint(iX, buffx, buffy);
    chi2manual+=buffy*buffy;
    }

    chi2manual=TMath::Sqrt(chi2manual);
    double chi2ndfmanual=chi2manual/double(ndof_LSig);
    cout<<"chi2manual "<<chi2manual<<endl;
    cout<<"ndof_LSig "<<ndof_LSig<<endl;
    cout<<"chi2ndfmanual "<<chi2ndfmanual<<endl;

  */

  //save pull distribution
  TH1F* pull = new TH1F("pull","pull distribution", 100,-10.,10.);
  gSystem->mkdir("Fit/pull",kTRUE);gSystem->mkdir("Fit/jpsiFit",kTRUE);
  //TFile *pullFile = new TFile(Form("Fit/pull/pullAverage_rap%d_pt%d_SR.root",rapBin,ptBin),"RECREATE");
  TFile *pullFile = new TFile(Form("Fit/pull/pull_rap%d_pt%d_SR.root",rapBin,ptBin),"RECREATE");

  hpull_ctauSig->SetMarkerSize(0.8);
  for(int i=0;i<hpull_ctauSig->GetN();i++){
    hpull_ctauSig->SetPointEYlow(i,0.);
    hpull_ctauSig->SetPointEYhigh(i,0.);
    double x,y;
    hpull_ctauSig->GetPoint(i,x,y);
    pull->Fill(y);
  }
  pullFile->cd();
  pull->Write();
  pullFile->Close();

  ctauFramePull->addPlotable(hpull_ctauSig,"P");


  //Prompt->plotOn(ctauFrame,
  //		ProjWData(*JpsictErr, *dataSRProj),
  //		Normalization(fPrompt),
  //		LineStyle(5),
  //		LineColor(onia::ColorPRJpsi),
  //		LineWidth(2), NumCPU(1));
  //
  //promptLifetime->plotOn(ctauFrame,
  //		ProjWData(*JpsictErr, *dataSRProj),
  //		Normalization(fPrompt*(1.-fracGauss2)),
  //		LineStyle(5),
  //		LineColor(onia::ColorPRJpsi),
  //		LineWidth(1.), NumCPU(1));
  //promptLifetime2->plotOn(ctauFrame,
  //		ProjWData(*JpsictErr, *dataSRProj),
  //		Normalization(fPrompt*fracGauss2),
  //		LineStyle(5),
  //		LineColor(onia::ColorPRJpsi),
  //		LineWidth(1.), NumCPU(1));
  //
  //nonPromptSSD->plotOn(ctauFrame,
  //		ProjWData(*JpsictErr, *dataSRProj),
  //		Normalization(fNonPrompt),
  //		LineStyle(2),
  //		LineColor(onia::ColorNPJpsi),
  //		LineWidth(2), NumCPU(1));
  //
  //backgroundlifetime->plotOn(ctauFrame,
  //		ProjWData(*JpsictErr, *dataSRProj),
  //		Normalization(fBkg),
  //		LineStyle(7),
  //		LineColor(onia::ColorMuMuBG),
  //		LineWidth(2), NumCPU(1));

  ModelLife->plotOn(ctauFrame,
                    Components("jpsi_TotalPromptLifetime"),
                    ProjWData(*JpsictErr, *dataSRProj),
                    LineStyle(5),
                    LineColor(onia::ColorPRJpsi),
                    LineWidth(2), NumCPU(1));

  ModelLife->plotOn(ctauFrame,
                    Components("jpsi_nonPromptSSD"),
                    ProjWData(*JpsictErr, *dataSRProj),
                    LineStyle(2),
                    LineColor(onia::ColorNPJpsi),
                    LineWidth(2), NumCPU(1));

  ModelLife->plotOn(ctauFrame,
                    Components("jpsi_backgroundlifetime"),
                    ProjWData(*JpsictErr, *dataSRProj),
                    LineStyle(7),
                    LineColor(onia::ColorMuMuBG),
                    LineWidth(2), NumCPU(1));



  ////// JpsictErr plotting test code
  //RooDataSet *dataSRcterr = (RooDataSet*)dataSR->reduce(SelectVars(RooArgSet(*JpsictErr)),Name("dataSRcterr"));
  //RooAbsReal *MLNLLSR = NULL;
  //
  //int nevt=nentries*fPrompt;
  //
  //
  //cout<<"generating datasets"<<endl;
  //RooDataSet *SignalPseudoData = Prompt->generate(Jpsict,ProtoData(*dataSRcterr),NumEvents(nevt));
  //cout<<"finished generating datasets"<<endl;
  //
  //MLNLLSR = (RooAbsReal *)Prompt->createNLL(*SignalPseudoData,
  //		ConditionalObservables(RooArgSet(*ws->var("JpsictErr"))),
  //		Extended(kFALSE),
  //		NumCPU(6));
  //
  //RooMinuit *lMinuit = new RooMinuit(*MLNLLSR);
  //
  //lMinuit->setStrategy(1);
  //lMinuit->setPrintEvalErrors(-1);
  //lMinuit->setEvalErrorWall(false);
  //lMinuit->setVerbose(false);
  //lMinuit->setPrintLevel(-1);
  //
  //lMinuit->migrad();
  //
  //RooFitResult *fitresult = (RooFitResult*)lMinuit->save("tmp");
  //fitresult->Print();
  //
  //
  //SignalPseudoData->plotOn(ctauFrame,MarkerSize(onia::markerSize_ML),MarkerColor(kRed+0));
  //Prompt->plotOn(ctauFrame,
  //		ProjWData(*JpsictErr, *dataSRProj),
  //		//Normalization(fPrompt),
  //		LineStyle(5),
  //		LineColor(kGreen+0),
  //		LineWidth(2), NumCPU(1));


  //int nbinsHists=10000;
  //TH1F* promptlifetime_asHist = (TH1F*)Prompt->createHistogram("promptlifetime_asHist",Jpsict,ConditionalObservables(*JpsictErr),Binning(nbinsHists, PlotMin, PlotMax));
  //RooDataHist* promptlifetime_asRooDataHist = new RooDataHist("promptlifetime_asRooDataHist","promptlifetime_asRooDataHist", RooArgList(Jpsict), promptlifetime_asHist);
  //RooHistPdf* jpsi_promptlifetime_asHistPdf = new RooHistPdf("jpsi_promptlifetime_asHistPdf","jpsi_promptlifetime_asHistPdf", RooArgSet(Jpsict), *promptlifetime_asRooDataHist, 1);
  //ws->import(*jpsi_promptlifetime_asHistPdf);
  //
  //ws->factory("SUM::jpsi_fulllifetimeSRnoC_PrHistPdf(jpsi_fBkg*jpsi_backgroundlifetime,jpsi_fPrompt*jpsi_promptlifetime_asHistPdf,jpsi_nonPromptSSD)");
  //
  //ws->pdf("jpsi_fulllifetimeSRnoC_PrHistPdf")->plotOn(ctauFrame,
  //		ProjWData(*JpsictErr, *dataSRProj),
  //		LineStyle(5),
  //		LineColor(kGreen),
  //		LineWidth(2), NumCPU(1));
  //
  //jpsi_promptlifetime_asHistPdf->plotOn(ctauFrame,
  //		Normalization(fPrompt),
  //		LineStyle(5),
  //		LineColor(onia::ColorSumJpsiSignal),
  //		LineWidth(2), NumCPU(1));








  double Ymax = ctauFrame->GetMaximum();
  cout<<"Ymax: "<<Ymax<<endl;
  ctauFrame->SetMaximum(3*Ymax);
  ctauFrame->SetMinimum(1.1);


  TH1* legendBlue = dataSR->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(JpsiFitColor) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
  TH1* legendBlueDash = dataSR->createHistogram("legendBlueDash",JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(onia::ColorPRJpsi) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
  TH1* legendRed = dataSR->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(onia::ColorNPJpsi) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
  TH1* legendBlack = dataSR->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
  TH1* legendGreen = dataSR->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;
  TH1* legendGreenDash = dataSR->createHistogram("legendGreenDash",JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(kGreen) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
  TH1* legendPink = dataSR->createHistogram("legendPink",JpsiMass,Binning(50)) ; legendPink->SetLineColor(onia::ColorMuMuBG) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;


  TLegend* LifetimeLegendSig=new TLegend(0.13,0.75,0.23,0.91);
  LifetimeLegendSig->SetFillColor(kWhite);
  LifetimeLegendSig->SetTextFont(42);
  LifetimeLegendSig->SetTextSize(0.035);
  LifetimeLegendSig->SetBorderSize(0.);
  LifetimeLegendSig->AddEntry(legendBlue,"sum","l");
  LifetimeLegendSig->AddEntry(legendBlueDash,"PR J/#psi","l");
  LifetimeLegendSig->AddEntry(legendRed,"NP J/#psi","l");
  LifetimeLegendSig->AddEntry(legendPink,"#mu#mu BG","l");

  double left=0.5, top=0.9, textSize=0.03;
  TLatex *latex=new TLatex();
  latex->SetTextFont(42);
  latex->SetNDC(kTRUE);
  latex->SetTextSize(textSize);
  double step=textSize*1.6;

  gStyle->SetPadBottomMargin(0.08); //0.12
  gStyle->SetPadLeftMargin(0.09); //0.12
  gStyle->SetPadRightMargin(0.035); //0.05
  gStyle->SetPadTopMargin(0.05); //0.05


  TCanvas *c1;
  TPad *pad1;
  TPad *pad2;

  ctauFrame->SetMinimum(minY);
  ctauFrame->SetMaximum(maxY);

  for(int LinLog=0; LinLog<2; LinLog++){
    cout<<"LinLog "<<LinLog<<endl;

    c1=new TCanvas("c1","",1000,900);

    c1->cd();
    pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
    pad1->SetGridy();
    pad1->SetBottomMargin(0.2);
    pad1->Draw();
    c1->cd();
    pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
    pad2->Draw();

    //Sig
    pad2->cd(0);
    if(LinLog==0) pad2->SetLogy(false);
    if(LinLog==1) pad2->SetLogy(true);

    ctauFrame->Draw();
    LifetimeLegendSig->Draw();

    left=0.54; top=0.885; textSize=0.030; latex->SetTextSize(textSize);
    //if(nState == 4) latex->DrawLatex(left,top,"J/#psi");
    //if(nState == 5) latex->DrawLatex(left,top,"#psi(2S)");
    //top-=step;
    textSize=0.03; latex->SetTextSize(textSize);

    if(rapBin==0) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[onia::kNbRapForPTBins]));
    else if(rapBin==1) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[rapBin]));
    else latex->DrawLatex(left,top,Form("%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin-1],onia::KinParticleChar,onia::rapForPTRange[rapBin]));
    top-=step;
    if(ptBin==0)
      latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTBins[rapBin]]));
    else
      latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]));

    top-=step;
    latex->SetTextColor(kRed);
    latex->DrawLatex(left,top,"J/#psi SR");
    latex->SetTextColor(kBlack);
    top-=step;
    latex->DrawLatex(left,top,Form("#chi^{2}/ndf = %.1f/%d",chi2_LSig,ndof_LSig));
    top-=step;

    latexFloatingLifetimePars(ws, latex, 0);

    pad1->cd(0); pad1->SetLogy(0);
    ctauFramePull->Draw();

    std::stringstream saveCtau;
    if(!zoom){
      if(LinLog==0) saveCtau << "Fit/jpsiFit/ctau_SR_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
      if(LinLog==1) saveCtau << "Fit/jpsiFit/ctau_SR_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
    }
    else{
      if(LinLog==0) saveCtau << "Fit/jpsiFit/ctau_SR_zoom_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
      if(LinLog==1) saveCtau << "Fit/jpsiFit/ctau_SR_zoom_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
    }

    c1->SaveAs(saveCtau.str().c_str());

  }

  if(correctResolutionForPlotting){
    ws->var("jpsi_ctResolution")->setVal(ws->var("jpsi_ctResolution")->getVal()/resCorrFactor);
    ws->var("jpsi_ctResolution2")->setVal(ws->var("jpsi_ctResolution2")->getVal()/resCorrFactor);
  }

  delete c1;
  delete legendBlue;
  delete legendBlueDash;
  delete legendRed;
  delete legendBlack;
  delete legendGreen;
  delete legendGreenDash;
  delete legendPink;
  return chi2_LSig/ndof_LSig;
}

//==============================================
double plotLifeBg(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom, int SB){
  int nbins=onia::LifetimePlotBins;
  TGaxis::SetMaxDigits(3);

  double PlotMin, PlotMax;
  if(zoom){
    PlotMin=onia::ctPlotMinZoom;
    PlotMax=onia::ctPlotMaxZoom;
  }
  else{
    PlotMin=onia::ctPlotMin;
    PlotMax=onia::ctPlotMax;
  }

  char SBchar[200];
  if(SB==1) sprintf(SBchar,"LSB");
  if(SB==2) sprintf(SBchar,"RSB");

  double binWidth=(PlotMax-PlotMin)/double(nbins)*1000;


  RooRealVar JpsiMass(*ws->var("JpsiMass"));
  RooRealVar Jpsict(*ws->var("Jpsict"));
  RooRealVar *JpsictErr = ws->var("JpsictErr");

  RooPlot *ctauFrame=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(PlotMin, PlotMax));
  ctauFrame->SetName(Form("ctausig_plot_rap%d_pt%d",rapBin,ptBin));
  ctauFrame->GetYaxis()->SetTitle(Form("Events / %1.0f micron",binWidth));
  ctauFrame->SetTitle("");
  //ctauFrame->SetXTitle("lifetime [mm]");
  ctauFrame->GetYaxis()->SetTitleOffset(1.3);

  RooPlot *ctauFramePull=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(PlotMin, PlotMax));
  ctauFramePull->SetName(Form("pullctausig_plot_rap%d_pt%d",rapBin,ptBin));
  ctauFramePull->SetTitle("");
  //ctauFramePull->SetXTitle("lifetime [mm]");
  ctauFramePull->GetYaxis()->SetTitle("pull");
  ctauFramePull->GetXaxis()->SetTitleSize(0.08);
  ctauFramePull->GetYaxis()->SetTitleSize(0.08);
  ctauFramePull->GetXaxis()->SetLabelSize(0.08);
  ctauFramePull->GetYaxis()->SetLabelSize(0.08);
  ctauFramePull->GetYaxis()->SetTitleOffset(0.4);
  ctauFramePull->GetYaxis()->SetRangeUser(-5.99,5.99);


  RooDataSet *dataSB=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_%s",rapBin,ptBin,SBchar));
  RooAbsData* dataSBProj;
  if(SpeedPlotting){
    dataSBProj = dataSB->reduce(EventRange(0,nSpeedPlotting));
  }
  else{
    dataSBProj = dataSB->reduce(EventRange(0,1e10));
  }
  cout<<"number of events in dataSB = "<<dataSB->numEntries()<<endl;
  cout<<"number of events in dataSBProj = "<<dataSBProj->numEntries()<<endl;

  double minY = 0.;
  double maxY = 0.;
  double enlargeYby=onia::enlargeYby_ML;


  RooFitResult* fitRlt = (RooFitResult*)ws->obj(Form("jpsi_l_fitresult_rap%d_pt%d",rapBin,ptBin));
  if(!fitRlt) { cout<<">>=======Error: no Fit result object in workspace=========="<<endl; return 0.; }
  fitRlt->Print();
  cout<<"dataSB->numEntries: "<<dataSB->numEntries()<<endl;

  bool correctResolutionForPlotting=false;
  double resCorrFactor=1.0525;
  if(ptBin==0) resCorrFactor=1.0825;
  if(ptBin==1) resCorrFactor=1.0725;
  if(correctResolutionForPlotting){
    ws->var("jpsi_ctResolution")->setVal(ws->var("jpsi_ctResolution")->getVal()*resCorrFactor);
    ws->var("jpsi_ctResolution2")->setVal(ws->var("jpsi_ctResolution2")->getVal()*resCorrFactor);
  }

  RooRealVar *promptMean_=(RooRealVar*)ws->var("jpsi_promptMean");
  RooRealVar *ctResolution_=(RooRealVar*)ws->var("jpsi_ctResolution");
  RooRealVar *ctResolution2_=(RooRealVar*)ws->var("jpsi_ctResolution2");
  RooRealVar *fracGauss2_=(RooRealVar*)ws->var("jpsi_fracGauss2");

  RooRealVar *nonPromptTau_=(RooRealVar*)ws->var("jpsi_nonPromptTau");

  char varchar[200];

  if(SB==1) sprintf(varchar,"jpsi_fBkgLSB");
  if(SB==2) sprintf(varchar,"jpsi_fBkgRSB");
  RooRealVar *fBkg_ = (RooRealVar*)ws->var(varchar);
  if(SB==1) sprintf(varchar,"jpsi_fPromptLSB");
  if(SB==2) sprintf(varchar,"jpsi_fPromptRSB");
  RooRealVar *fPrompt_ = (RooRealVar*)ws->function(varchar);



  double promptMean = promptMean_->getVal();
  double promptMeanErr = promptMean_->getError();
  double promptCtRe = ctResolution_->getVal();
  double promptCtReErr = ctResolution_->getError();
  double promptCtRe2 = ctResolution2_->getVal();
  double promptCtRe2Err = ctResolution2_->getError();
  double fracGauss2 = fracGauss2_->getVal();
  double fracGauss2Err = fracGauss2_->getError();
  double nonPromptTau = nonPromptTau_->getVal();
  double nonPromptTauErr = nonPromptTau_->getError();

  double fBkg = fBkg_->getVal();
  double fBkgErr = fBkg_->getError();
  double fPrompt = fPrompt_->getVal();
  double fPromptErr = fPrompt_->getError();
  double fNonPrompt =  1.-fBkg-fPrompt;

  cout<<"promptMean: "<<promptMean<<endl;
  cout<<"promptCtRe: "<<promptCtRe<<endl;
  cout<<"nonPromptTau: "<<nonPromptTau<<endl;
  cout<<"fBkg: "<<fBkg<<endl;
  cout<<"fBkgErr: "<<fBkgErr<<endl;
  cout<<"fPrompt: "<<fPrompt<<endl;
  cout<<"fNonPrompt: "<<fNonPrompt<<endl;

  char modelchar[200];
  if(SB==1) sprintf(modelchar,"jpsi_backgroundlifetimeLnoC");
  if(SB==2) sprintf(modelchar,"jpsi_backgroundlifetimeRnoC");

  RooAddPdf *ModelLife = (RooAddPdf*)ws->pdf(modelchar);

  RooAddModel *Prompt = (RooAddModel*)ws->pdf("jpsi_TotalPromptLifetime");
  RooAbsPdf *nonPromptSSD = (RooAbsPdf*)ws->pdf("jpsi_nonPromptSSD");

  if(SB==1) sprintf(modelchar,"jpsi_backgroundlifetimeLpre");
  if(SB==2) sprintf(modelchar,"jpsi_backgroundlifetimeRpre");
  RooAbsPdf *backgroundlifetime = (RooAbsPdf*)ws->pdf(modelchar);

  int parsFit;
  //ploting signal region
  dataSB->plotOn(ctauFrame,MarkerSize(onia::markerSize_ML), Name("myHist"));

  maxY = ctauFrame->GetMaximum()*enlargeYby;
  minY = 5e-1;
  ctauFrame->GetYaxis()->SetRangeUser(minY,maxY);

  int JpsiFitColor=633;

  ModelLife->plotOn(ctauFrame,
                    ProjWData(*JpsictErr, *dataSBProj),
                    LineWidth(2),
                    LineColor(JpsiFitColor),
                    NumCPU(1),
                    Name("myCurve"));

  //------get chi2------------SHOULD DONE after PLOTTING------
  parsFit=(fitRlt->floatParsFinal()).getSize(); //this used the full p.d.f.

  int nBins_LSig=ctauFrame->GetNbinsX();
  double chi2Pre_LSig=ctauFrame->chiSquare(parsFit);  //reduced chi-squared = chi2/ndof
  int ndof_LSig=nBins_LSig-parsFit;  //num of degree of freedom
  double chi2_LSig=chi2Pre_LSig*ndof_LSig;

  //RooHist* hpull_ctauSig = ctauFrame->pullHist() ;
  RooHist* hpull_ctauSig = ctauFrame->pullHist("myHist","myCurve",kTRUE);

  //save pull distribution
  TH1F* pull = new TH1F("pull","pull distribution", 100,-10.,10.);
  gSystem->mkdir("Fit/pull",kTRUE);gSystem->mkdir("Fit/jpsiFit",kTRUE);
  //TFile *pullFile = new TFile(Form("Fit/pull/pullAverage_rap%d_pt%d_SB.root",rapBin,ptBin),"RECREATE");
  TFile *pullFile = new TFile(Form("Fit/pull/pull_rap%d_pt%d_%s.root",rapBin,ptBin,SBchar),"RECREATE");

  hpull_ctauSig->SetMarkerSize(0.8);
  for(int i=0;i<hpull_ctauSig->GetN();i++){
    hpull_ctauSig->SetPointEYlow(i,0.);
    hpull_ctauSig->SetPointEYhigh(i,0.);
    double x,y;
    hpull_ctauSig->GetPoint(i,x,y);
    pull->Fill(y);
  }
  pullFile->cd();
  pull->Write();
  pullFile->Close();

  ctauFramePull->addPlotable(hpull_ctauSig,"P");

  //Prompt->plotOn(ctauFrame,
  //		ProjWData(*JpsictErr, *dataSBProj),
  //		Normalization(fPrompt),
  //		LineStyle(5),
  //		LineColor(onia::ColorPRJpsi),
  //		LineWidth(2), NumCPU(1));
  //
  //nonPromptSSD->plotOn(ctauFrame,
  //		ProjWData(*JpsictErr, *dataSBProj),
  //		Normalization(fNonPrompt),
  //		LineStyle(2),
  //		LineColor(onia::ColorNPJpsi),
  //		LineWidth(2), NumCPU(1));
  //
  //backgroundlifetime->plotOn(ctauFrame,
  //		ProjWData(*JpsictErr, *dataSBProj),
  //		Normalization(fBkg),
  //		LineStyle(7),
  //		LineColor(onia::ColorMuMuBG),
  //		LineWidth(2), NumCPU(1));

  ModelLife->plotOn(ctauFrame,
                    Components("jpsi_TotalPromptLifetime"),
                    ProjWData(*JpsictErr, *dataSBProj),
                    LineStyle(5),
                    LineColor(onia::ColorPRJpsi),
                    LineWidth(2), NumCPU(1));

  ModelLife->plotOn(ctauFrame,
                    Components("jpsi_nonPromptSSD"),
                    ProjWData(*JpsictErr, *dataSBProj),
                    LineStyle(2),
                    LineColor(onia::ColorNPJpsi),
                    LineWidth(2), NumCPU(1));

  if(SB==1) sprintf(modelchar,"jpsi_backgroundlifetimeLpre");
  if(SB==2) sprintf(modelchar,"jpsi_backgroundlifetimeRpre");
  ModelLife->plotOn(ctauFrame,
                    Components(modelchar),
                    ProjWData(*JpsictErr, *dataSBProj),
                    LineStyle(7),
                    LineColor(onia::ColorMuMuBG),
                    LineWidth(2), NumCPU(1));


  double Ymax = ctauFrame->GetMaximum();
  cout<<"Ymax: "<<Ymax<<endl;
  ctauFrame->SetMaximum(3*Ymax);
  ctauFrame->SetMinimum(1.1);


  TH1* legendBlue = dataSB->createHistogram("legendBlue",JpsiMass,Binning(50)) ; legendBlue->SetLineColor(JpsiFitColor) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
  TH1* legendBlueDash = dataSB->createHistogram("legendBlueDash",JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(onia::ColorPRJpsi) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
  TH1* legendRed = dataSB->createHistogram("legendRed",JpsiMass,Binning(50)) ; legendRed->SetLineColor(onia::ColorNPJpsi) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
  TH1* legendBlack = dataSB->createHistogram("legendBlack",JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
  TH1* legendGreen = dataSB->createHistogram("legendGreen",JpsiMass,Binning(50)) ; legendGreen->SetLineColor(kGreen) ; legendGreen->SetLineStyle(kSolid) ; legendGreen->SetLineWidth(2) ;
  TH1* legendGreenDash = dataSB->createHistogram("legendGreenDash",JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(kGreen) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
  TH1* legendPink = dataSB->createHistogram("legendPink",JpsiMass,Binning(50)) ; legendPink->SetLineColor(onia::ColorMuMuBG) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;


  TLegend* LifetimeLegendSig=new TLegend(0.13,0.75,0.23,0.91);
  LifetimeLegendSig->SetFillColor(kWhite);
  LifetimeLegendSig->SetTextFont(42);
  LifetimeLegendSig->SetTextSize(0.035);
  LifetimeLegendSig->SetBorderSize(0.);
  LifetimeLegendSig->AddEntry(legendBlue,"sum","l");
  LifetimeLegendSig->AddEntry(legendBlueDash,"PR J/#psi","l");
  LifetimeLegendSig->AddEntry(legendRed,"NP J/#psi","l");
  LifetimeLegendSig->AddEntry(legendPink,"#mu#mu BG","l");

  double left=0.5, top=0.9, textSize=0.03;
  TLatex *latex=new TLatex();
  latex->SetTextFont(42);
  latex->SetNDC(kTRUE);
  latex->SetTextSize(textSize);
  double step=textSize*1.6;

  gStyle->SetPadBottomMargin(0.08); //0.12
  gStyle->SetPadLeftMargin(0.09); //0.12
  gStyle->SetPadRightMargin(0.035); //0.05
  gStyle->SetPadTopMargin(0.05); //0.05


  TCanvas *c1;
  TPad *pad1;
  TPad *pad2;

  ctauFrame->SetMinimum(minY);
  ctauFrame->SetMaximum(maxY);

  for(int LinLog=0; LinLog<2; LinLog++){
    cout<<"LinLog "<<LinLog<<endl;

    c1=new TCanvas("c1","",1000,900);

    c1->cd();
    pad1 = new TPad("pad1","pad1",0.,0.,1.,0.3);
    pad1->SetGridy();
    pad1->SetBottomMargin(0.2);
    pad1->Draw();
    c1->cd();
    pad2 = new TPad("pad2","pad2",0.,0.3,1.,1.);
    pad2->Draw();

    //Sig
    pad2->cd(0);
    if(LinLog==0) pad2->SetLogy(false);
    if(LinLog==1) pad2->SetLogy(true);

    ctauFrame->Draw();
    LifetimeLegendSig->Draw();

    left=0.54; top=0.885; textSize=0.030; latex->SetTextSize(textSize);
    if(rapBin==0) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[onia::kNbRapForPTBins]));
    else if(rapBin==1) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[rapBin]));
    else latex->DrawLatex(left,top,Form("%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin-1],onia::KinParticleChar,onia::rapForPTRange[rapBin]));
    top-=step;
    if(ptBin==0)
      latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTBins[rapBin]]));
    else
      latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]));

    top-=step;
    latex->SetTextColor(kRed);
    latex->DrawLatex(left,top,Form("J/#psi %s",SBchar));
    latex->SetTextColor(kBlack);
    top-=step;
    latex->DrawLatex(left,top,Form("#chi^{2}/ndf = %.1f/%d",chi2_LSig,ndof_LSig));
    top-=step;


    latexFloatingLifetimePars(ws, latex, SB);


    pad1->cd(0); pad1->SetLogy(0);
    ctauFramePull->Draw();

    std::stringstream saveCtau;
    if(!zoom){
      if(LinLog==0) saveCtau << "Fit/jpsiFit/ctau_" << SBchar << "_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
      if(LinLog==1) saveCtau << "Fit/jpsiFit/ctau_" << SBchar << "_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
    }
    else{
      if(LinLog==0) saveCtau << "Fit/jpsiFit/ctau_" << SBchar << "_zoom_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
      if(LinLog==1) saveCtau << "Fit/jpsiFit/ctau_" << SBchar << "_zoom_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
    }

    c1->SaveAs(saveCtau.str().c_str());

  }

  if(correctResolutionForPlotting){
    ws->var("jpsi_ctResolution")->setVal(ws->var("jpsi_ctResolution")->getVal()/resCorrFactor);
    ws->var("jpsi_ctResolution2")->setVal(ws->var("jpsi_ctResolution2")->getVal()/resCorrFactor);
  }

  delete c1;
  delete legendBlue;
  delete legendBlueDash;
  delete legendRed;
  delete legendBlack;
  delete legendGreen;
  delete legendGreenDash;
  delete legendPink;
  return chi2_LSig/ndof_LSig;
}


//==============================================
double plotJpsiPedagogical(RooWorkspace *ws, int rapBin, int ptBin, int nState, bool SpeedPlotting, int nSpeedPlotting, bool zoom, int LinLog){

  bool plotModels=true;

  int  nbins=60; //0.005 bin size
  TGaxis::SetMaxDigits(3);

  double PlotMin;
  double PlotMax;

  PlotMin=2.875;
  PlotMax=3.3;


  gSystem->mkdir("Fit/jpsiFit",kTRUE);

  double binWidth=(PlotMax-PlotMin)/double(nbins)*1000;

  RooRealVar *JpsiMass = ws->var("JpsiMass");
  assert( 0 != JpsiMass );
  RooRealVar *JpsiRap = ws->var("JpsiRap");
  assert( 0 != JpsiRap );

  RooPlot *massFrame = JpsiMass->frame(Bins(nbins), Range(PlotMin, PlotMax));
  assert ( 0 != massFrame );
  massFrame->SetName(Form("mass_plot_rap%d_pt%d",rapBin,ptBin));
  massFrame->SetTitle("");
  massFrame->GetYaxis()->SetTitle(Form("Events / %1.0f MeV",binWidth));
  massFrame->GetYaxis()->SetTitleOffset(0.7);


  RooAbsData *data= ws->data(Form("jpsi_data_rap%d_pt%d",rapBin,ptBin));
  assert ( 0 != data );
  RooFitResult* fitRlt = dynamic_cast<RooFitResult*>(ws->obj(Form("m_fitresult_rap%d_pt%d",rapBin,ptBin)));
  assert ( 0 != fitRlt);


  double Mean = -1.0, MeanErr = -1.0;
  getVarFromWorkspace(ws, "CBmass_p0_jpsi", Mean, MeanErr);
  double Sigma = -1.0, SigmaErr = -1.0;
  getVarFromWorkspace(ws, "var_massres", Sigma, SigmaErr);
  double fracBkg = -1.0, fracBkgErr = -1.0;
  getVarFromWorkspace(ws, "fracBkg_jpsi", fracBkg, fracBkgErr);

  std::stringstream masssnapshotname;
  masssnapshotname << "m_snapshot_rap" << rapBin << "_pt" << ptBin;
  ws->loadSnapshot(masssnapshotname.str().c_str());

  RooAbsPdf *massPdf = ws->pdf("massModel_jpsi");
  assert ( 0 != massPdf );
  RooAbsPdf *bkgMassShape = ws->pdf("bkgMassShape_jpsi");
  assert ( 0 != bkgMassShape );
  RooAbsPdf *sigMassShape_jpsi = ws->pdf("sigMassShape_jpsi");
  assert ( 0 != sigMassShape_jpsi );
  RooAbsPdf *gaussMassShape_jpsi = ws->pdf("gaussMassShape_jpsi");
  assert ( 0 != gaussMassShape_jpsi );
  RooAbsPdf *FullmassPdf = ws->pdf("FullmassPdf");
  assert ( 0 != FullmassPdf );

  int nEntries = data->numEntries();

  int JpsiFitColor=633;

  cout<<"plot data"<<endl;
  data->plotOn(massFrame,MarkerSize(onia::markerSize_ML));
  cout<<"plot FullmassPdf"<<endl;

  if(plotModels){
    FullmassPdf->plotOn(massFrame,
			Normalization(nEntries,2),
			LineWidth(2),
			LineColor(JpsiFitColor)
			);


    cout<<"plot sigMassShape_jpsi"<<endl;
    sigMassShape_jpsi->plotOn(massFrame,
                              Normalization(nEntries*(1.-fracBkg),2),
                              LineStyle(2),
                              LineColor(onia::ColorSumJpsiSignal),
                              LineWidth(2),
                              ProjWData(*JpsiRap, *data));
    cout<<"plot bkgMassShape"<<endl;
    bkgMassShape->plotOn(massFrame,
                         Normalization(nEntries*fracBkg,2),
                         LineStyle(2),
                         LineColor(onia::ColorMuMuBG),
                         LineWidth(2),
                         ProjWData(*JpsiRap, *data));
  }

  double minY = 0.;

  double sbLowMass=Mean-onia::nSigBkgLow*Sigma;
  double sbHighMass=Mean+onia::nSigBkgHigh*Sigma;
  double sigMinMass=Mean-onia::nSigMass*Sigma;
  double sigMaxMass=Mean+onia::nSigMass*Sigma;

  double maxY = 0.;
  maxY = massFrame->GetMaximum()*0.3;
  double lineWidth = 2.0;
  TLine *lineSBLow = new TLine(sbLowMass, minY, sbLowMass, maxY);
  TLine *lineSBHigh = new TLine(sbHighMass, minY, sbHighMass, maxY);
  TLine *lineSigLow = new TLine(sigMinMass, minY, sigMinMass, maxY);
  TLine *lineSigHigh = new TLine(sigMaxMass, minY, sigMaxMass, maxY);
  lineSBLow->SetLineWidth(lineWidth);lineSBHigh->SetLineWidth(lineWidth);
  lineSigLow->SetLineWidth(lineWidth);lineSigHigh->SetLineWidth(lineWidth);
  lineSBLow->SetLineColor(onia::ColorMuMuBG);lineSBHigh->SetLineColor(onia::ColorMuMuBG);
  lineSigLow->SetLineColor(onia::ColorSumJpsiSignal);lineSigHigh->SetLineColor(onia::ColorSumJpsiSignal);
  lineSBLow->SetLineStyle(7);lineSBHigh->SetLineStyle(7);
  lineSigLow->SetLineStyle(5);lineSigHigh->SetLineStyle(5);

  TH1* legendBlue = data->createHistogram("legendBlue",*JpsiMass,Binning(50)) ; legendBlue->SetLineColor(JpsiFitColor) ; legendBlue->SetLineStyle(kSolid) ; legendBlue->SetLineWidth(2) ;
  TH1* legendBlueDash = data->createHistogram("legendBlueDash",*JpsiMass,Binning(50)) ; legendBlueDash->SetLineColor(onia::ColorPRJpsi) ; legendBlueDash->SetLineStyle(5) ; legendBlueDash->SetLineWidth(2) ;
  TH1* legendRed = data->createHistogram("legendRed",*JpsiMass,Binning(50)) ; legendRed->SetLineColor(onia::ColorNPJpsi) ; legendRed->SetLineStyle(2) ; legendRed->SetLineWidth(2) ;
  TH1* legendBlack = data->createHistogram("legendBlack",*JpsiMass,Binning(50)) ; legendBlack->SetLineColor(kBlack) ; legendBlack->SetLineStyle(kSolid) ; legendBlack->SetLineWidth(2) ;
  TH1* legendGreen = data->createHistogram("legendGreen",*JpsiMass,Binning(50)) ; legendGreen->SetLineColor(onia::ColorSumJpsiSignal) ; legendGreen->SetLineStyle(2) ; legendGreen->SetLineWidth(2) ;
  TH1* legendGreenDash = data->createHistogram("legendGreenDash",*JpsiMass,Binning(50)) ; legendGreenDash->SetLineColor(onia::ColorSumJpsiSignal) ; legendGreenDash->SetLineStyle(2) ; legendGreenDash->SetLineWidth(2) ;
  TH1* legendPink = data->createHistogram("legendPink",*JpsiMass,Binning(50)) ; legendPink->SetLineColor(onia::ColorMuMuBG) ; legendPink->SetLineStyle(7) ; legendPink->SetLineWidth(2) ;

  TLegend* MassLegend=new TLegend(0.1,0.585,0.3,0.925);
  MassLegend->SetFillColor(kWhite);
  MassLegend->SetTextFont(42);
  MassLegend->SetTextSize(0.055);
  MassLegend->SetBorderSize(0.);
  MassLegend->AddEntry(legendBlue,"sum","l");
  MassLegend->AddEntry(legendGreen,"PR+NP J/#psi","l");
  MassLegend->AddEntry(legendBlueDash,"PR J/#psi","l");
  MassLegend->AddEntry(legendRed,"NP J/#psi","l");
  MassLegend->AddEntry(legendPink,"#mu#mu BG","l");














  nbins=onia::LifetimePlotBins/3.;
  TGaxis::SetMaxDigits(3);


  if(zoom){
    PlotMin=onia::ctPlotMinZoom;
    PlotMax=onia::ctPlotMaxZoom;
  }
  else{
    PlotMin=onia::ctPlotMin;
    PlotMax=onia::ctPlotMax;
  }

  binWidth=(PlotMax-PlotMin)/double(nbins)*1000;

  RooRealVar Jpsict(*ws->var("Jpsict"));
  RooRealVar *JpsictErr = ws->var("JpsictErr");

  double enlargeYby=onia::enlargeYby_ML;

  RooFitResult* fitRltCtau = (RooFitResult*)ws->obj(Form("jpsi_l_fitresult_rap%d_pt%d",rapBin,ptBin));
  if(!fitRltCtau) { cout<<">>=======Error: no Fit result object in workspace=========="<<endl; return 0.; }
  fitRltCtau->Print();

  RooAddModel *Prompt = (RooAddModel*)ws->pdf("jpsi_TotalPromptLifetime");
  RooAbsPdf *nonPromptSSD = (RooAbsPdf*)ws->pdf("jpsi_nonPromptSSD");

  RooAddPdf *ModelLifeLSB = (RooAddPdf*)ws->pdf("jpsi_backgroundlifetimeLnoC");
  RooAddPdf *ModelLifeRSB = (RooAddPdf*)ws->pdf("jpsi_backgroundlifetimeRnoC");
  RooAddPdf *ModelLifeSR = (RooAddPdf*)ws->pdf("jpsi_fulllifetimeSR");

  RooAbsPdf *backgroundlifetimeLSB = (RooAbsPdf*)ws->pdf("jpsi_backgroundlifetimeLpre");
  RooAbsPdf *backgroundlifetimeRSB = (RooAbsPdf*)ws->pdf("jpsi_backgroundlifetimeLpre");
  RooAbsPdf *backgroundlifetimeSR = (RooAbsPdf*)ws->pdf("jpsi_backgroundlifetime");




  RooRealVar *fBkgLSB_ = (RooRealVar*)ws->var("jpsi_fBkgLSB");
  RooRealVar *fBkgRSB_ = (RooRealVar*)ws->var("jpsi_fBkgRSB");
  RooRealVar *fBkgSR_ = (RooRealVar*)ws->var("jpsi_fBkg");
  RooRealVar *fPromptLSB_ = (RooRealVar*)ws->function("jpsi_fPromptLSB");
  RooRealVar *fPromptRSB_ = (RooRealVar*)ws->function("jpsi_fPromptRSB");
  RooRealVar *fPromptSR_ = (RooRealVar*)ws->function("jpsi_fPrompt");

  double fBkgLSB = fBkgLSB_->getVal();
  double fPromptLSB = fPromptLSB_->getVal();
  double fNonPromptLSB =  1.-fBkgLSB-fPromptLSB;

  double fBkgRSB = fBkgRSB_->getVal();
  double fPromptRSB = fPromptRSB_->getVal();
  double fNonPromptRSB =  1.-fBkgRSB-fPromptRSB;

  double fBkgSR = fBkgSR_->getVal();
  double fPromptSR = fPromptSR_->getVal();
  double fNonPromptSR =  1.-fBkgSR-fPromptSR;

  double Ymax=0;


  double ctauYOffset=1.75;

  RooPlot *ctauFrameLSB=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(PlotMin, PlotMax));
  ctauFrameLSB->SetName(Form("ctausig_plot_rap%d_pt%d",rapBin,ptBin));
  ctauFrameLSB->GetYaxis()->SetTitle(Form("Events / %1.0f micron",binWidth));
  ctauFrameLSB->SetTitle("");
  ctauFrameLSB->GetYaxis()->SetTitleOffset(ctauYOffset);

  RooDataSet *dataLSB=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_LSB",rapBin,ptBin));
  RooAbsData* dataLSBProj;
  if(SpeedPlotting){
    dataLSBProj = dataLSB->reduce(EventRange(0,nSpeedPlotting));
  }
  else{
    dataLSBProj = dataLSB->reduce(EventRange(0,1e10));
  }
  cout<<"number of events in dataLSB = "<<dataLSB->numEntries()<<endl;
  cout<<"number of events in dataLSBProj = "<<dataLSBProj->numEntries()<<endl;

  cout<<"dataLSB->numEntries: "<<dataLSB->numEntries()<<endl;


  //ploting signal region
  dataLSB->plotOn(ctauFrameLSB,MarkerSize(onia::markerSize_ML), Name("myHist"));

  maxY = ctauFrameLSB->GetMaximum()*enlargeYby;
  minY = 5e-1;
  ctauFrameLSB->GetYaxis()->SetRangeUser(minY,maxY);

  if(plotModels){
    ModelLifeLSB->plotOn(ctauFrameLSB,
                         ProjWData(*JpsictErr, *dataLSBProj),
                         LineWidth(2),
                         LineColor(JpsiFitColor),
                         NumCPU(1),
                         Name("myCurve"));


    ModelLifeLSB->plotOn(ctauFrameLSB,
                         Components("jpsi_TotalPromptLifetime"),
                         ProjWData(*JpsictErr, *dataLSBProj),
                         LineStyle(5),
                         LineColor(onia::ColorPRJpsi),
                         LineWidth(2), NumCPU(1));

    ModelLifeLSB->plotOn(ctauFrameLSB,
                         Components("jpsi_nonPromptSSD"),
                         ProjWData(*JpsictErr, *dataLSBProj),
                         LineStyle(2),
                         LineColor(onia::ColorNPJpsi),
                         LineWidth(2), NumCPU(1));

    ModelLifeLSB->plotOn(ctauFrameLSB,
                         Components("jpsi_backgroundlifetimeLpre"),
                         ProjWData(*JpsictErr, *dataLSBProj),
                         LineStyle(7),
                         LineColor(onia::ColorMuMuBG),
                         LineWidth(2), NumCPU(1));
  }

  Ymax = ctauFrameLSB->GetMaximum();
  cout<<"Ymax: "<<Ymax<<endl;
  ctauFrameLSB->SetMaximum(Ymax);
  ctauFrameLSB->SetMinimum(1.1);




  RooPlot *ctauFrameRSB=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(PlotMin, PlotMax));
  ctauFrameRSB->SetName(Form("ctausig_plot_rap%d_pt%d",rapBin,ptBin));
  ctauFrameRSB->GetYaxis()->SetTitle(Form("Events / %1.0f micron",binWidth));
  ctauFrameRSB->SetTitle("");
  ctauFrameRSB->GetYaxis()->SetTitleOffset(ctauYOffset);

  RooDataSet *dataRSB=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_RSB",rapBin,ptBin));
  RooAbsData* dataRSBProj;
  if(SpeedPlotting){
    dataRSBProj = dataRSB->reduce(EventRange(0,nSpeedPlotting));
  }
  else{
    dataRSBProj = dataRSB->reduce(EventRange(0,1e10));
  }
  cout<<"number of events in dataRSB = "<<dataRSB->numEntries()<<endl;
  cout<<"number of events in dataRSBProj = "<<dataRSBProj->numEntries()<<endl;

  cout<<"dataRSB->numEntries: "<<dataRSB->numEntries()<<endl;


  //ploting signal region
  dataRSB->plotOn(ctauFrameRSB,MarkerSize(onia::markerSize_ML), Name("myHist"));

  maxY = ctauFrameRSB->GetMaximum()*enlargeYby;
  minY = 5e-1;
  ctauFrameRSB->GetYaxis()->SetRangeUser(minY,maxY);

  if(plotModels){
    ModelLifeRSB->plotOn(ctauFrameRSB,
                         ProjWData(*JpsictErr, *dataRSBProj),
                         LineWidth(2),
                         LineColor(JpsiFitColor),
                         NumCPU(1),
                         Name("myCurve"));


    ModelLifeRSB->plotOn(ctauFrameRSB,
                         Components("jpsi_TotalPromptLifetime"),
                         ProjWData(*JpsictErr, *dataRSBProj),
                         LineStyle(5),
                         LineColor(onia::ColorPRJpsi),
                         LineWidth(2), NumCPU(1));

    ModelLifeRSB->plotOn(ctauFrameRSB,
                         Components("jpsi_nonPromptSSD"),
                         ProjWData(*JpsictErr, *dataRSBProj),
                         LineStyle(2),
                         LineColor(onia::ColorNPJpsi),
                         LineWidth(2), NumCPU(1));

    ModelLifeRSB->plotOn(ctauFrameRSB,
                         Components("jpsi_backgroundlifetimeRpre"),
                         ProjWData(*JpsictErr, *dataRSBProj),
                         LineStyle(7),
                         LineColor(onia::ColorMuMuBG),
                         LineWidth(2), NumCPU(1));
  }

  Ymax = ctauFrameRSB->GetMaximum();
  cout<<"Ymax: "<<Ymax<<endl;
  ctauFrameRSB->SetMaximum(Ymax);
  ctauFrameRSB->SetMinimum(1.1);




  RooPlot *ctauFrameSR=((RooRealVar*)ws->var("Jpsict"))->frame(Bins(nbins), Range(PlotMin, PlotMax));
  ctauFrameSR->SetName(Form("ctausig_plot_rap%d_pt%d",rapBin,ptBin));
  ctauFrameSR->GetYaxis()->SetTitle(Form("Events / %1.0f micron",binWidth));
  ctauFrameSR->SetTitle("");
  ctauFrameSR->GetYaxis()->SetTitleOffset(ctauYOffset);

  RooDataSet *dataSR=(RooDataSet *)ws->data(Form("data_rap%d_pt%d_SR",rapBin,ptBin));
  RooAbsData* dataSRProj;
  if(SpeedPlotting){
    dataSRProj = dataSR->reduce(EventRange(0,nSpeedPlotting));
  }
  else{
    dataSRProj = dataSR->reduce(EventRange(0,1e10));
  }
  cout<<"number of events in dataSR = "<<dataSR->numEntries()<<endl;
  cout<<"number of events in dataSRProj = "<<dataSRProj->numEntries()<<endl;

  cout<<"dataSR->numEntries: "<<dataSR->numEntries()<<endl;


  //ploting signal region
  dataSR->plotOn(ctauFrameSR,MarkerSize(onia::markerSize_ML), Name("myHist"));

  maxY = ctauFrameSR->GetMaximum()*enlargeYby;
  minY = 5e-1;
  ctauFrameSR->GetYaxis()->SetRangeUser(minY,maxY);

  if(plotModels){
    ModelLifeSR->plotOn(ctauFrameSR,
			ProjWData(*JpsictErr, *dataSRProj),
			LineWidth(2),
			LineColor(JpsiFitColor),
			NumCPU(1),
			Name("myCurve"));


    ModelLifeSR->plotOn(ctauFrameSR,
			Components("jpsi_TotalPromptLifetime"),
			ProjWData(*JpsictErr, *dataSRProj),
			LineStyle(5),
			LineColor(onia::ColorPRJpsi),
			LineWidth(2), NumCPU(1));

    ModelLifeSR->plotOn(ctauFrameSR,
			Components("jpsi_nonPromptSSD"),
			ProjWData(*JpsictErr, *dataSRProj),
			LineStyle(2),
			LineColor(onia::ColorNPJpsi),
			LineWidth(2), NumCPU(1));

    ModelLifeSR->plotOn(ctauFrameSR,
			Components("jpsi_backgroundlifetime"),
			ProjWData(*JpsictErr, *dataSRProj),
			LineStyle(7),
			LineColor(onia::ColorMuMuBG),
			LineWidth(2), NumCPU(1));
  }

  Ymax = ctauFrameSR->GetMaximum();
  cout<<"Ymax: "<<Ymax<<endl;
  ctauFrameSR->SetMaximum(Ymax);
  ctauFrameSR->SetMinimum(1.1);

















  if(LinLog==1) ctauFrameSR->SetMinimum(0.5);
  if(LinLog==1) ctauFrameLSB->SetMinimum(0.5);
  if(LinLog==1) ctauFrameRSB->SetMinimum(0.5);
  if(LinLog==1) ctauFrameSR->SetMaximum(ctauFrameSR->GetMaximum()*3.);
  if(LinLog==1) ctauFrameLSB->SetMaximum(ctauFrameLSB->GetMaximum()*3.);
  if(LinLog==1) ctauFrameRSB->SetMaximum(ctauFrameRSB->GetMaximum()*3.);

  if(LinLog==1) massFrame->SetMinimum(massFrame->GetMaximum()*4e-3);
  if(LinLog==1) massFrame->SetMaximum(massFrame->GetMaximum()*3.);

  double left=0.5, top=0.9, textSize=0.03;
  TLatex *latex=new TLatex();
  latex->SetTextFont(42);
  latex->SetNDC(kTRUE);
  latex->SetTextSize(textSize);
  double step=textSize*1.6;


  TCanvas *c1;
  TPad *pad1;
  TPad *pad2;
  TPad *pad3;
  TPad *pad4;


  c1=new TCanvas("c1","",1200,900);

  c1->cd();
  pad1 = new TPad("pad1","pad1",0.,0.5,1.,1.);
  pad1->Draw();
  c1->cd();
  pad2 = new TPad("pad2","pad2",0.,0.,1./3.,0.5);
  pad2->Draw();
  c1->cd();
  pad3 = new TPad("pad3","pad3",1./3.,0.,2./3.,0.5);
  pad3->Draw();
  c1->cd();
  pad4 = new TPad("pad4","pad4",2./3.,0.,3./3.,0.5);
  pad4->Draw();

  double bottomMargin=0.08;
  double leftMargin=0.05;
  double rightMargin=0.01;
  double topMargin=0.05;

  pad1->SetBottomMargin(bottomMargin);
  pad1->SetLeftMargin(leftMargin);
  pad1->SetRightMargin(rightMargin);
  pad1->SetTopMargin(topMargin);

  pad2->SetBottomMargin(bottomMargin);
  pad2->SetLeftMargin(leftMargin*3.);
  pad2->SetRightMargin(rightMargin*3.);
  pad2->SetTopMargin(topMargin);

  pad3->SetBottomMargin(bottomMargin);
  pad3->SetLeftMargin(leftMargin*3.);
  pad3->SetRightMargin(rightMargin*3.);
  pad3->SetTopMargin(topMargin);

  pad4->SetBottomMargin(bottomMargin);
  pad4->SetLeftMargin(leftMargin*3.);
  pad4->SetRightMargin(rightMargin*3.);
  pad4->SetTopMargin(topMargin);


  pad1->cd(0);
  if(LinLog==0) pad1->SetLogy(false);
  if(LinLog==1) pad1->SetLogy(true);
  pad2->cd(0);
  if(LinLog==0) pad2->SetLogy(false);
  if(LinLog==1) pad2->SetLogy(true);
  pad3->cd(0);
  if(LinLog==0) pad3->SetLogy(false);
  if(LinLog==1) pad3->SetLogy(true);
  pad4->cd(0);
  if(LinLog==0) pad4->SetLogy(false);
  if(LinLog==1) pad4->SetLogy(true);

  pad1->cd(0);
  massFrame->Draw();
  MassLegend->Draw();
  lineSBLow->Draw("same");
  lineSBHigh->Draw("same");
  lineSigLow->Draw("same");
  lineSigHigh->Draw("same");

  left=0.8; top=0.855; textSize=0.060; latex->SetTextSize(textSize);
  if(rapBin==0) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[onia::kNbRapForPTBins]));
  else if(rapBin==1) latex->DrawLatex(left,top,Form("|y%s| < %.1f",onia::KinParticleChar,onia::rapForPTRange[rapBin]));
  else latex->DrawLatex(left,top,Form("%.1f < |y%s| < %.1f",onia::rapForPTRange[rapBin-1],onia::KinParticleChar,onia::rapForPTRange[rapBin]));
  step=textSize*1.6;
  top-=step;
  if(ptBin==0)
    latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin],onia::KinParticleChar,onia::pTRange[rapBin][onia::kNbPTBins[rapBin]]));
  else
    latex->DrawLatex(left,top,Form("%.0f < p%s_{T} < %.0f GeV",onia::pTRange[rapBin][ptBin-1],onia::KinParticleChar,onia::pTRange[rapBin][ptBin]));


  pad2->cd(0);
  ctauFrameLSB->Draw();

  left=0.7; top=0.855; textSize=0.065; latex->SetTextSize(textSize);
  latex->SetTextColor(kRed);
  latex->DrawLatex(left,top,Form("J/#psi LSB"));

  pad3->cd(0);
  ctauFrameSR->Draw();

  left=0.7; top=0.855; textSize=0.065; latex->SetTextSize(textSize);
  latex->SetTextColor(kRed);
  latex->DrawLatex(left,top,Form("J/#psi SR"));

  pad4->cd(0);
  ctauFrameRSB->Draw();

  left=0.7; top=0.855; textSize=0.065; latex->SetTextSize(textSize);
  latex->SetTextColor(kRed);
  latex->DrawLatex(left,top,Form("J/#psi RSB"));


  std::stringstream saveCtau;
  if(!zoom){
    if(LinLog==0) saveCtau << "Fit/jpsiFit/pedagogical_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
    if(LinLog==1) saveCtau << "Fit/jpsiFit/pedagogical_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
  }
  else{
    if(LinLog==0) saveCtau << "Fit/jpsiFit/pedagogical_zoom_lin_rap" << rapBin << "_pt" << ptBin << ".pdf";
    if(LinLog==1) saveCtau << "Fit/jpsiFit/pedagogical_zoom_log_rap" << rapBin << "_pt" << ptBin << ".pdf";
  }

  c1->SaveAs(saveCtau.str().c_str());


  delete c1;
  delete legendBlue;
  delete legendBlueDash;
  delete legendRed;
  delete legendBlack;
  delete legendGreen;
  delete legendGreenDash;
  delete legendPink;
}

void latexFloatingLifetimePars(RooWorkspace *ws, TLatex* latex, int region){


  double textSize, left, top;

  textSize=0.03; latex->SetTextSize(textSize);
  left=0.725; top=0.885;
  double stepsizeTimes=1.9;
  latex->SetTextSize(textSize);
  if(region==0){


    if(!ws->var("jpsi_promptMean")->isConstant()){
      latex->DrawLatex(left,top,Form("#mu_{c#tau}  =  %.3f #pm %.3f",ws->var("jpsi_promptMean")->getVal(), ws->var("jpsi_promptMean")->getError()));
      top-=textSize*stepsizeTimes;
    }
    if(!ws->var("jpsi_ctResolution")->isConstant()){
      latex->DrawLatex(left,top,Form("#sigma^{scale}_{l}  =  %.3f #pm %.3f",ws->var("jpsi_ctResolution")->getVal(), ws->var("jpsi_ctResolution")->getError()));
      top-=textSize*stepsizeTimes;
    }
    if(!ws->var("jpsi_ctResolution2")->isConstant()){
      latex->DrawLatex(left,top,Form("#sigma^{scale2}_{l}  =  %.3f #pm %.3f",ws->var("jpsi_ctResolution2")->getVal(), ws->var("jpsi_ctResolution2")->getError()));
      top-=textSize*stepsizeTimes;
    }
    if(!ws->var("jpsi_fracGauss2")->isConstant()){
      latex->DrawLatex(left,top,Form("f_{G_{2}}  =  %.3f #pm %.3f",ws->var("jpsi_fracGauss2")->getVal(), ws->var("jpsi_fracGauss2")->getError()));
      top-=textSize*stepsizeTimes;
    }
    if(!ws->var("jpsi_fPrompt")->isConstant()){
      latex->DrawLatex(left,top,Form("f^{#psiSR}_{PR}  =  %.3f #pm %.3f",ws->var("jpsi_fPrompt")->getVal(), ws->var("jpsi_fPrompt")->getError()));
      top-=textSize*stepsizeTimes;
    }
    if(!ws->var("jpsi_fBkg")->isConstant()){
      latex->DrawLatex(left,top,Form("f^{#psiSR}_{BG}  =  %.3f #pm %.3f",ws->var("jpsi_fBkg")->getVal(), ws->var("jpsi_fBkg")->getError()));
      top-=textSize*stepsizeTimes;
    }
    if(!ws->var("jpsi_nonPromptTau")->isConstant()){
      latex->DrawLatex(left,top,Form("#tau_{NP}  =  %.3f #pm %.3f mm",ws->var("jpsi_nonPromptTau")->getVal(), ws->var("jpsi_nonPromptTau")->getError()));
      top-=textSize*stepsizeTimes;
    }
    latex->DrawLatex(left,top,Form("f_{B}   =  %.3f ",(1.-ws->var("jpsi_fPrompt")->getVal()-ws->var("jpsi_fBkg")->getVal())/(1.-ws->var("jpsi_fBkg")->getVal())));

  }
  else if(region>0){
    if(!ws->var("jpsi_fBkgSSDR")->isConstant()){
      latex->DrawLatex(left,top,Form("f^{RS}_{BG}  =  %.3f #pm %.3f",ws->var("jpsi_fBkgSSDR")->getVal(), ws->var("jpsi_fBkgSSDR")->getError()));
      top-=textSize*stepsizeTimes;
    }
    //if(!ws->var("jpsi_fBkgDSD")->isConstant()){
    //latex->DrawLatex(left,top,Form("f^{DS}_{BG}  =  %.3f #pm %.3f",ws->var("jpsi_fBkgDSD")->getVal(), ws->var("jpsi_fBkgDSD")->getError()));
    //top-=textSize*stepsizeTimes;
    //}
    if(!ws->var("jpsi_fBkgSSDL")->isConstant()){
      latex->DrawLatex(left,top,Form("f^{LS}_{BG}  =  %.3f #pm %.3f",ws->var("jpsi_fBkgSSDL")->getVal(), ws->var("jpsi_fBkgSSDL")->getError()));
      top-=textSize*stepsizeTimes;
    }
    if(!ws->var("jpsi_bkgTauSSD")->isConstant()){
      latex->DrawLatex(left,top,Form("#tau^{RS}_{BG}  =  %.3f #pm %.3f mm",ws->var("jpsi_bkgTauSSD")->getVal(), ws->var("jpsi_bkgTauSSD")->getError()));
      top-=textSize*stepsizeTimes;
    }
    if(!ws->var("jpsi_bkgTauDSD")->isConstant()){
      latex->DrawLatex(left,top,Form("#tau^{DS}_{BG}  =  %.3f #pm %.3f mm",ws->var("jpsi_bkgTauDSD")->getVal(), ws->var("jpsi_bkgTauDSD")->getError()));
      top-=textSize*stepsizeTimes;
    }
    if(!ws->var("jpsi_bkgTauFD")->isConstant()){
      latex->DrawLatex(left,top,Form("#tau^{LS}_{BG}  =  %.3f #pm %.3f mm",ws->var("jpsi_bkgTauFD")->getVal(), ws->var("jpsi_bkgTauFD")->getError()));
      top-=textSize*stepsizeTimes;
    }

  }

  return;

}
