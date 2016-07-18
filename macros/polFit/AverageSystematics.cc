/*
 * AverageSystematics.cc
 *
 *  Created on: Dec 5, 2011
 *      Author: valentinknuenz
 */

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "rootIncludes.inc"
#include "commonVar.h"
#include "ToyMC.h"
#include "clarg_parsing.h"

#include "/afs/hephy.at/user/t/tmadlener/snippets/vector_stuff.h"

using namespace std;

/** small helper struct to clean up code. */
struct SystematicInput {
  /** constructor. Gets The TFile. Initializes the graph to nullptr. */
  SystematicInput(const std::string& JobID, const std::string& basedir, const std::string& SystID, const std::string& state);
  ~SystematicInput(); /**< Destructor. Closes the TFile. Delets the graph (if present) */
  void loadGraph(const std::string& name);

  std::string jobID; /**< JobID of the input. */
  TGraphAsymmErrors* graph; /**< TGraphAsymmErrors associated with input file. */
  TFile* inFile; /**< Input file. */
};

SystematicInput::SystematicInput(const std::string& JobID, const std::string& basedir, const std::string& SystID, const std::string& state) :
  jobID(JobID),
  graph(NULL)
{
  std::stringstream filename;
  filename << basedir << "/macros/polFit/Systematics/" << SystID << "/" << JobID << "/TGraphResults_" << state << ".root";
  inFile = TFile::Open(filename.str().c_str(), "READ");
}

SystematicInput::~SystematicInput()
{
  if (graph) delete graph;
  // if (inFile && inFile->IsOpen()) inFile->Close(); // DO NOT do this. Cannot be stored in a vector then.
}

void SystematicInput::loadGraph(const std::string& name)
{
  std::cout << jobID << " loading graph " << name.c_str() << std::endl;
  if (graph) delete graph; // do some cleanup before loading the new graph
  graph = dynamic_cast<TGraphAsymmErrors*>(inFile->Get(name.c_str()));
}

/**
 * subtract two TGraphAsymmErrors from each other, returning a new one.
 * subtracts central values of g2 from g1. Errors in x are not touched,
 * errors in y is the absolute value of the difference of the errors of g1 and g2.
 */
TGraphAsymmErrors* subtractTGAE(const TGraphAsymmErrors* g1, const TGraphAsymmErrors* g2)
{
  if (g1->GetN() != g2->GetN()) {
    std::cerr << "cannot subtract TGraphAsymmErrors with different numbers of points!" << std::endl;
    return NULL;
  }

  int nBins = g1->GetN();
  std::vector<double> centValsY; centValsY.reserve(nBins);
  std::vector<double> centValsX; centValsX.reserve(nBins);
  std::vector<double> highErrY; highErrY.reserve(nBins);
  std::vector<double> lowErrY; lowErrY.reserve(nBins);
  std::vector<double> highErrX; highErrX.reserve(nBins);
  std::vector<double> lowErrX; lowErrX.reserve(nBins);

  for (int iBin = 0; iBin < nBins; ++iBin) {
    double c1, c2; // central values
    // double el1, eh1, el2, eh2; // low and high error values
    double x; // dummy value to read in the x-coordinates, which are not used

    if (g1->GetPoint(iBin, x, c1) == iBin && g2->GetPoint(iBin, x, c2) == iBin) {
      std::cout << "point " << iBin << ": x = " << x << ", y1 = " << c1 << ", y2 = " << c2 << std::endl;
      centValsY.push_back(c2 - c1);
      lowErrY.push_back( TMath::Abs(g1->GetErrorYlow(iBin) - TMath::Abs(g2->GetErrorYhigh(iBin))) );
      highErrY.push_back( TMath::Abs(g1->GetErrorYhigh(iBin) - TMath::Abs(g2->GetErrorYhigh(iBin))) );

      // no subtraction in X direction, assume errors are the same and get the errors of g1
      centValsX.push_back(x);
      lowErrX.push_back(g1->GetErrorXlow(iBin));
      highErrX.push_back(g1->GetErrorXhigh(iBin));
    } else {
      std::cerr << "Could not get point " << iBin << " in on of the TGraphsAsymmErrors" << std::endl;
      return NULL;
    }
  }

  return new TGraphAsymmErrors(nBins, centValsX.data(), centValsY.data(),
                               lowErrX.data(), highErrX.data(), lowErrY.data(), highErrY.data());
}

/** multiply all values of the y-coordinate with a factor. */
void multiplyFactor(TGraphAsymmErrors* g, const double factor)
{
  for (int iBin = 0; iBin < g->GetN(); ++iBin) {
    double x, y;
    if (g->GetPoint(iBin, x, y) == iBin) { // only reset the values if retrieval is succesful
      y *= factor;
      double el = g->GetErrorYlow(iBin) * factor;
      double eh = g->GetErrorYhigh(iBin) * factor;

      g->SetPoint(iBin, x, y);
      g->SetPointError(iBin, g->GetErrorXlow(iBin), g->GetErrorXhigh(iBin), el, eh);
    }
  }
}

/** main. (Comment merely fore orientation) */
int main(int argc, char** argv) {
  Char_t *storagedir = "Default"; //Storage Directory
  Char_t *basedir = "Default"; //Storage Directory
  Char_t *SystID = "Default";

  std::vector<SystematicInput> inputs;
  std::vector<std::string> JobIDs;

  int ptBinMin=1;
  int ptBinMax=1;
  int nState=1;
  int nSystematics=1;

  bool subtractGraphs = false;

  for( int i=0;i < argc; ++i ) {
    std::string arg = argv[i];
    fromSplit("subtractGraphs", arg, subtractGraphs);

    std::string tmp = "";
    fromSplit("JobID1", arg, tmp); if (!tmp.empty()) JobIDs.push_back(tmp);
    tmp="";
    fromSplit("JobID2", arg, tmp); if (!tmp.empty()) JobIDs.push_back(tmp);
    tmp="";
    fromSplit("JobID3", arg, tmp); if (!tmp.empty()) JobIDs.push_back(tmp);
    tmp="";
    fromSplit("JobID4", arg, tmp); if (!tmp.empty()) JobIDs.push_back(tmp);
    tmp="";
    fromSplit("JobID5", arg, tmp); if (!tmp.empty()) JobIDs.push_back(tmp);
    tmp="";
    fromSplit("JobID6", arg, tmp); if (!tmp.empty()) JobIDs.push_back(tmp);
    tmp="";
    fromSplit("JobID7", arg, tmp); if (!tmp.empty()) JobIDs.push_back(tmp);
    tmp="";
    fromSplit("JobID8", arg, tmp); if (!tmp.empty()) JobIDs.push_back(tmp);
    tmp="";
    fromSplit("JobID9", arg, tmp); if (!tmp.empty()) JobIDs.push_back(tmp);

    if(std::string(argv[i]).find("SystID") != std::string::npos) {char* SystIDchar = argv[i]; char* SystIDchar2 = strtok (SystIDchar, "="); SystID = SystIDchar2; cout<<"SystID = "<<SystID<<endl;}
    if(std::string(argv[i]).find("storagedir") != std::string::npos) {char* storagedirchar = argv[i]; char* storagedirchar2 = strtok (storagedirchar, "="); storagedir = storagedirchar2; cout<<"storagedir = "<<storagedir<<endl;}
    if(std::string(argv[i]).find("basedir") != std::string::npos) {char* basedirchar = argv[i]; char* basedirchar2 = strtok (basedirchar, "="); basedir = basedirchar2; cout<<"basedir = "<<basedir<<endl;}

    if(std::string(argv[i]).find("ptBinMin") != std::string::npos) {char* ptBinMinchar = argv[i]; char* ptBinMinchar2 = strtok (ptBinMinchar, "p"); ptBinMin = atof(ptBinMinchar2); cout<<"ptBinMin = "<<ptBinMin<<endl;}
    if(std::string(argv[i]).find("ptBinMax") != std::string::npos) {char* ptBinMaxchar = argv[i]; char* ptBinMaxchar2 = strtok (ptBinMaxchar, "p"); ptBinMax = atof(ptBinMaxchar2); cout<<"ptBinMax = "<<ptBinMax<<endl;}
    if(std::string(argv[i]).find("nState") != std::string::npos) {char* nStatechar = argv[i]; char* nStatechar2 = strtok (nStatechar, "p"); nState = atof(nStatechar2); cout<<"nState = "<<nState<<endl;}
    if(std::string(argv[i]).find("nSystematics") != std::string::npos) {char* nSystematicschar = argv[i]; char* nSystematicschar2 = strtok (nSystematicschar, "p"); nSystematics = atof(nSystematicschar2); cout<<"nSystematics = "<<nSystematics<<endl;}

  }

  std::stringstream stateStr;
  if (nState < 6) stateStr << "Psi" << nState - 3 << "S";
  else stateStr << "chic" << nState - 5;

  char tmpfilename[200];
  sprintf(tmpfilename,"%s/macros/polFit/Systematics/%s/AverageSyst/TGraphResults_%s.root",basedir,SystID,stateStr.str().c_str());
  gSystem->Unlink(tmpfilename);

  for (size_t iID = 0; iID < JobIDs.size(); ++iID) {
    inputs.push_back(SystematicInput(JobIDs[iID], basedir, SystID, stateStr.str()));
  }

  char GraphName[200];
  char filename[200];
  int nRapBins = 1;
  if(nState==5) nRapBins =  3;
  cout<<"nRapBins: "<<nRapBins<<endl;

  for(int iLam = 1; iLam<19; iLam++){

    for(int rapBin = 1; rapBin <= nRapBins; rapBin++){

      sprintf(filename,"%s/macros/polFit/Systematics/%s/AverageSyst/TGraphResults_%s.root",basedir,SystID,stateStr.str().c_str());
      TFile *outfile = new TFile(filename,"UPDATE");


      if(iLam==1)  sprintf(GraphName,"lth_CS_rap%d",rapBin);
      if(iLam==2)  sprintf(GraphName,"lph_CS_rap%d",rapBin);
      if(iLam==3)  sprintf(GraphName,"ltp_CS_rap%d",rapBin);
      if(iLam==4)  sprintf(GraphName,"lthstar_CS_rap%d",rapBin);
      if(iLam==5)  sprintf(GraphName,"lphstar_CS_rap%d",rapBin);
      if(iLam==6)  sprintf(GraphName,"ltilde_CS_rap%d",rapBin);

      if(iLam==7)  sprintf(GraphName,"lth_HX_rap%d",rapBin);
      if(iLam==8)  sprintf(GraphName,"lph_HX_rap%d",rapBin);
      if(iLam==9)  sprintf(GraphName,"ltp_HX_rap%d",rapBin);
      if(iLam==10) sprintf(GraphName,"lthstar_HX_rap%d",rapBin);
      if(iLam==11) sprintf(GraphName,"lphstar_HX_rap%d",rapBin);
      if(iLam==12) sprintf(GraphName,"ltilde_HX_rap%d",rapBin);

      if(iLam==13) sprintf(GraphName,"lth_PX_rap%d",rapBin);
      if(iLam==14) sprintf(GraphName,"lph_PX_rap%d",rapBin);
      if(iLam==15) sprintf(GraphName,"ltp_PX_rap%d",rapBin);
      if(iLam==16) sprintf(GraphName,"lthstar_PX_rap%d",rapBin);
      if(iLam==17) sprintf(GraphName,"lphstar_PX_rap%d",rapBin);
      if(iLam==18) sprintf(GraphName,"ltilde_PX_rap%d",rapBin);
      cout<<"GraphName: "<<GraphName<<endl;

      // new way of input handling
      for (size_t iInput = 0; iInput < inputs.size(); ++iInput) {
        inputs[iInput].loadGraph(GraphName);
      }

      int nBinspT=ptBinMax-ptBinMin+1;
      double ptCentre_[nBinspT];
      double ptCentreErr_low[nBinspT];
      double ptCentreErr_high[nBinspT];
      double lmean[nBinspT];
      double errlmean[nBinspT];

      double fit_lmean[nSystematics];
      double fit_errlmean[nSystematics];
      double fit_X[nSystematics];

      double ptCentre__[nSystematics][nBinspT];
      double ptCentreErr_low_[nSystematics][nBinspT];
      double ptCentreErr_high_[nSystematics][nBinspT];

      double lmean_[nSystematics][nBinspT];

      TGraphAsymmErrors* graph_;

      if(!subtractGraphs) {
        int pt=0;
        for(int ptBin = ptBinMin; ptBin < ptBinMax+1; ptBin++) {

          double lmean_Buffer=0;
          double ptCentre_Buffer=0;
          double ptCentreErr_low_Buffer=0;
          double ptCentreErr_high_Buffer=0;

          for(int iSyst = 0; iSyst < nSystematics && iSyst < inputs.size(); iSyst++) {
            graph_ = inputs[iSyst].graph;
            graph_->GetPoint(pt,ptCentre__[iSyst][pt],lmean_[iSyst][pt]);
            ptCentreErr_high_[iSyst][pt]=graph_->GetErrorXhigh(pt);
            ptCentreErr_low_[iSyst][pt]=graph_->GetErrorXlow(pt);

            fit_errlmean[iSyst]=graph_->GetErrorYhigh(pt);
            fit_lmean[iSyst]=lmean_[iSyst][pt];
            fit_X[iSyst]=0.;
          }

          double pTcentreReal=0;
          double pTcentreReallow=0;
          double pTcentreRealhigh=0;

          double lmeanHighest = 0.;
          for(int iSyst=0;iSyst<nSystematics;iSyst++){
            //lmean_Buffer=lmean_Buffer+TMath::Abs(lmean_[iSyst][pt]); //ifMean
            //lmean_Buffer=lmean_Buffer+lmean_[iSyst][pt]; //ifChange
            //lmean_[iSyst][pt] = lmean_[iSyst][pt] / 2. ;  //fracBg0_TO_default_half
            lmean_Buffer=lmean_Buffer+lmean_[iSyst][pt]*lmean_[iSyst][pt]; //ifSquare
            //if(ptBin<10 && iSyst==0) lmean_Buffer = lmean_Buffer + lmean_[iSyst][pt]*lmean_[iSyst][pt];
            //if(ptBin>=10 && iSyst==1)  lmean_Buffer = lmean_Buffer + lmean_[iSyst][pt]*lmean_[iSyst][pt];

            cout<<"lmean["<<iSyst<<"]["<<pt<<"]: "<<lmean_[iSyst][pt]<<endl;
            //if(lmean_[iSyst][pt]*lmean_[iSyst][pt]>lmeanHighest){
            //	lmeanHighest = lmean_[iSyst][pt]*lmean_[iSyst][pt] ;
            //	lmean_Buffer = lmeanHighest ;
            //	cout<<"lmeanHighest: "<<sqrt(lmeanHighest)<<endl;
            //}

            ptCentre_Buffer=ptCentre_Buffer+ptCentre__[iSyst][pt];
            ptCentreErr_low_Buffer=ptCentreErr_low_Buffer+ptCentreErr_low_[iSyst][pt];
            ptCentreErr_high_Buffer=ptCentreErr_high_Buffer+ptCentreErr_high_[iSyst][pt];

            if(iSyst==0) {//delete loop
              pTcentreReal=ptCentre__[iSyst][pt];
              pTcentreReallow=ptCentreErr_low_[iSyst][pt];
              pTcentreRealhigh=ptCentreErr_high_[iSyst][pt];
            }

          }

          ptCentre_[pt]=ptCentre_Buffer/nSystematics;
          ptCentreErr_low[pt]=ptCentreErr_low_Buffer/nSystematics;
          ptCentreErr_high[pt]=ptCentreErr_high_Buffer/nSystematics;

          ptCentre_[pt]=pTcentreReal;//delete
          ptCentreErr_low[pt]=pTcentreReallow;//delete
          ptCentreErr_high[pt]=pTcentreRealhigh;//delete

          //lmean[pt]=lmean_Buffer/nSystematics; //ifMean
          //if(pt>3) {lmean[pt]=0; std::cout << "point 4 set to 0" << std::endl;} // ifrho
          lmean[pt]=TMath::Sqrt(lmean_Buffer); //ifSquare
          //lmean[pt]=lmean_Buffer;

          std::cout << pt << ": pT = " << ptCentre_[pt] << ", lambda = " << lmean[pt] << std::endl;

          // IF FIT THE NSYSTEMATIC VALUES INSTEAD OF USING THE MEAN:::

          //if(pt==9) {
          //      fit_X[2]=-999.;
          //}

          //TGraphAsymmErrors *fitGraph = new TGraphAsymmErrors(nSystematics,fit_X,fit_lmean,0,0,fit_errlmean,fit_errlmean);
          //TF1* fConst = new TF1("fConst","pol0",-1,1);
          //fitGraph->Fit("fConst","EFNR");
          //
          //lmean[pt]=fConst->GetParameter(0);
          //errlmean[pt]=fConst->GetParError(0);

          // END Fit

          pt++;
        } // pT-Loop
      } else { // i.e. if subtractGraph
        TGraphAsymmErrors* diffSyst = subtractTGAE(inputs[0].graph, inputs[1].graph);
        multiplyFactor(diffSyst, TMath::Sqrt(12.));
        for (int i = 0; i < nBinspT; ++i) {
          diffSyst->GetPoint(i, ptCentre_[i], lmean[i]); // WARNING: assuming nBinspT and diffSyst->GetN() are the same
          ptCentreErr_low[i] = diffSyst->GetErrorXlow(i);
          ptCentreErr_high[i] = diffSyst->GetErrorXhigh(i);
        }
      }
      //////////////// Change TGraphs
      /*
        double ptCentre_Change[nBinspT-2];
        double ptCentreErr_low_Change[nBinspT-2];
        double ptCentreErr_highChange[nBinspT-2];
        double lmeanChange[nBinspT-2];

        pt=0;
        for(int ptBin = ptBinMin; ptBin < ptBinMax+1-2; ptBin++) {
        ptCentre_Change[pt]=ToyMC::ptCentre[rapBin-1][pt];
        ptCentreErr_low_Change[pt]=TMath::Abs(ToyMC::ptCentre[rapBin-1][pt]-onia::pTRange[rapBin][pt]);
        ptCentreErr_highChange[pt]=TMath::Abs(ToyMC::ptCentre[rapBin-1][pt]-onia::pTRange[rapBin][pt+1]);

        if(pt<5)lmeanChange[pt]=lmean_[0][pt];
        if(pt==5)lmeanChange[pt]=(lmean_[0][5]+lmean_[0][6])/2.;
        if(pt==6)lmeanChange[pt]=(lmean_[0][7]+lmean_[0][8])/2.;
        if(pt==7)lmeanChange[pt]=lmean_[0][9];
        if(pt==8)lmeanChange[pt]=lmean_[0][10];
        if(pt==9)lmeanChange[pt]=lmean_[0][11];

        pt++;
        }


        TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT-2,ptCentre_Change,lmeanChange,ptCentreErr_low_Change,ptCentreErr_highChange,0,0);
        graphSyst->SetMarkerColor(ToyMC::MarkerColor[rapBin]);
        graphSyst->SetLineColor(ToyMC::MarkerColor[rapBin]);
        //              graphSyst->SetMarkerStyle(ToyMC::MarkerStyle[rapBin]);
        graphSyst->SetMarkerSize(2);
        graphSyst->SetName(GraphName);
      */
      //////////////// END Change TGraphs

      TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean,ptCentreErr_low,ptCentreErr_high,0,0);//Original
      //TGraphAsymmErrors *graphSyst = new TGraphAsymmErrors(nBinspT,ptCentre_,lmean,ptCentreErr_low,ptCentreErr_high,errlmean,errlmean);//If Fit the nSyst values with constant
      graphSyst->SetMarkerColor(ToyMC::MarkerColor[rapBin]);
      graphSyst->SetLineColor(ToyMC::MarkerColor[rapBin]);
      //graphSyst->SetMarkerStyle(ToyMC::MarkerStyle[rapBin]);
      graphSyst->SetMarkerSize(2);
      graphSyst->SetName(GraphName);

      outfile->cd();
      graphSyst->Draw("P");
      graphSyst->Write();

      outfile->Write();
      outfile->Close();
      delete outfile;
      outfile = NULL;


    }


  }

  return 0;
}
