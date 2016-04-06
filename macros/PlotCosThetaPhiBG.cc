#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"

#include "TH2F.h"

enum{L,R};
const char *bgLabel[2] = {"L", "R"};
TH2F *hCosThetaPhi[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::kNbFrames][2][2];
TH2F *hCosThetaPhiHighct[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::kNbFrames][2][2];
TH2F *hCosThetaPhiNPBG[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::kNbFrames][2];
TH2F *hTCosThetaPhi[onia::kNbRapForPTBins][onia::kNbPTMaxBins][onia::kNbFrames][2];
TH1F *events_PRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];
TH1F *events_NPSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];
TH1F *hMeanPt[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];
TH1F *hMeanY[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];
TH1F *hfBGinPRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];
TH1F *hfNPinPRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];

int evtPinPRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2],
  evtNPinPRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2],
  evtBGinPRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2],
  evtPinNPSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2],
  evtNPinNPSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2],
  evtBGinNPSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];

double fracBGinPRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2],
  fracNPinPRSR[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];
double meanPt[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2],
  meanRap[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];

void LoadHistos(int iRapBin, int iPTBin, int nState);
void PlotHistos(int iRapBin, int iPTBin, int nState, int iFrame, int iWindow);

//========================================================
// code to read input arguments
template<typename T>
void fromSplit(const std::string& key, const std::string &arg, T& out)
{
  const char delim = '=';
  // Skip if key or delimiter not there
  if ((arg.find(key) == std::string::npos) ||
      (arg.find(delim) == std::string::npos))
    return;

  std::string skey, sval;
  std::stringstream sstr(arg);
  std::getline(sstr, skey, delim); // Dummy read to skip key
  std::getline(sstr, sval, delim); // Get value
  T tout;
  if (!(std::istringstream(sval) >> std::boolalpha >> tout))
    return;
  out = tout;
  std::cout << std::boolalpha << skey << ": "  << out << std::endl;
}
// Special version for string without the conversion
template<>
void fromSplit(const std::string& key, const std::string &arg, std::string &out)
{
  const char delim = '=';
  // Skip if key or delimiter not there
  if ((arg.find(key) == std::string::npos) ||
      (arg.find(delim) == std::string::npos))
    return;
  std::string skey, sval;
  std::stringstream sstr(arg);
  std::getline(sstr, skey, delim); // Dummy read to skip key
  std::getline(sstr, sval, delim); // Get value
  out = sval;
  std::cout << skey << ": "  << out << std::endl;
}

//=================== main =====================
int main(int argc, char* argv[]){

  // set default values
  int nState = 999;

  // Loop over argument list
  for (int i=1; i < argc; i++){
    std::string arg = argv[i];
    fromSplit("nState", arg, nState);
  }

  gStyle->SetPalette(1);

  for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
    for(int iPT = 0; iPT < onia::kNbPTBins[iRap+1]; iPT++){
      LoadHistos(iRap, iPT, nState);
      for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
        if(iFrame<3){
          PlotHistos(iRap, iPT, nState, iFrame, L);
          PlotHistos(iRap, iPT, nState, iFrame, R);
        }
      }
    }
  }


  FILE *NumFile, *NumFile1, *NumFilePt, *NumFilePt1;

  bool SaveTables=true;
  if(SaveTables){
    for(int iRap = 0; iRap < onia::kNbRapForPTBins; iRap++){
      for(int iPT = 0; iPT < onia::kNbPTBins[iRap+1]; iPT++){
        if(nState == 6){
          std::cout << "--- chic1 ---"
                    <<"evtPinPRSR  ["<<iRap<<"]["<<iPT<<"]: "<<evtPinPRSR[iRap][iPT][0] << "\n"
                    <<"evtNPinPRSR ["<<iRap<<"]["<<iPT<<"]: "<<evtNPinPRSR[iRap][iPT][0] << "\n"
                    <<"evtBGinPRSR ["<<iRap<<"]["<<iPT<<"]: "<<evtBGinPRSR[iRap][iPT][0] << "\n"
                    <<"fracBGinPRSR["<<iRap<<"]["<<iPT<<"]: "<<fracBGinPRSR[iRap][iPT][0] << "\n"
                    <<"evtPinNPSR  ["<<iRap<<"]["<<iPT<<"]: "<<evtPinNPSR[iRap][iPT][0] << "\n"
                    <<"evtNPinNPSR ["<<iRap<<"]["<<iPT<<"]: "<<evtNPinNPSR[iRap][iPT][0] << "\n"
                    <<"evtBGinNPSR ["<<iRap<<"]["<<iPT<<"]: "<<evtBGinNPSR[iRap][iPT][0] << std::endl;
          std::cout << "--- chic2 ---"
                    <<"evtPinPRSR  ["<<iRap<<"]["<<iPT<<"]: "<<evtPinPRSR[iRap][iPT][1] << "\n"
                    <<"evtNPinPRSR ["<<iRap<<"]["<<iPT<<"]: "<<evtNPinPRSR[iRap][iPT][1] << "\n"
                    <<"evtBGinPRSR ["<<iRap<<"]["<<iPT<<"]: "<<evtBGinPRSR[iRap][iPT][1] << "\n"
                    <<"fracBGinPRSR["<<iRap<<"]["<<iPT<<"]: "<<fracBGinPRSR[iRap][iPT][1] << "\n"
                    <<"evtPinNPSR  ["<<iRap<<"]["<<iPT<<"]: "<<evtPinNPSR[iRap][iPT][1] << "\n"
                    <<"evtNPinNPSR ["<<iRap<<"]["<<iPT<<"]: "<<evtNPinNPSR[iRap][iPT][1] << "\n"
                    <<"evtBGinNPSR ["<<iRap<<"]["<<iPT<<"]: "<<evtBGinNPSR[iRap][iPT][1] << std::endl;
        } else {
          cout<<"evtPinPRSR  ["<<iRap<<"]["<<iPT<<"]: "<<evtPinPRSR[iRap][iPT][0]<<endl;
          cout<<"evtNPinPRSR ["<<iRap<<"]["<<iPT<<"]: "<<evtNPinPRSR[iRap][iPT][0]<<endl;
          cout<<"evtBGinPRSR ["<<iRap<<"]["<<iPT<<"]: "<<evtBGinPRSR[iRap][iPT][0]<<endl;
          cout<<"fracBGinPRSR["<<iRap<<"]["<<iPT<<"]: "<<fracBGinPRSR[iRap][iPT][0]<<endl;

          cout<<"evtPinNPSR  ["<<iRap<<"]["<<iPT<<"]: "<<evtPinNPSR[iRap][iPT][0]<<endl;
          cout<<"evtNPinNPSR ["<<iRap<<"]["<<iPT<<"]: "<<evtNPinNPSR[iRap][iPT][0]<<endl;
          cout<<"evtBGinNPSR ["<<iRap<<"]["<<iPT<<"]: "<<evtBGinNPSR[iRap][iPT][0]<<endl;
        }

      }
    }

    std::cout << "write numbers to tex file" << std::endl;
    /// write estimated number signal events and Bg fractions to file
    std::stringstream numFileName, numFileName1;
    if(nState == 6){
      numFileName << "NumEventsFracBG_chic1.tex";
      numFileName1 << "NumEventsFracBG_chic2.tex";
    } else
      numFileName << "NumEventsFracBG_Psi" << nState-3 << ".tex";

    char framerap[200];

    NumFile = fopen(numFileName.str().c_str(),"w");
    fprintf(NumFile, "\n");
    fprintf(NumFile,"\\documentclass{article}\n\\usepackage[applemac]{inputenc}\n\\usepackage{amsmath}\n\\usepackage{textcomp}\n\\pagestyle{plain}\n\\usepackage{graphicx}\n\\usepackage{multicol}\n\\usepackage{geometry}\n\\usepackage{subfigure}\n\\usepackage{booktabs}\n\\usepackage{setspace}\n\n\n\n\\begin{document}\n");

    fprintf(NumFile, "\n\n\n\n");

    fprintf(NumFile, "\n\n\n\n");

    fprintf(NumFile, "\\begin{table}[!H]\n\\centering\n \\caption{Number of prompt signal events, non-prompt fraction and background fraction in the PRSR, as a function of $p_{T}$ and $y$. The fractions are in \\% . } \n");

    // Psi states
    if(nState != 6){
      fprintf(NumFile, "\\begin{tabular}{|c|ccc|ccc|ccc|}\n\\hline\n");
      fprintf(NumFile, "$p_{T}$ [GeV] & $N_{PR}^{PRSR}$ & $f_{NP}^{PRSR} $ & $f_{BG}^{PRSR} $ & $N_{PR}^{PRSR}$ & $f_{NP}^{PRSR} $ & $f_{BG}^{PRSR} $ & $N_{PR}^{PRSR}$ & $f_{NP}^{PRSR} $ & $f_{BG}^{PRSR} $  \\\\\n");

      if(nState==4){
        sprintf(framerap,"\\hline \\multicolumn{1}{|c|}{} & \\multicolumn{3}{|c|}{$%1.1f < |y| < %1.1f$ } & \\multicolumn{3}{|c|}{$%1.1f < |y| < %1.1f$ } & \\multicolumn{3}{|c|}{$%1.1f < |y| < %1.1f$ } \\\\\n \\hline \n",0.0, 0.6, 0.6, 1.2, 1.2, 1.5);
        fprintf(NumFile,framerap);
        fprintf(NumFile, "\\multicolumn{10}{|c|}{$\\Psi(1S)$} \\\\\n \\hline \n \\rule{0pt}{4mm} \n");
        int pt=0;
        for(int ptBin = 1; ptBin < onia::kNbPTBins[1]+1; ptBin++) {
          fprintf(NumFile, "%1.0f--%1.0f   &  $%d  $  & $%1.1f $  &  $%1.1f $ &  $%d $ &  $%1.1f $ &  $%1.1f $ & -- & -- & -- \\\\\n",
                  onia::pTRange[1][ptBin-1], onia::pTRange[1][ptBin],
                  evtPinPRSR[0][ptBin-1][0], 100.*fracNPinPRSR[0][ptBin-1][0], 100*fracBGinPRSR[0][ptBin-1][0],
                  evtPinPRSR[1][ptBin-1][0], 100.*fracNPinPRSR[1][ptBin-1][0], 100*fracBGinPRSR[1][ptBin-1][0]);
          pt++;
        }
      }

      else if(nState==5){
        sprintf(framerap,"\\hline \\multicolumn{1}{|c|}{} & \\multicolumn{3}{|c|}{$%1.1f < |y| < %1.1f$ } & \\multicolumn{3}{|c|}{$%1.1f < |y| < %1.1f$} & \\multicolumn{3}{|c|}{$%1.1f < |y| < %1.1f$} \\\\\n \\hline \n",onia::rapForPTRange[0],onia::rapForPTRange[1], onia::rapForPTRange[1], onia::rapForPTRange[2], onia::rapForPTRange[2], onia::rapForPTRange[3]);
        fprintf(NumFile,framerap);
        fprintf(NumFile, "\\multicolumn{10}{|c|}{$\\Psi(2S)$} \\\\\n \\hline \n \\rule{0pt}{4mm} \n");

        int pt=0;
        for(int ptBin = 1; ptBin < onia::kNbPTBins[1]+1; ptBin++) {

          fprintf(NumFile, "%1.0f--%1.0f   &  $%d  $  & $%1.1f $  &  $%1.1f $ &  $%d $ &  $%1.1f $ &  $%1.1f $ &  $%d $ &  $%1.1f $ &  $%1.1f $ \\\\\n",
                  onia::pTRange[1][ptBin-1], onia::pTRange[1][ptBin],
                  evtPinPRSR[0][ptBin-1][0], 100.*fracNPinPRSR[0][ptBin-1][0], 100*fracBGinPRSR[0][ptBin-1][0],
                  evtPinPRSR[1][ptBin-1][0], 100.*fracNPinPRSR[1][ptBin-1][0], 100*fracBGinPRSR[1][ptBin-1][0],
                  evtPinPRSR[2][ptBin-1][0], 100.*fracNPinPRSR[2][ptBin-1][0], 100*fracBGinPRSR[2][ptBin-1][0] );
          pt++;
        }
      }
    }
            // chi1 state
    else if(nState==6){
      fprintf(NumFile, "\\begin{tabular}{|c|c|c|c|}\n\\hline\n");
      fprintf(NumFile, "\\multicolumn{4}{|c|}{$\\chi_{c1}$} \\\\\n \\hline \n");
      sprintf(framerap,"\\multicolumn{4}{|c|}{$%1.1f < |y| < %1.1f$ } \\\\\n \\hline \n",0.0,1.2);
      fprintf(NumFile, framerap);
      fprintf(NumFile, "$p_{T}$ [GeV] & $N_{PR}^{PRSR}$ & $f_{NP}^{PRSR} $ & $f_{BG}^{PRSR} $  \\\\\n \\hline \n \\rule{0pt}{4mm} \n");
      int pt=0;
      for(int ptBin = 1; ptBin < onia::kNbPTBins[1]+1; ptBin++) {
        fprintf(NumFile, "%1.0f--%1.0f   &  $%d  $  & $%1.1f $  &  $%1.1f $ \\\\\n",
                onia::pTRange[1][ptBin-1], onia::pTRange[1][ptBin],
                evtPinPRSR[0][ptBin-1][0], 100.*fracNPinPRSR[0][ptBin-1][0], 100*fracBGinPRSR[0][ptBin-1][0]);
        pt++;
      }
    }
    fprintf(NumFile, "\\hline\n");
    fprintf(NumFile, "\\end{tabular}\n");
    fprintf(NumFile, "\\label{tab:NumEventsFracBG}\n");
    fprintf(NumFile, "\\end{table}\n");
    fprintf(NumFile, "\n");

    fprintf(NumFile, "\\end{document}");

    fclose(NumFile);
    std::cout << "first file written" << std::endl;

    // chi2 state
    if(nState == 6){
      std::cout << "open second file" << std::endl;
      NumFile1 = fopen(numFileName1.str().c_str(),"w");
      fprintf(NumFile1, "\n");
      fprintf(NumFile1,"\\documentclass{article}\n\\usepackage[applemac]{inputenc}\n\\usepackage{amsmath}\n\\usepackage{textcomp}\n\\pagestyle{plain}\n\\usepackage{graphicx}\n\\usepackage{multicol}\n\\usepackage{geometry}\n\\usepackage{subfigure}\n\\usepackage{booktabs}\n\\usepackage{setspace}\n\n\n\n\\begin{document}\n");
      fprintf(NumFile1, "\n\n\n\n");
      fprintf(NumFile1, "\n\n\n\n");
      fprintf(NumFile1, "\\begin{table}[!H]\n\\centering\n \\caption{Number of prompt signal events, non-prompt fraction and background fraction in the PRSR, as a function of $p_{T}$ and $y$. The fractions are in \\%. }\n");
      fprintf(NumFile1, "\\begin{tabular}{|c|c|c|c|}\n\\hline\n");
      fprintf(NumFile1, "\\multicolumn{4}{|c|}{$\\chi_{c2}$} \\\\\n \\hline \n");
      sprintf(framerap,"\\multicolumn{4}{|c|}{$%1.1f < |y| < %1.1f$ }\\\\\n \\hline \n",0.0,1.2);
      fprintf(NumFile1,framerap);
      fprintf(NumFile1, "$p_{T}$ [GeV] & $N_{PR}^{PRSR}$ & $f_{NP}^{PRSR} $ & $f_{BG}^{PRSR} $  \\\\\n \\hline \n \\rule{0pt}{4mm} \n");
      int pt=0;
      for(int ptBin = 1; ptBin < onia::kNbPTBins[1]+1; ptBin++) {
        fprintf(NumFile1, "%1.0f--%1.0f   &  $%d  $  & $%1.1f $  &  $%1.1f $ \\\\\n",
                onia::pTRange[1][ptBin-1], onia::pTRange[1][ptBin],
                evtPinPRSR[0][ptBin-1][1], 100.*fracNPinPRSR[0][ptBin-1][1], 100*fracBGinPRSR[0][ptBin-1][1]);
        pt++;
      }
      fprintf(NumFile1, "\\hline\n");
      fprintf(NumFile1, "\\end{tabular}\n");
      fprintf(NumFile1, "\\label{tab:NumEventsFracBG}\n");
      fprintf(NumFile1, "\\end{table}\n");
      fprintf(NumFile1, "\n");
      fprintf(NumFile1, "\\end{document}");
      fclose(NumFile1);
      std::cout << "second file written" << std::endl;
    }

    //// write estimated mean pT and y to file
    std::stringstream numFileNamePt, numFileNamePt1;
    if(nState == 6){
      numFileNamePt << "meanPt_chic1.tex";
      numFileNamePt1 << "meanPt_chic2.tex";
    } else
      numFileNamePt << "meanPt_Psi" << nState-3 << ".tex";

    NumFilePt = fopen(numFileNamePt.str().c_str(),"w");
    fprintf(NumFilePt, "\n");
    fprintf(NumFilePt,"\\documentclass{article}\n\\usepackage[applemac]{inputenc}\n\\usepackage{amsmath}\n\\usepackage{textcomp}\n\\pagestyle{plain}\n\\usepackage{graphicx}\n\\usepackage{multicol}\n\\usepackage{geometry}\n\\usepackage{subfigure}\n\\usepackage{booktabs}\n\\usepackage{setspace}\n\n\n\n\\begin{document}\n");
    fprintf(NumFilePt, "\n\n\n\n");
    fprintf(NumFilePt, "\n\n\n\n");
    fprintf(NumFilePt, "\\begin{table}[!H]\n\\centering\n \\caption{Estimated mean $p_{T}$ and $|y|$, for each $\\psi(nS)$ kinematic bin} \n");



    // psi states
    if(nState != 6){
      fprintf(NumFilePt, "\\begin{tabular}{|c|cc|cc|cc|}\n\\hline\n");
      fprintf(NumFilePt, "$p_{T}$ [GeV] & $\\hat{p_{T}}$ [GeV] & $\\hat{|y|}$ & $\\hat{p_{T}}$ [GeV] & $\\hat{|y|}$ & $\\hat{p_{T}}$ [GeV] & $\\hat{|y|}$ \\\\\n");

      if(nState==4){
        sprintf(framerap,"\\hline \\multicolumn{1}{|c|}{} & \\multicolumn{2}{|c|}{$%1.1f < |y| < %1.1f$ } & \\multicolumn{2}{|c|}{$%1.1f < |y| < %1.1f$ } & \\multicolumn{2}{|c|}{$%1.1f < |y| < %1.1f$ } \\\\\n \\hline \n",0.0, 0.6, 0.6, 1.2, 1.2, 1.5);
        fprintf(NumFilePt,framerap);
        fprintf(NumFilePt, "\\multicolumn{7}{|c|}{$\\Psi(1S)$} \\\\\n \\hline \n \\rule{0pt}{4mm} \n");
        int pt=0;
        for(int ptBin = 1; ptBin < onia::kNbPTBins[1]+1; ptBin++) {

          fprintf(NumFilePt, "%1.0f--%1.0f & $%1.3f$ & $%1.3f$ & $%1.3f$ & $%1.3f$ & -- & -- \\\\\n",
                  onia::pTRange[1][ptBin-1], onia::pTRange[1][ptBin],
                  meanPt[0][ptBin-1][0], meanRap[0][ptBin-1][0],
                  meanPt[1][ptBin-1][0], meanRap[1][ptBin-1][0]);
          pt++;
        }
      }

      if(nState==5){
        sprintf(framerap,"\\hline \\multicolumn{1}{|c|}{} & \\multicolumn{2}{|c|}{$%1.1f < |y| < %1.1f$ } & \\multicolumn{2}{|c|}{$%1.1f < |y| < %1.1f$} & \\multicolumn{2}{|c|}{$%1.1f < |y| < %1.1f$} \\\\\n \\hline \n",onia::rapForPTRange[0],onia::rapForPTRange[1], onia::rapForPTRange[1], onia::rapForPTRange[2], onia::rapForPTRange[2], onia::rapForPTRange[3]);
        fprintf(NumFilePt,framerap);
        fprintf(NumFilePt, "\\multicolumn{7}{|c|}{$\\Psi(2S)$} \\\\\n \\hline \n \\rule{0pt}{4mm} \n");

        int pt=0;
        for(int ptBin = 1; ptBin < onia::kNbPTBins[1]+1; ptBin++) {

          fprintf(NumFilePt, "%1.0f--%1.0f & $%1.3f$ & $%1.3f$ & $%1.3f$ & $%1.3f$ & $%1.3f$ & $%1.3f$ \\\\\n",
                  onia::pTRange[1][ptBin-1], onia::pTRange[1][ptBin],
                  meanPt[0][ptBin-1][0], meanRap[0][ptBin-1][0],
                  meanPt[1][ptBin-1][0], meanRap[1][ptBin-1][0],
                  meanPt[2][ptBin-1][0], meanRap[2][ptBin-1][0]);
          pt++;
        }
      }
    }
    //chi1 state
    else if(nState == 6){
      fprintf(NumFilePt, "\\begin{tabular}{|c|c|c|}\n\\hline\n");
      fprintf(NumFilePt, "\\multicolumn{3}{|c|}{$\\chi_{c1}$} \\\\\n \\hline \n");
      sprintf(framerap,"\\multicolumn{3}{|c|}{$%1.1f < |y| < %1.1f$ }\\\\\n \\hline \n",0.0,1.2);
      fprintf(NumFilePt,framerap);
      fprintf(NumFilePt, "$p_{T}$ [GeV] & $\\hat{p_{T}}$ [GeV] & $\\hat{|y|}$ \\\\\n \\hline \n \\rule{0pt}{4mm} \n");
      int pt=0;
      for(int ptBin = 1; ptBin < onia::kNbPTBins[1]+1; ptBin++) {
        fprintf(NumFilePt, "%1.0f--%1.0f & $%1.3f$ & $%1.3f$ \\\\\n",
                onia::pTRange[1][ptBin-1], onia::pTRange[1][ptBin],
                meanPt[0][ptBin-1][0], meanRap[0][ptBin-1][0]);
        pt++;
      }
    }
    fprintf(NumFilePt, "\\hline\n");
    fprintf(NumFilePt, "\\end{tabular}\n");
    fprintf(NumFilePt, "\\label{tab:NumEventsFracBG}\n");
    fprintf(NumFilePt, "\\end{table}\n");
    fprintf(NumFilePt, "\n");
    fprintf(NumFilePt, "\\end{document}");
    fclose(NumFilePt);
    std::cout << "mean pt and y file written" << std::endl;

  // chi2 state
  if(nState == 6){
    NumFilePt1 = fopen(numFileNamePt1.str().c_str(),"w");
    fprintf(NumFilePt1, "\n");
    fprintf(NumFilePt1,"\\documentclass{article}\n\\usepackage[applemac]{inputenc}\n\\usepackage{amsmath}\n\\usepackage{textcomp}\n\\pagestyle{plain}\n\\usepackage{graphicx}\n\\usepackage{multicol}\n\\usepackage{geometry}\n\\usepackage{subfigure}\n\\usepackage{booktabs}\n\\usepackage{setspace}\n\n\n\n\\begin{document}\n");
    fprintf(NumFilePt1, "\n\n\n\n");
    fprintf(NumFilePt1, "\n\n\n\n");
    fprintf(NumFilePt1, "\\begin{table}[!H]\n\\centering\n \\caption{Estimated mean $p_{T}$ and $|y|$, for each $\\psi(nS)$ kinematic bin} \n");
    fprintf(NumFilePt1, "\\begin{tabular}{|c|c|c|}\n\\hline\n");
    fprintf(NumFilePt1, "\\multicolumn{3}{|c|}{$\\chi_{c2}$} \\\\\n \\hline \n");
    sprintf(framerap,"\\multicolumn{3}{|c|}{$%1.1f < |y| < %1.1f$ }\\\\\n \\hline \n",0.0,1.2);
    fprintf(NumFilePt1,framerap);
    fprintf(NumFilePt1, "$p_{T}$ [GeV] & $\\hat{p_{T}}$ [GeV] & $\\hat{|y|}$ \\\\\n \\hline \n \\rule{0pt}{4mm} \n");
    int pt=0;
    for(int ptBin = 1; ptBin < onia::kNbPTBins[1]+1; ptBin++) {
      fprintf(NumFilePt1, "%1.0f--%1.0f & $%1.3f$ & $%1.3f$ \\\\\n",
              onia::pTRange[1][ptBin-1], onia::pTRange[1][ptBin],
              meanPt[0][ptBin-1][1], meanRap[0][ptBin-1][1]);
      pt++;
    }
    fprintf(NumFilePt1, "\\hline\n");
    fprintf(NumFilePt1, "\\end{tabular}\n");
    fprintf(NumFilePt1, "\\label{tab:NumEventsFracBG}\n");
    fprintf(NumFilePt1, "\\end{table}\n");
    fprintf(NumFilePt1, "\n");
    fprintf(NumFilePt1, "\\end{document}");
    fclose(NumFilePt1);
    std::cout << "second mean pt and y file written" << std::endl;
  }
  std::cout << "tex files written" << std::endl;
}
std::cout << "end of program" << std::endl;
return 0;
}

//===========================
void PlotHistos(int iRapBin, int iPTBin, int nState, int iFrame, int iWindow){

  TGaxis::SetMaxDigits(3);

  double lvalue = 0.28, tvalue = 0.92;
  double left=lvalue, top=tvalue, textSize=0.035;
  TLatex *latex=new TLatex();
  latex->SetTextFont(42);
  latex->SetNDC(kTRUE);
  latex->SetTextSize(textSize);
  double step=textSize*1.3;

  gStyle->SetPadRightMargin(0.12);
  gStyle->SetOptStat(0);
  gStyle->SetFrameBorderMode(0);

  TCanvas *c1 = new TCanvas("", "", 500, 500);
  gStyle->SetPalette(1);
  gPad->SetFillColor(kWhite);
  gPad->SetLeftMargin(0.125);
  gPad->SetRightMargin(0.15);

  double yOffset=1.4;
  //if(iWindow > 0){
  // comb. background PRSR1
  hCosThetaPhi[iRapBin][iPTBin][iFrame][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
  hCosThetaPhi[iRapBin][iPTBin][iFrame][iWindow][0]->Draw("colz");
  if(iRapBin==0)
    latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  else
    latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));

  std::stringstream nameBG;
  if(nState == 6)
    nameBG << "Figures/cosThetaPhi_chic1_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
  else
    nameBG << "Figures/cosThetaPhi_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
  c1->SaveAs(nameBG.str().c_str());

  // comb. background PRSR2
  hCosThetaPhi[iRapBin][iPTBin][iFrame][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
  hCosThetaPhi[iRapBin][iPTBin][iFrame][iWindow][1]->Draw("colz");
  if(iRapBin==0)
    latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  else
    latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));

  std::stringstream nameBG1;
  if(nState == 6){
    nameBG1 << "Figures/cosThetaPhi_chic2_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
    c1->SaveAs(nameBG1.str().c_str());
  }

  // comb. background NPSR1
  hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
  hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][iWindow][0]->Draw("colz");

  if(iRapBin==0)
    latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  else
    latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));

  std::stringstream nameBGinNP;
  if(nState == 6)
    nameBGinNP << "Figures/cosThetaPhi_chic1_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
  else
    nameBGinNP << "Figures/cosThetaPhi_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
  c1->SaveAs(nameBGinNP.str().c_str());


  // comb. background NPSR2
  hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
  hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][iWindow][1]->Draw("colz");
  if(iRapBin==0)
    latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  else
    latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  std::stringstream nameBGinNP1;
  if(nState == 6){
    nameBGinNP1 << "Figures/cosThetaPhi_chic2_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
    c1->SaveAs(nameBGinNP1.str().c_str());
  }
  //}// if(iWindow)
  //else {
  // non prompt background SR1
  hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame][0]->GetYaxis()->SetTitleOffset(yOffset);
  hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame][0]->Draw("colz");
  if(iRapBin==0) 
    latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  else
    latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  std::stringstream nameNPBG;
  if(nState == 6)
    nameNPBG << "Figures/cosThetaPhiNPBG_chic1_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << ".pdf";
  else
    nameNPBG << "Figures/cosThetaPhiNPBG_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << ".pdf";
  c1->SaveAs(nameNPBG.str().c_str());

  // non prompt background SR2
  hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame][1]->GetYaxis()->SetTitleOffset(yOffset);
  hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame][1]->Draw("colz");

  if(iRapBin==0)
    latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  else
    latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));


  std::stringstream nameNPBG1;
  if(nState == 6){
    nameNPBG1 << "Figures/cosThetaPhiNPBG_chic2_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << ".pdf";
    c1->SaveAs(nameNPBG1.str().c_str());
  }

  // total background in PRSR1
  hTCosThetaPhi[iRapBin][iPTBin][iFrame][0]->GetYaxis()->SetTitleOffset(yOffset);
  hTCosThetaPhi[iRapBin][iPTBin][iFrame][0]->Draw("colz");
  if(iRapBin==0) 
    latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  else
    latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  std::stringstream name;
  if(nState == 6)
    name << "Figures/cosThetaPhi_chic1_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_total.pdf";
  else
    name << "Figures/cosThetaPhi_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_total.pdf";
  c1->SaveAs(name.str().c_str());

  // total background in PRSR2
  hTCosThetaPhi[iRapBin][iPTBin][iFrame][1]->GetYaxis()->SetTitleOffset(yOffset);
  hTCosThetaPhi[iRapBin][iPTBin][iFrame][1]->Draw("colz");
  if(iRapBin==0)
    latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  else
    latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                   onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                   onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
  std::stringstream name1;
  if(nState == 6){
    name1 << "Figures/cosThetaPhi_chic2_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_total.pdf";
    c1->SaveAs(name1.str().c_str());
  }
  //}// else

  delete c1;
}

//===========================
void LoadHistos(int iRapBin, int iPTBin, int nState){

  std::cout << "rap " << iRapBin+1 << " pt " << iPTBin+1 << std::endl;
  std::stringstream name, name1;
  if(nState == 6){
    name << "tmpFiles/data_chic1_rap" << iRapBin+1 << "_pT" << iPTBin+1 << ".root";
    name1 << "tmpFiles/data_chic2_rap" << iRapBin+1 << "_pT" << iPTBin+1 << ".root";
  } else {
    name << "tmpFiles/data_Psi" << nState-3 << "_rap" << iRapBin+1 << "_pT" << iPTBin+1 << ".root";
    name1 << "tmpFiles/data_Psi" << nState-3 << "_rap" << iRapBin+1 << "_pT" << iPTBin+1 << ".root";
  }

  std::cout << "open files" << std::endl;
  TFile *fIn = new TFile(name.str().c_str());
  TFile* fIn1 = new TFile(name1.str().c_str());

  // get mean pt and y
  std::cout << "get mean pt and y" << std::endl;
  hMeanPt[iRapBin][iPTBin][0] = (TH1F*)fIn->Get("mean_pT");
  meanPt[iRapBin][iPTBin][0] = hMeanPt[iRapBin][iPTBin][0]->GetBinContent(1);
  hMeanY[iRapBin][iPTBin][0] = (TH1F*)fIn->Get("mean_y");
  meanRap[iRapBin][iPTBin][0] = hMeanY[iRapBin][iPTBin][0]->GetBinContent(1);
  hMeanPt[iRapBin][iPTBin][1] = (TH1F*)fIn1->Get("mean_pT");
  meanPt[iRapBin][iPTBin][1] = hMeanPt[iRapBin][iPTBin][1]->GetBinContent(1);
  hMeanY[iRapBin][iPTBin][1] = (TH1F*)fIn1->Get("mean_y");
  meanRap[iRapBin][iPTBin][1] = hMeanY[iRapBin][iPTBin][1]->GetBinContent(1);

  // get number of events in PRSR and NPSR
  std::cout << "get number of events" << std::endl;
  events_PRSR[iRapBin][iPTBin][0] = (TH1F*) fIn->Get("events_promptSR");
  events_NPSR[iRapBin][iPTBin][0] = (TH1F*) fIn->Get("events_nonpromptSR");
  evtPinPRSR [iRapBin][iPTBin][0] = (int)events_PRSR[iRapBin][iPTBin][0]->GetBinContent(1);
  evtNPinPRSR[iRapBin][iPTBin][0] = (int)events_PRSR[iRapBin][iPTBin][0]->GetBinContent(2);
  evtBGinPRSR[iRapBin][iPTBin][0] = (int)events_PRSR[iRapBin][iPTBin][0]->GetBinContent(3);
  evtPinNPSR [iRapBin][iPTBin][0] = (int)events_NPSR[iRapBin][iPTBin][0]->GetBinContent(1);
  evtNPinNPSR[iRapBin][iPTBin][0] = (int)events_NPSR[iRapBin][iPTBin][0]->GetBinContent(2);
  evtBGinNPSR[iRapBin][iPTBin][0] = (int)events_NPSR[iRapBin][iPTBin][0]->GetBinContent(3);
  events_PRSR[iRapBin][iPTBin][1] = (TH1F*) fIn1->Get("events_promptSR");
  events_NPSR[iRapBin][iPTBin][1] = (TH1F*) fIn1->Get("events_nonpromptSR");
  evtPinPRSR [iRapBin][iPTBin][1] = (int)events_PRSR[iRapBin][iPTBin][1]->GetBinContent(1);
  evtNPinPRSR[iRapBin][iPTBin][1] = (int)events_PRSR[iRapBin][iPTBin][1]->GetBinContent(2);
  evtBGinPRSR[iRapBin][iPTBin][1] = (int)events_PRSR[iRapBin][iPTBin][1]->GetBinContent(3);
  evtPinNPSR [iRapBin][iPTBin][1] = (int)events_NPSR[iRapBin][iPTBin][1]->GetBinContent(1);
  evtNPinNPSR[iRapBin][iPTBin][1] = (int)events_NPSR[iRapBin][iPTBin][1]->GetBinContent(2);
  evtBGinNPSR[iRapBin][iPTBin][1] = (int)events_NPSR[iRapBin][iPTBin][1]->GetBinContent(3);

  if(nState == 6){
    std::cout << "get fractions from file" << std::endl;
    hfBGinPRSR[iRapBin][iPTBin][0] = (TH1F*)fIn->Get("comb_background_fraction");
    hfNPinPRSR[iRapBin][iPTBin][0] = (TH1F*)fIn->Get("nonprompt_background_fraction");
    hfBGinPRSR[iRapBin][iPTBin][1] = (TH1F*)fIn1->Get("comb_background_fraction");
    hfNPinPRSR[iRapBin][iPTBin][1] = (TH1F*)fIn1->Get("nonprompt_background_fraction");
    fracBGinPRSR[iRapBin][iPTBin][0] = hfBGinPRSR[iRapBin][iPTBin][0]->GetBinContent(1);
    fracNPinPRSR[iRapBin][iPTBin][0] = hfNPinPRSR[iRapBin][iPTBin][0]->GetBinContent(1);
    fracBGinPRSR[iRapBin][iPTBin][1] = hfBGinPRSR[iRapBin][iPTBin][1]->GetBinContent(1);
    fracNPinPRSR[iRapBin][iPTBin][1] = hfNPinPRSR[iRapBin][iPTBin][1]->GetBinContent(1);
  } else {
    fracBGinPRSR[iRapBin][iPTBin][0] = (double)evtBGinPRSR[iRapBin][iPTBin][0] /
      ( (double)evtPinPRSR [iRapBin][iPTBin][0] + (double)evtNPinPRSR[iRapBin][iPTBin][0] + (double)evtBGinPRSR[iRapBin][iPTBin][0] ) ;
    fracNPinPRSR[iRapBin][iPTBin][0] = (double)evtNPinPRSR[iRapBin][iPTBin][0] /
      ( (double)evtPinPRSR [iRapBin][iPTBin][0] + (double)evtNPinPRSR[iRapBin][iPTBin][0] + (double)evtBGinPRSR[iRapBin][iPTBin][0] ) ;
    fracBGinPRSR[iRapBin][iPTBin][1] = fracBGinPRSR[iRapBin][iPTBin][0];
    fracNPinPRSR[iRapBin][iPTBin][1] = fracNPinPRSR[iRapBin][iPTBin][0];
  }

  std::cout << "loop through frames" << std::endl;
  for(int iFrame = 0; iFrame < onia::kNbFrames; iFrame++){
    std::cout << "----------------" << onia::frameLabel[iFrame] << " frame " << "----------------" << std::endl;
    std::stringstream nameL, nameR, nameBG, nameBGinNPL, nameBGinNPR, nameNPB;
    nameL << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
    nameR << "hBG_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
    nameBGinNPL << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_L";
    nameBGinNPR << "hBGinNP_cosThetaPhi_" << onia::frameLabel[iFrame] << "_R";
    nameNPB << "hNPBG_cosThetaPhi_" << onia::frameLabel[iFrame];
    nameBG << "background_costhphi" << onia::frameLabel[iFrame];
    std::stringstream nameL1, nameR1, nameBG1, nameBGinNPL1, nameBGinNPR1, nameNPB1;
    nameL1 << "hCosThetaPhi_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" << iPTBin+1 << "_L";
    nameR1 << "hCosThetaPhi_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" << iPTBin+1 << "_R";
    nameBGinNPL1 << "hCosThetaPhi_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" << iPTBin+1 << "_highct_L";
    nameBGinNPR1 << "hCosThetaPhi_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" << iPTBin+1 << "_highct_R";
    nameNPB1 << "hCosThetaPhiNPBG_" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" << iPTBin+1;
    nameBG1 << "hTCosThetaPhi" << onia::frameLabel[iFrame] << "_rap" << iRapBin+1 << "_pT" << iPTBin+1;
    
    hCosThetaPhi[iRapBin][iPTBin][iFrame][L][0] = (TH2F *) fIn->Get(nameL.str().c_str());
    hCosThetaPhi[iRapBin][iPTBin][iFrame][L][1] = (TH2F *) fIn1->Get(nameL.str().c_str());
    hCosThetaPhi[iRapBin][iPTBin][iFrame][L][0]->SetName(nameL1.str().c_str());
    hCosThetaPhi[iRapBin][iPTBin][iFrame][L][1]->SetName(nameL1.str().c_str());
    hCosThetaPhi[iRapBin][iPTBin][iFrame][R][0] = (TH2F *) fIn->Get(nameR.str().c_str());
    hCosThetaPhi[iRapBin][iPTBin][iFrame][R][1] = (TH2F *) fIn1->Get(nameR.str().c_str());
    hCosThetaPhi[iRapBin][iPTBin][iFrame][R][0]->SetName(nameR1.str().c_str());
    hCosThetaPhi[iRapBin][iPTBin][iFrame][R][1]->SetName(nameR1.str().c_str());
    
    hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][L][0] = (TH2F *) fIn->Get(nameBGinNPL.str().c_str());
    hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][L][1] = (TH2F *) fIn1->Get(nameBGinNPL.str().c_str());
    hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][L][0]->SetName(nameBGinNPL1.str().c_str());
    hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][L][1]->SetName(nameBGinNPL1.str().c_str());
    hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][R][0] = (TH2F *) fIn->Get(nameBGinNPR.str().c_str());
    hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][R][1] = (TH2F *) fIn1->Get(nameBGinNPR.str().c_str());
    hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][R][0]->SetName(nameBGinNPR1.str().c_str());
    hCosThetaPhiHighct[iRapBin][iPTBin][iFrame][R][1]->SetName(nameBGinNPR1.str().c_str());
    
    hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame][0] = (TH2F *) fIn->Get(nameNPB.str().c_str());
    hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame][1] = (TH2F *) fIn1->Get(nameNPB.str().c_str());
    hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame][0]->SetName(nameNPB1.str().c_str());
    hCosThetaPhiNPBG[iRapBin][iPTBin][iFrame][1]->SetName(nameNPB1.str().c_str());
    
    hTCosThetaPhi[iRapBin][iPTBin][iFrame][0] = (TH2F *) fIn->Get(nameBG.str().c_str());
    hTCosThetaPhi[iRapBin][iPTBin][iFrame][1] = (TH2F *) fIn1->Get(nameBG.str().c_str());
    hTCosThetaPhi[iRapBin][iPTBin][iFrame][0]->SetName(nameBG1.str().c_str());
    hTCosThetaPhi[iRapBin][iPTBin][iFrame][1]->SetName(nameBG1.str().c_str());
    
    if(iFrame==2){
      std::cout << "PrintBin rap " << iRapBin+1 << " pt " << iPTBin+1 << " nBinsCosthTBG: " <<
        hTCosThetaPhi[iRapBin][iPTBin][iFrame][0]->GetNbinsX() << " nBinsPhiTBG:   "<<
        hTCosThetaPhi[iRapBin][iPTBin][iFrame][0]->GetNbinsY() << std::endl;
    }

  }

}
