#include <iostream>
#include <string>
#include <sstream>
using namespace std;

#include "rootIncludes.inc"
#include "commonVar.h"

#include "TH1.h"
#include "TH3.h"

enum{L,R};
const char *bgLabel[2] = {"L", "R"};
TH3D *hMassRapPt[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH3D *hJpsiMassRapPt[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH3D *hMassRapPtHighct[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH3D *hJpsiMassRapPtHighct[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH3D *hMassRapPtNPBG[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];
TH3D *hTMassRapPt[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];
TH1D *hMass[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH1D *hJpsiMass[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH1D *hMassHighct[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH1D *hJpsiMassHighct[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH1D *hMassNPBG[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];
TH1D *hTMass[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];
TH1D *hPt[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH1D *hJpsiPt[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH1D *hPtHighct[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH1D *hJpsiPtHighct[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH1D *hPtNPBG[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];
TH1D *hTPt[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];
TH1D *hRap[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH1D *hJpsiRap[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH1D *hRapHighct[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH1D *hJpsiRapHighct[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2][2];
TH1D *hRapNPBG[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];
TH1D *hTRap[onia::kNbRapForPTBins][onia::kNbPTMaxBins][2];

void LoadHistos(int iRapBin, int iPTBin, int nState);
void PlotHistos(int iRapBin, int iPTBin, int nState, int iWindow);

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
            PlotHistos(iRap, iPT, nState, L);
            PlotHistos(iRap, iPT, nState, R);
        }
    }

    return 0;

}

//===========================
void PlotHistos(int iRapBin, int iPTBin, int nState, int iWindow){

    std::cout << "-------------------------------\n Plot histos" << std::endl;
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
    // comb. background chic1
    // mass
    hMass[iRapBin][iPTBin][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
    hMass[iRapBin][iPTBin][iWindow][0]->Draw("");
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream massBG;
    if(nState==6)
        massBG << "Figures/mass_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
    else
        massBG << "Figures/mass_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
    c1->SaveAs(massBG.str().c_str());
    //--
    hJpsiMass[iRapBin][iPTBin][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
    hJpsiMass[iRapBin][iPTBin][iWindow][0]->Draw("");
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState==6){
        std::stringstream massJpsiBG;
        massJpsiBG << "Figures/mass_JpsiBG_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(massJpsiBG.str().c_str());
    }

    // rapidity
    hRap[iRapBin][iPTBin][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
    hRap[iRapBin][iPTBin][iWindow][0]->Draw("");
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream rapBG;
    if(nState == 6)
        rapBG << "Figures/rap_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
    else
        rapBG << "Figures/rap_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
    c1->SaveAs(rapBG.str().c_str());
    //--
    hJpsiRap[iRapBin][iPTBin][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
    hJpsiRap[iRapBin][iPTBin][iWindow][0]->Draw("");
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream rapJpsiBG;
        rapJpsiBG << "Figures/rap_JpsiBG_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(rapJpsiBG.str().c_str());
    }

    // pt
    hPt[iRapBin][iPTBin][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
    hPt[iRapBin][iPTBin][iWindow][0]->Draw("");
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream ptBG;
    if(nState == 6)
        ptBG << "Figures/pT_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
    else
        ptBG << "Figures/pT_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
    c1->SaveAs(ptBG.str().c_str());
    //---
    hJpsiPt[iRapBin][iPTBin][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
    hJpsiPt[iRapBin][iPTBin][iWindow][0]->Draw("");
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream ptJpsiBG;
        ptJpsiBG << "Figures/pT_JpsiBG_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(ptJpsiBG.str().c_str());
    }

    // comb. background chic2
    // mass
    hMass[iRapBin][iPTBin][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
    hMass[iRapBin][iPTBin][iWindow][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream massBG1;
        massBG1 << "Figures/mass_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(massBG1.str().c_str());
    }
    //--
    hJpsiMass[iRapBin][iPTBin][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
    hJpsiMass[iRapBin][iPTBin][iWindow][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream massBG1;
        massBG1 << "Figures/mass_JpsiBG_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(massBG1.str().c_str());
    }

    // rapidity
    hRap[iRapBin][iPTBin][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
    hRap[iRapBin][iPTBin][iWindow][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream rapBG1;
        rapBG1 << "Figures/rap_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(rapBG1.str().c_str());
    }
    //---
    hJpsiRap[iRapBin][iPTBin][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
    hJpsiRap[iRapBin][iPTBin][iWindow][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
    std::stringstream rapBG1;
        rapBG1 << "Figures/rap_JpsiBG_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(rapBG1.str().c_str());
    }

    // pt
    hPt[iRapBin][iPTBin][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
    hPt[iRapBin][iPTBin][iWindow][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream ptBG1;
        ptBG1 << "Figures/pT_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(ptBG1.str().c_str());
    }
    //---
    hJpsiPt[iRapBin][iPTBin][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
    hJpsiPt[iRapBin][iPTBin][iWindow][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream ptBG1;
        ptBG1 << "Figures/pT_JpsiBG_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(ptBG1.str().c_str());
    }

    // non prompt comb. background chic1
    // mass
    hMassHighct[iRapBin][iPTBin][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
    hMassHighct[iRapBin][iPTBin][iWindow][0]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream massBGinNP;
    if(nState == 6)
        massBGinNP << "Figures/mass_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
    else
        massBGinNP << "Figures/mass_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
    c1->SaveAs(massBGinNP.str().c_str());
    //---
    hJpsiMassHighct[iRapBin][iPTBin][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
    hJpsiMassHighct[iRapBin][iPTBin][iWindow][0]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream massJpsiBGinNP;
        massJpsiBGinNP << "Figures/mass_JpsiBG_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(massJpsiBGinNP.str().c_str());
    }    

    // rapidity
    hRapHighct[iRapBin][iPTBin][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
    hRapHighct[iRapBin][iPTBin][iWindow][0]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream rapBGinNP;
    if(nState == 6)
        rapBGinNP << "Figures/rap_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
    else
        rapBGinNP << "Figures/rap_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
    c1->SaveAs(rapBGinNP.str().c_str());
    //---
    hJpsiRapHighct[iRapBin][iPTBin][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
    hJpsiRapHighct[iRapBin][iPTBin][iWindow][0]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream rapJpsiBGinNP;
        rapJpsiBGinNP << "Figures/rap_JpsiBG_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(rapJpsiBGinNP.str().c_str());
    }

    //pt
    hPtHighct[iRapBin][iPTBin][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
    hPtHighct[iRapBin][iPTBin][iWindow][0]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream ptBGinNP;
    if(nState == 6)
        ptBGinNP << "Figures/pT_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
    else
        ptBGinNP << "Figures/pT_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
    c1->SaveAs(ptBGinNP.str().c_str());
    //---
    hJpsiPtHighct[iRapBin][iPTBin][iWindow][0]->GetYaxis()->SetTitleOffset(yOffset);
    hJpsiPtHighct[iRapBin][iPTBin][iWindow][0]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream ptJpsiBGinNP;
        ptJpsiBGinNP << "Figures/pT_JpsiBG_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(ptJpsiBGinNP.str().c_str());
    }

    // non prompt comb. background chic1
    // mass
    hMassHighct[iRapBin][iPTBin][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
    hMassHighct[iRapBin][iPTBin][iWindow][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
    std::stringstream massBGinNP1;
        massBGinNP1 << "Figures/mass_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(massBGinNP1.str().c_str());
    }
    //---
    hMassHighct[iRapBin][iPTBin][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
    hMassHighct[iRapBin][iPTBin][iWindow][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream massBGinNP1;
        massBGinNP1 << "Figures/mass_JpsiBG_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(massBGinNP1.str().c_str());
    }

    // rapidity
    hRapHighct[iRapBin][iPTBin][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
    hRapHighct[iRapBin][iPTBin][iWindow][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream rapBGinNP1;
        rapBGinNP1 << "Figures/rap_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(rapBGinNP1.str().c_str());
    }
    //---
    hJpsiRapHighct[iRapBin][iPTBin][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
    hJpsiRapHighct[iRapBin][iPTBin][iWindow][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream rapBGinNP1;
        rapBGinNP1 << "Figures/rap_JpsiBG_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(rapBGinNP1.str().c_str());
    }

    // pt
    hPtHighct[iRapBin][iPTBin][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
    hPtHighct[iRapBin][iPTBin][iWindow][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream ptBGinNP1;
        ptBGinNP1 << "Figures/pT_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(ptBGinNP1.str().c_str());
    }
    //---
    hJpsiPtHighct[iRapBin][iPTBin][iWindow][1]->GetYaxis()->SetTitleOffset(yOffset);
    hJpsiPtHighct[iRapBin][iPTBin][iWindow][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    if(nState == 6){
        std::stringstream ptBGinNP1;
        ptBGinNP1 << "Figures/pT_JpsiBG_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_highct_" << bgLabel[iWindow] << ".pdf";
        c1->SaveAs(ptBGinNP1.str().c_str());
    }

    // non prompt background SR1
    hMassNPBG[iRapBin][iPTBin][0]->GetYaxis()->SetTitleOffset(yOffset);
    hMassNPBG[iRapBin][iPTBin][0]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream massNPBG;
    if(nState == 6)
        massNPBG << "Figures/massNPBG_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << ".pdf";
    else
        massNPBG << "Figures/massNPBG_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << ".pdf";
    c1->SaveAs(massNPBG.str().c_str());

    hRapNPBG[iRapBin][iPTBin][0]->GetYaxis()->SetTitleOffset(yOffset);
    hRapNPBG[iRapBin][iPTBin][0]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream rapNPBG;
    if(nState == 6)
        rapNPBG << "Figures/rapNPBG_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << ".pdf";
    else
        rapNPBG << "Figures/rapNPBG_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << ".pdf";
    c1->SaveAs(rapNPBG.str().c_str());

    hPtNPBG[iRapBin][iPTBin][0]->GetYaxis()->SetTitleOffset(yOffset);
    hPtNPBG[iRapBin][iPTBin][0]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream ptNPBG;
    if(nState == 6)
        ptNPBG << "Figures/pTNPBG_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << ".pdf";
    else
        ptNPBG << "Figures/pTNPBG_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << ".pdf";
    c1->SaveAs(ptNPBG.str().c_str());

    // non prompt background SR2
    hMassNPBG[iRapBin][iPTBin][1]->GetYaxis()->SetTitleOffset(yOffset);
    hMassNPBG[iRapBin][iPTBin][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream massNPBG1;
    if(nState == 6){
        massNPBG1 << "Figures/massNPBG_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << ".pdf";
        c1->SaveAs(massNPBG1.str().c_str());
    }

    hRapNPBG[iRapBin][iPTBin][1]->GetYaxis()->SetTitleOffset(yOffset);
    hRapNPBG[iRapBin][iPTBin][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream rapNPBG1;
    if(nState == 6){
        rapNPBG1 << "Figures/rapNPBG_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << ".pdf";
        c1->SaveAs(rapNPBG1.str().c_str());
    }

    hPtNPBG[iRapBin][iPTBin][1]->GetYaxis()->SetTitleOffset(yOffset);
    hPtNPBG[iRapBin][iPTBin][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream ptNPBG1;
    if(nState == 6){
        ptNPBG1 << "Figures/pTNPBG_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << ".pdf";
        c1->SaveAs(ptNPBG1.str().c_str());
    }

    // total background in PRSR1
    hTMass[iRapBin][iPTBin][0]->GetYaxis()->SetTitleOffset(yOffset);
    hTMass[iRapBin][iPTBin][0]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream mass;
    if(nState == 6)
        mass << "Figures/mass_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_total.pdf";
    else
        mass << "Figures/mass_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_total.pdf";
    c1->SaveAs(mass.str().c_str());

    hTRap[iRapBin][iPTBin][0]->GetYaxis()->SetTitleOffset(yOffset);
    hTRap[iRapBin][iPTBin][0]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream rap;
    if(nState == 6)
        rap << "Figures/rap_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_total.pdf";
    else
        rap << "Figures/rap_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_total.pdf";
    c1->SaveAs(rap.str().c_str());

    hTPt[iRapBin][iPTBin][0]->GetYaxis()->SetTitleOffset(yOffset);
    hTPt[iRapBin][iPTBin][0]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream pt;
    if(nState == 6)
        pt << "Figures/pT_chic1_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_total.pdf";
    else
        pt << "Figures/pT_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_total.pdf";
    c1->SaveAs(pt.str().c_str());

    // total background in PRSR2
    hTMass[iRapBin][iPTBin][1]->GetYaxis()->SetTitleOffset(yOffset);
    hTMass[iRapBin][iPTBin][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream mass1;
    if(nState == 6){
        mass1 << "Figures/mass_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_total.pdf";
        c1->SaveAs(mass1.str().c_str());
    }

    hTRap[iRapBin][iPTBin][1]->GetYaxis()->SetTitleOffset(yOffset);
    hTRap[iRapBin][iPTBin][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream rap1;
    if(nState == 6){
        rap1 << "Figures/rap_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_total.pdf";
        c1->SaveAs(rap1.str().c_str());
    }

    hTPt[iRapBin][iPTBin][1]->GetYaxis()->SetTitleOffset(yOffset);
    hTPt[iRapBin][iPTBin][1]->Draw();
    if(iRapBin==0)
        latex->DrawLatex(left,top,Form("|y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    else
        latex->DrawLatex(left,top,Form("%.1f < |y| < %.1f, %.1f < p_{T} < %.1f GeV",
                                       onia::rapForPTRange[iRapBin],onia::rapForPTRange[iRapBin+1],
                                       onia::pTRange[iRapBin+1][iPTBin],onia::pTRange[iRapBin+1][iPTBin+1]));
    std::stringstream pt1;
    if(nState == 6){
        pt1 << "Figures/pT_chic2_rap" << iRapBin+1 << "_pT" <<  iPTBin+1 << "_total.pdf";
        c1->SaveAs(pt1.str().c_str());
    }

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
    std::cout << name.str() << ", " << name1.str() << std::endl;
    TFile *fIn = new TFile(name.str().c_str());
    TFile *fIn1 = new TFile(name1.str().c_str());
    std::cout << "files opened" << std::endl;

    std::string nameL = "hBG_pTRapMass_L";
    std::string nameR = "hBG_pTRapMass_R";
    std::string nameJpsiL = "hBG_Jpsi_pTRapMass_L";
    std::string nameJpsiR = "hBG_Jpsi_pTRapMass_R";
    std::string nameBGinNPL = "hBG_pTRapMass_highct_L";
    std::string nameBGinNPR = "hBG_pTRapMass_highct_R";
    std::string nameJpsiBGinNPL = "hBG_Jpsi_pTRapMass_highct_L";
    std::string nameJpsiBGinNPR = "hBG_Jpsi_pTRapMass_highct_R";
    std::string nameNPB = "NPS_highct_pTrapMass";
    std::string nameBG = "background_pTrapMass";

    std::cout << "left sideband" << std::endl;
    hMassRapPt[iRapBin][iPTBin][L][0] = (TH3D *) fIn->Get(nameL.c_str());
    hPt[iRapBin][iPTBin][L][0] = hMassRapPt[iRapBin][iPTBin][L][0]->ProjectionX();
    hRap[iRapBin][iPTBin][L][0] = hMassRapPt[iRapBin][iPTBin][L][0]->ProjectionY();
    hMass[iRapBin][iPTBin][L][0] = hMassRapPt[iRapBin][iPTBin][L][0]->ProjectionZ();
    hMassRapPt[iRapBin][iPTBin][L][1] = (TH3D *) fIn1->Get(nameL.c_str());
    hPt[iRapBin][iPTBin][L][1] = hMassRapPt[iRapBin][iPTBin][L][1]->ProjectionX();
    hRap[iRapBin][iPTBin][L][1] = hMassRapPt[iRapBin][iPTBin][L][1]->ProjectionY();
    hMass[iRapBin][iPTBin][L][1] = hMassRapPt[iRapBin][iPTBin][L][1]->ProjectionZ();
    std::cout << "left sideband Jpsi bkg" << std::endl;
    hJpsiMassRapPt[iRapBin][iPTBin][L][0] = (TH3D *) fIn->Get(nameJpsiL.c_str());
    // std::cout << hJpsiMassRapPt[iRapBin][iPTBin][L][0] << ", " << iRapBin << " " << iPTBin << " " << L << std::endl;
    hJpsiPt[iRapBin][iPTBin][L][0] = hJpsiMassRapPt[iRapBin][iPTBin][L][0]->ProjectionX();
    hJpsiRap[iRapBin][iPTBin][L][0] = hJpsiMassRapPt[iRapBin][iPTBin][L][0]->ProjectionY();
    hJpsiMass[iRapBin][iPTBin][L][0] = hJpsiMassRapPt[iRapBin][iPTBin][L][0]->ProjectionZ();
    hJpsiMassRapPt[iRapBin][iPTBin][L][1] = (TH3D *) fIn1->Get(nameJpsiL.c_str());
    hJpsiPt[iRapBin][iPTBin][L][1] = hJpsiMassRapPt[iRapBin][iPTBin][L][1]->ProjectionX();
    hJpsiRap[iRapBin][iPTBin][L][1] = hJpsiMassRapPt[iRapBin][iPTBin][L][1]->ProjectionY();
    hJpsiMass[iRapBin][iPTBin][L][1] = hJpsiMassRapPt[iRapBin][iPTBin][L][1]->ProjectionZ();
    std::cout << "right sideband" << std::endl;
    hMassRapPt[iRapBin][iPTBin][R][0] = (TH3D *) fIn->Get(nameR.c_str());
    hPt[iRapBin][iPTBin][R][0] = hMassRapPt[iRapBin][iPTBin][R][0]->ProjectionX();
    hRap[iRapBin][iPTBin][R][0] = hMassRapPt[iRapBin][iPTBin][R][0]->ProjectionY();
    hMass[iRapBin][iPTBin][R][0] = hMassRapPt[iRapBin][iPTBin][R][0]->ProjectionZ();
    hMassRapPt[iRapBin][iPTBin][R][1] = (TH3D *) fIn1->Get(nameR.c_str());
    hPt[iRapBin][iPTBin][R][1] = hMassRapPt[iRapBin][iPTBin][R][1]->ProjectionX();
    hRap[iRapBin][iPTBin][R][1] = hMassRapPt[iRapBin][iPTBin][R][1]->ProjectionY();
    hMass[iRapBin][iPTBin][R][1] = hMassRapPt[iRapBin][iPTBin][R][1]->ProjectionZ();
    std::cout << "right sideband Jpsi bkg" << std::endl;
    hJpsiMassRapPt[iRapBin][iPTBin][R][0] = (TH3D *) fIn->Get(nameJpsiR.c_str());
    hJpsiPt[iRapBin][iPTBin][R][0] = hJpsiMassRapPt[iRapBin][iPTBin][R][0]->ProjectionX();
    hJpsiRap[iRapBin][iPTBin][R][0] = hJpsiMassRapPt[iRapBin][iPTBin][R][0]->ProjectionY();
    hJpsiMass[iRapBin][iPTBin][R][0] = hJpsiMassRapPt[iRapBin][iPTBin][R][0]->ProjectionZ();
    hJpsiMassRapPt[iRapBin][iPTBin][R][1] = (TH3D *) fIn1->Get(nameJpsiR.c_str());
    hJpsiPt[iRapBin][iPTBin][R][1] = hJpsiMassRapPt[iRapBin][iPTBin][R][1]->ProjectionX();
    hJpsiRap[iRapBin][iPTBin][R][1] = hJpsiMassRapPt[iRapBin][iPTBin][R][1]->ProjectionY();
    hJpsiMass[iRapBin][iPTBin][R][1] = hJpsiMassRapPt[iRapBin][iPTBin][R][1]->ProjectionZ();

    std::cout << "left sideband high ctau" << std::endl;
    hMassRapPtHighct[iRapBin][iPTBin][L][0] = (TH3D *) fIn->Get(nameBGinNPL.c_str());
    hPtHighct[iRapBin][iPTBin][L][0] = hMassRapPtHighct[iRapBin][iPTBin][L][0]->ProjectionX();
    hRapHighct[iRapBin][iPTBin][L][0] = hMassRapPtHighct[iRapBin][iPTBin][L][0]->ProjectionY();
    hMassHighct[iRapBin][iPTBin][L][0] = hMassRapPtHighct[iRapBin][iPTBin][L][0]->ProjectionZ();
    hMassRapPtHighct[iRapBin][iPTBin][L][1] = (TH3D *) fIn1->Get(nameBGinNPL.c_str());
    hPtHighct[iRapBin][iPTBin][L][1] = hMassRapPtHighct[iRapBin][iPTBin][L][1]->ProjectionX();
    hRapHighct[iRapBin][iPTBin][L][1] = hMassRapPtHighct[iRapBin][iPTBin][L][1]->ProjectionY();
    hMassHighct[iRapBin][iPTBin][L][1] = hMassRapPtHighct[iRapBin][iPTBin][L][1]->ProjectionZ();
    std::cout << "left sideband high ctau Jpsi bkg" << std::endl;
    hJpsiMassRapPtHighct[iRapBin][iPTBin][L][0] = (TH3D *) fIn->Get(nameJpsiBGinNPL.c_str());
    hJpsiPtHighct[iRapBin][iPTBin][L][0] = hJpsiMassRapPtHighct[iRapBin][iPTBin][L][0]->ProjectionX();
    hJpsiRapHighct[iRapBin][iPTBin][L][0] = hJpsiMassRapPtHighct[iRapBin][iPTBin][L][0]->ProjectionY();
    hJpsiMassHighct[iRapBin][iPTBin][L][0] = hJpsiMassRapPtHighct[iRapBin][iPTBin][L][0]->ProjectionZ();
    hJpsiMassRapPtHighct[iRapBin][iPTBin][L][1] = (TH3D *) fIn1->Get(nameJpsiBGinNPL.c_str());
    hJpsiPtHighct[iRapBin][iPTBin][L][1] = hJpsiMassRapPtHighct[iRapBin][iPTBin][L][1]->ProjectionX();
    hJpsiRapHighct[iRapBin][iPTBin][L][1] = hJpsiMassRapPtHighct[iRapBin][iPTBin][L][1]->ProjectionY();
    hJpsiMassHighct[iRapBin][iPTBin][L][1] = hJpsiMassRapPtHighct[iRapBin][iPTBin][L][1]->ProjectionZ();
    std::cout << "right sideband high ctau" << std::endl;
    hMassRapPtHighct[iRapBin][iPTBin][R][0] = (TH3D *) fIn->Get(nameBGinNPR.c_str());
    hPtHighct[iRapBin][iPTBin][R][0] = hMassRapPtHighct[iRapBin][iPTBin][R][0]->ProjectionX();
    hRapHighct[iRapBin][iPTBin][R][0] = hMassRapPtHighct[iRapBin][iPTBin][R][0]->ProjectionY();
    hMassHighct[iRapBin][iPTBin][R][0] = hMassRapPtHighct[iRapBin][iPTBin][R][0]->ProjectionZ();
    hMassRapPtHighct[iRapBin][iPTBin][R][1] = (TH3D *) fIn1->Get(nameBGinNPR.c_str());
    hPtHighct[iRapBin][iPTBin][R][1] = hMassRapPtHighct[iRapBin][iPTBin][R][1]->ProjectionX();
    hRapHighct[iRapBin][iPTBin][R][1] = hMassRapPtHighct[iRapBin][iPTBin][R][1]->ProjectionY();
    hMassHighct[iRapBin][iPTBin][R][1] = hMassRapPtHighct[iRapBin][iPTBin][R][1]->ProjectionZ();
    std::cout << "right sideband high ctau Jpsi bkg" << std::endl;
    hJpsiMassRapPtHighct[iRapBin][iPTBin][R][0] = (TH3D *) fIn->Get(nameJpsiBGinNPR.c_str());
    hJpsiPtHighct[iRapBin][iPTBin][R][0] = hJpsiMassRapPtHighct[iRapBin][iPTBin][R][0]->ProjectionX();
    hJpsiRapHighct[iRapBin][iPTBin][R][0] = hJpsiMassRapPtHighct[iRapBin][iPTBin][R][0]->ProjectionY();
    hJpsiMassHighct[iRapBin][iPTBin][R][0] = hJpsiMassRapPtHighct[iRapBin][iPTBin][R][0]->ProjectionZ();
    hJpsiMassRapPtHighct[iRapBin][iPTBin][R][1] = (TH3D *) fIn1->Get(nameJpsiBGinNPR.c_str());
    hJpsiPtHighct[iRapBin][iPTBin][R][1] = hJpsiMassRapPtHighct[iRapBin][iPTBin][R][1]->ProjectionX();
    hJpsiRapHighct[iRapBin][iPTBin][R][1] = hJpsiMassRapPtHighct[iRapBin][iPTBin][R][1]->ProjectionY();
    hJpsiMassHighct[iRapBin][iPTBin][R][1] = hJpsiMassRapPtHighct[iRapBin][iPTBin][R][1]->ProjectionZ();

    std::cout << "non prompt background" << std::endl;
    hMassRapPtNPBG[iRapBin][iPTBin][0] = (TH3D *) fIn->Get(nameNPB.c_str());
    hPtNPBG[iRapBin][iPTBin][0] = hMassRapPtNPBG[iRapBin][iPTBin][0]->ProjectionX();
    hRapNPBG[iRapBin][iPTBin][0] = hMassRapPtNPBG[iRapBin][iPTBin][0]->ProjectionY();
    hMassNPBG[iRapBin][iPTBin][0] = hMassRapPtNPBG[iRapBin][iPTBin][0]->ProjectionZ();
    hMassRapPtNPBG[iRapBin][iPTBin][1] = (TH3D *) fIn1->Get(nameNPB.c_str());
    hPtNPBG[iRapBin][iPTBin][1] = hMassRapPtNPBG[iRapBin][iPTBin][1]->ProjectionX();
    hRapNPBG[iRapBin][iPTBin][1] = hMassRapPtNPBG[iRapBin][iPTBin][1]->ProjectionY();
    hMassNPBG[iRapBin][iPTBin][1] = hMassRapPtNPBG[iRapBin][iPTBin][1]->ProjectionZ();

    std::cout << "total background" << std::endl;
    hTMassRapPt[iRapBin][iPTBin][0] = (TH3D *) fIn->Get(nameBG.c_str());
    hTPt[iRapBin][iPTBin][0] = hTMassRapPt[iRapBin][iPTBin][0]->ProjectionX();
    hTRap[iRapBin][iPTBin][0] = hTMassRapPt[iRapBin][iPTBin][0]->ProjectionY();
    hTMass[iRapBin][iPTBin][0] = hTMassRapPt[iRapBin][iPTBin][0]->ProjectionZ();
    hTMassRapPt[iRapBin][iPTBin][1] = (TH3D *) fIn1->Get(nameBG.c_str());
    hTPt[iRapBin][iPTBin][1] = hTMassRapPt[iRapBin][iPTBin][1]->ProjectionX();
    hTRap[iRapBin][iPTBin][1] = hTMassRapPt[iRapBin][iPTBin][1]->ProjectionY();
    hTMass[iRapBin][iPTBin][1] = hTMassRapPt[iRapBin][iPTBin][1]->ProjectionZ();
    std::cout << "histograms loaded" << std::endl;

}
