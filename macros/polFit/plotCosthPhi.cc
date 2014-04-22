#include "../../interface/rootIncludes.inc"

#include <string>
#include <iostream>
#include <sstream>
#include "commonVar.h"

using namespace std;

void plotCosthPhi() {

	gStyle->SetOptStat(0);

	char fileName[500];
	TFile *infile[5];

	TH1D  *PDF_cth_HX_0[5];
	TH1D  *PDF_ph_HX_0[5];

	sprintf(fileName,"/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/Psi/Data/Psi1S_ctauScen5_FracLSB-1_26May2013_1FracBg_1BgModel_1RhoAll/Figures/costh_phi_HX_zeroCurve_Psi1S_rap1_pT7.root");
	infile[0] = new TFile(fileName,"R");
	PDF_cth_HX_0[0] = (TH1D*) infile[0] -> Get("PDF_cth_HX_0");
	PDF_ph_HX_0[0] = (TH1D*) infile[0] -> Get("PDF_ph_HX_0");

	sprintf(fileName,"/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/Psi/Data/Psi1S_ctauScen5_FracLSB-1_26May2013_1FracBg_1BgModel_1RhoAll/Figures/costh_phi_HX_zeroCurve_Psi1S_rap1_pT8.root");
	infile[1] = new TFile(fileName,"R");
	PDF_cth_HX_0[1] = (TH1D*) infile[1] -> Get("PDF_cth_HX_0");
	PDF_ph_HX_0[1] = (TH1D*) infile[1] -> Get("PDF_ph_HX_0");

	sprintf(fileName,"/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/Psi/Data/Psi1S_ctauScen5_FracLSB-1_26May2013_1FracBg_1BgModel_1RhoAll/Figures/costh_phi_HX_zeroCurve_Psi1S_rap1_pT9.root");
	infile[2] = new TFile(fileName,"R");
	PDF_cth_HX_0[2] = (TH1D*) infile[2] -> Get("PDF_cth_HX_0");
	PDF_ph_HX_0[2] = (TH1D*) infile[2] -> Get("PDF_ph_HX_0");

	sprintf(fileName,"/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/Psi/Data/Psi1S_ctauScen5_FracLSB-1_26May2013_1FracBg_1BgModel_1RhoAll/Figures/costh_phi_HX_zeroCurve_Psi1S_rap1_pT10.root");
	infile[3] = new TFile(fileName,"R");
	PDF_cth_HX_0[3] = (TH1D*) infile[3] -> Get("PDF_cth_HX_0");
	PDF_ph_HX_0[3] = (TH1D*) infile[3] -> Get("PDF_ph_HX_0");
	
	sprintf(fileName,"/afs/ihep.ac.cn/users/z/zhangll/fs/work/polarization/PsiPol2011/Psi/Data/Psi1S_ctauScen5_FracLSB-1_26May2013_1FracBg_1BgModel_1RhoAll/Figures/costh_phi_HX_zeroCurve_Psi1S_rap1_pT11.root");
	infile[4] = new TFile(fileName,"R");
	PDF_cth_HX_0[4] = (TH1D*) infile[4] -> Get("PDF_cth_HX_0");
	PDF_ph_HX_0[4] = (TH1D*) infile[4] -> Get("PDF_ph_HX_0");


	PDF_cth_HX_0[0]->Scale(1./PDF_cth_HX_0[0]->Integral());
	PDF_cth_HX_0[1]->Scale(1./PDF_cth_HX_0[1]->Integral());
	PDF_cth_HX_0[2]->Scale(1./PDF_cth_HX_0[2]->Integral());
	PDF_cth_HX_0[3]->Scale(1./PDF_cth_HX_0[3]->Integral());
	PDF_cth_HX_0[4]->Scale(1./PDF_cth_HX_0[4]->Integral());

	PDF_ph_HX_0[0]->Scale(1./PDF_ph_HX_0[0]->Integral());
	PDF_ph_HX_0[1]->Scale(1./PDF_ph_HX_0[1]->Integral());
	PDF_ph_HX_0[2]->Scale(1./PDF_ph_HX_0[2]->Integral());
	PDF_ph_HX_0[3]->Scale(1./PDF_ph_HX_0[3]->Integral());
	PDF_ph_HX_0[4]->Scale(1./PDF_ph_HX_0[4]->Integral());
	
	//PDF_cth_HX_0[0]->Scale(1./PDF_cth_HX_0[0]->GetMaximum());
	//PDF_cth_HX_0[1]->Scale(1./PDF_cth_HX_0[1]->GetMaximum());
	//PDF_cth_HX_0[2]->Scale(1./PDF_cth_HX_0[2]->GetMaximum());
	//PDF_cth_HX_0[3]->Scale(1./PDF_cth_HX_0[3]->GetMaximum());
	//PDF_cth_HX_0[4]->Scale(1./PDF_cth_HX_0[4]->GetMaximum());

	//PDF_ph_HX_0[0]->Scale(1./PDF_ph_HX_0[0]->GetMaximum());
	//PDF_ph_HX_0[1]->Scale(1./PDF_ph_HX_0[1]->GetMaximum());
	//PDF_ph_HX_0[2]->Scale(1./PDF_ph_HX_0[2]->GetMaximum());

	TCanvas* c3 = new TCanvas("c3", "c3", 10, 28, 588,563);
	c3->Range(-1.370833,-114.012,1.029167,745.8285);
	c3->SetFillColor(0);
	c3->SetBorderMode(0);
	c3->SetBorderSize(0);
	c3->SetLeftMargin(0.1545139);
	c3->SetRightMargin(0.01215278);
	c3->SetTopMargin(0.01841621);
	c3->SetBottomMargin(0.1325967);
	c3->SetFrameBorderMode(0);

	PDF_cth_HX_0[0]->GetXaxis()->SetTitle("cos#vartheta_{HX}");
	PDF_cth_HX_0[0]->GetXaxis()->SetLabelOffset(0.028);
	PDF_cth_HX_0[0]->GetXaxis()->SetTitleSize(0.05);
	PDF_cth_HX_0[0]->GetXaxis()->SetTickLength(-0.03);
	PDF_cth_HX_0[0]->GetXaxis()->SetTitleOffset(1.20);
	PDF_cth_HX_0[0]->GetYaxis()->SetTitle("");//event PDF [a.u.]");
	PDF_cth_HX_0[0]->GetYaxis()->SetLabelOffset(0.032);
	PDF_cth_HX_0[0]->GetYaxis()->SetTitleSize(0.05);
	PDF_cth_HX_0[0]->GetYaxis()->SetTickLength(-0.03);
	PDF_cth_HX_0[0]->GetYaxis()->SetTitleOffset(1.55);
	PDF_cth_HX_0[0]->SetMinimum(0.);
	PDF_cth_HX_0[0]->SetLineColor(kRed);
	PDF_cth_HX_0[0]->SetLineStyle(1);
	PDF_cth_HX_0[0]->SetLineWidth(2);
	PDF_cth_HX_0[0]->SetMaximum(0.022);
	PDF_cth_HX_0[0]->Draw("L");

	PDF_cth_HX_0[1]->SetLineColor(kBlue);
	PDF_cth_HX_0[1]->Draw("L same");
	PDF_cth_HX_0[2]->SetLineColor(kBlack);
	PDF_cth_HX_0[2]->Draw("L same");
	PDF_cth_HX_0[3]->SetLineColor(kMagenta);
	PDF_cth_HX_0[3]->Draw("L same");
	PDF_cth_HX_0[4]->SetLineColor(kOrange);
	PDF_cth_HX_0[4]->Draw("L same");

	TLegend* plotLegend2;

	plotLegend2=new TLegend(0.3,0.15,0.8,0.3);//0.75,0.95);
	plotLegend2->SetFillColor(kWhite);
	//plotLegend2->SetTextFont(72);
	//plotLegend2->SetTextSize(0.035);
	plotLegend2->SetBorderSize(1);

	plotLegend2->AddEntry(PDF_cth_HX_0[0],"|y| < 0.6, 22 < p_{T} < 25 GeV");
	plotLegend2->AddEntry(PDF_cth_HX_0[1],"|y| < 0.6, 25 < p_{T} < 30 GeV");
	plotLegend2->AddEntry(PDF_cth_HX_0[2],"|y| < 0.6, 30 < p_{T} < 35 GeV");
	plotLegend2->AddEntry(PDF_cth_HX_0[3],"|y| < 0.6, 35 < p_{T} < 40 GeV");
	plotLegend2->AddEntry(PDF_cth_HX_0[4],"|y| < 0.6, 40 < p_{T} < 50 GeV");
	plotLegend2->Draw("same");

	c3->Print("costh_0curve.pdf");


	PDF_ph_HX_0[0]->GetXaxis()->SetTitle("#varphi_{HX}");
	PDF_ph_HX_0[0]->GetXaxis()->SetLabelOffset(0.028);
	PDF_ph_HX_0[0]->GetXaxis()->SetTitleSize(0.05);
	PDF_ph_HX_0[0]->GetXaxis()->SetTickLength(-0.03);
	PDF_ph_HX_0[0]->GetXaxis()->SetTitleOffset(1.20);
	PDF_ph_HX_0[0]->GetYaxis()->SetTitle("");//event PDF [a.u.]");
	PDF_ph_HX_0[0]->GetYaxis()->SetLabelOffset(0.032);
	PDF_ph_HX_0[0]->GetYaxis()->SetTitleSize(0.05);
	PDF_ph_HX_0[0]->GetYaxis()->SetTickLength(-0.03);
	PDF_ph_HX_0[0]->GetYaxis()->SetTitleOffset(1.55);
	PDF_ph_HX_0[0]->SetMinimum(0.011);
	PDF_ph_HX_0[0]->SetLineColor(kRed);
	PDF_ph_HX_0[0]->SetLineStyle(1);
	PDF_ph_HX_0[0]->SetLineWidth(2);
	PDF_ph_HX_0[0]->SetMaximum(0.016);
	PDF_ph_HX_0[0]->Draw("L");

	PDF_ph_HX_0[1]->SetLineColor(kBlue);
	PDF_ph_HX_0[1]->Draw("L same");
	PDF_ph_HX_0[2]->SetLineColor(kBlack);
	PDF_ph_HX_0[2]->Draw("L same");
	PDF_ph_HX_0[3]->SetLineColor(kMagenta);
	PDF_ph_HX_0[3]->Draw("L same");
	PDF_ph_HX_0[4]->SetLineColor(kOrange);
	PDF_ph_HX_0[4]->Draw("L same");

	//TLegend* plotLegend2;

	plotLegend2=new TLegend(0.2,0.15,0.75,0.3);//0.75,0.95);
	plotLegend2->SetFillColor(kWhite);
	//plotLegend2->SetTextFont(72);
	//plotLegend2->SetTextSize(0.035);
	plotLegend2->SetBorderSize(1);

	plotLegend2->AddEntry(PDF_ph_HX_0[0],"|y| < 0.6, 22 < p_{T} < 25 GeV");
	plotLegend2->AddEntry(PDF_ph_HX_0[1],"|y| < 0.6, 25 < p_{T} < 30 GeV");
	plotLegend2->AddEntry(PDF_ph_HX_0[2],"|y| < 0.6, 30 < p_{T} < 35 GeV");
	plotLegend2->AddEntry(PDF_ph_HX_0[3],"|y| < 0.6, 35 < p_{T} < 40 GeV");
	plotLegend2->AddEntry(PDF_ph_HX_0[4],"|y| < 0.6, 40 < p_{T} < 50 GeV");
	plotLegend2->Draw("same");
	

	c3->Print("phi_0curve.pdf");


}
