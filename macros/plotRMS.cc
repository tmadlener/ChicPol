#include "../interface/rootIncludes.inc"
#include "../interface/RooUtils.h"
#include "../interface/commonVar_Psi1S.h"

void plotRMS(){

	gStyle->SetPadBottomMargin(0.11);
	gStyle->SetPadLeftMargin(0.08); //0.12
	gStyle->SetPadRightMargin(0.02); //0.05
	gStyle->SetPadTopMargin(0.05); //0.05

	TCanvas *c1=new TCanvas("c1","");
	c1->SetTickx();
	c1->SetTicky();
	//c1->SetGridx();
	//c1->SetGridy();
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetOptFit(0);            
	gStyle->SetOptStat(1);
	gStyle->SetTitleFont(22);        
	gStyle->SetStatFont(22); 
	gStyle->SetStatColor(10);  
	gStyle->SetStatBorderSize(1);
	gStyle->SetLabelFont(22,"X");
	gStyle->SetLabelFont(22,"Y");
	gStyle->SetTitleXOffset(1.2);
	gStyle->SetTitleYOffset(1.2);    
	gStyle->SetHistLineWidth(2);
	gStyle->SetStatX(0.9);
	gStyle->SetStatY(0.9);
	gStyle->SetTitleX(0.15);
	gStyle->SetTitleY(0.96);

	char dataID[500];
  //sprintf(dataID,"SetOfCuts11_ctauScen0_FracLSB-1_newMLfit_4Mar2013");
	sprintf(dataID,"SetOfCuts11_ctauScen0_FracLSB-1_newMLfit_correctCtau_11April2013");

	/////////////////////// both 1S and 2S
	TFile *infile_1S = new TFile(
			Form("DataFiles/%s/Psi1S/Fit/parameter/evaluateCtau/PR_2.5sigma/graphFrac_2.5sigma.root",dataID),"R");
	TFile *infile_2S = new TFile(
			Form("DataFiles/%s/Psi2S/Fit/parameter/evaluateCtau/PR_2sigma/graphFrac_2.0sigma.root",dataID),"R");
	if(!infile_1S){cout<<"no input file"<<endl; return;}
	if(!infile_2S){cout<<"no input file"<<endl; return;}
	TGraphAsymmErrors *graph_sigmaP_1S[2], *graph_sigmaP_2S[4];
	graph_sigmaP_1S[0] = (TGraphAsymmErrors*)infile_1S->Get("graph_sigmaP_0");
	graph_sigmaP_1S[1] = (TGraphAsymmErrors*)infile_1S->Get("graph_sigmaP_1");
	graph_sigmaP_2S[0] = (TGraphAsymmErrors*)infile_2S->Get("graph_sigmaP_0");
	graph_sigmaP_2S[1] = (TGraphAsymmErrors*)infile_2S->Get("graph_sigmaP_1");
	graph_sigmaP_2S[2] = (TGraphAsymmErrors*)infile_2S->Get("graph_sigmaP_2");

	//remove first pT bin for 2S
	//graph_sigmaP_2S[0] -> RemovePoint(0);
	//graph_sigmaP_2S[1] -> RemovePoint(0);
	//graph_sigmaP_2S[2] -> RemovePoint(0);

	// Psi 2S 
	graph_sigmaP_2S[3] = (TGraphAsymmErrors*)graph_sigmaP_2S[0]->Clone("graph_sigmaP_3");
	int NumP =  graph_sigmaP_2S[0]->GetN();
	double sigmaRap0=0, sigmaRap1=0, sigmaRap2=0;
	cout<<"NumP: "<<NumP<<endl;
	for(int i=0; i<NumP; i++){
		double x0,y0,x1,y1,x2,y2;
		graph_sigmaP_2S[0]->GetPoint(i,x0,y0);
		graph_sigmaP_2S[1]->GetPoint(i,x1,y1);
		graph_sigmaP_2S[2]->GetPoint(i,x2,y2);
		double xNew = (x0+x1+x2)/3.;
		double yNew = (y0+y1+y2)/3.;
		graph_sigmaP_2S[3]->SetPoint(i,xNew,yNew);

		//if(i>0){
			cout<<"i: "<<i<<endl;
			sigmaRap0 += y0;
			sigmaRap1 += y1;
			sigmaRap2 += y2;
		//}
	}
	sigmaRap0 = sigmaRap0/(NumP-0.);
	sigmaRap1 = sigmaRap1/(NumP-0.);
	sigmaRap2 = sigmaRap2/(NumP-0.);
	cout<<"sigmaRap0: "<<sigmaRap0<<endl;
	cout<<"sigmaRap1: "<<sigmaRap1<<endl;
	cout<<"sigmaRap2: "<<sigmaRap2<<endl;

	// Psi 1S
	NumP =  graph_sigmaP_1S[0]->GetN();
	sigmaRap0 = 0; sigmaRap1 = 0;
	cout<<"NumP: "<<NumP<<endl;
	for(int i=0; i<NumP; i++){
		double x0,y0,x1,y1;
		graph_sigmaP_1S[0]->GetPoint(i,x0,y0);
		graph_sigmaP_1S[1]->GetPoint(i,x1,y1);

		if(i<11){
			cout<<"i: "<<i<<endl;
			sigmaRap0 += y0;
			sigmaRap1 += y1;
		}
	}
	sigmaRap0 = sigmaRap0/(NumP-1.);
	sigmaRap1 = sigmaRap1/(NumP-1.);
	cout<<"sigmaRap0: "<<sigmaRap0<<endl;
	cout<<"sigmaRap1: "<<sigmaRap1<<endl;

	graph_sigmaP_1S[0]->SetMarkerStyle(20);
	graph_sigmaP_1S[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_sigmaP_1S[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_sigmaP_1S[1]->SetMarkerStyle(24);
	graph_sigmaP_1S[1]->SetMarkerSize(1.2);
	graph_sigmaP_1S[1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_sigmaP_1S[1]->SetLineColor(onia::colour_rapForPTBins[2]);

	graph_sigmaP_2S[0]->SetMarkerStyle(21);
	graph_sigmaP_2S[0]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_sigmaP_2S[0]->SetLineColor(onia::colour_rapForPTBins[3]);
	graph_sigmaP_2S[1]->SetMarkerStyle(25);
	graph_sigmaP_2S[1]->SetMarkerSize(1.2);
	graph_sigmaP_2S[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_sigmaP_2S[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	graph_sigmaP_2S[2]->SetMarkerStyle(22);
	graph_sigmaP_2S[2]->SetMarkerSize(1.2);
	graph_sigmaP_2S[2]->SetMarkerColor(onia::colour_rapForPTBins[5]);
	graph_sigmaP_2S[2]->SetLineColor(onia::colour_rapForPTBins[5]);

	graph_sigmaP_2S[3]->SetMarkerStyle(33);
	graph_sigmaP_2S[3]->SetMarkerSize(1.5);
	graph_sigmaP_2S[3]->SetMarkerColor(kGreen);
	graph_sigmaP_2S[3]->SetLineColor(kGreen);

	//double blX = 0.6, blY = 0.63, trX = 0.83, trY = 0.84;
	double blX = 0.7, blY = 0.65, trX = 0.95, trY = 0.93;
	TLegend* legend_sigmaP_1S_2S=new TLegend(blX,blY,trX,trY);
	legend_sigmaP_1S_2S->SetFillColor(kWhite);
	legend_sigmaP_1S_2S->SetTextFont(42);
	legend_sigmaP_1S_2S->SetTextSize(0.035);
	legend_sigmaP_1S_2S->SetBorderSize(0.);
	legend_sigmaP_1S_2S->AddEntry(graph_sigmaP_1S[0],"J/#psi |y| < 0.6","lp");
	legend_sigmaP_1S_2S->AddEntry(graph_sigmaP_1S[1],"J/#psi 0.6 < |y| < 1.2","lp");
	legend_sigmaP_1S_2S->AddEntry(graph_sigmaP_2S[0],"#psi(2S) |y| < 0.6","lp");
	legend_sigmaP_1S_2S->AddEntry(graph_sigmaP_2S[1],"#psi(2S) 0.6 < |y| < 1.2","lp");
	legend_sigmaP_1S_2S->AddEntry(graph_sigmaP_2S[2],"#psi(2S) 1.2 < |y| < 1.5","lp");
	//legend_sigmaP_1S_2S->AddEntry(graph_sigmaP_2S[3],"#psi(2S) average","lp");

	graph_sigmaP_1S[0]->Draw("AP");
	graph_sigmaP_1S[1]->Draw("P");
	graph_sigmaP_2S[0]->Draw("P");
	graph_sigmaP_2S[1]->Draw("P");
	graph_sigmaP_2S[2]->Draw("P");
	//graph_sigmaP_2S[3]->Draw("P");
	legend_sigmaP_1S_2S->Draw();
	c1->SaveAs(Form("DataFiles/%s/RMS_1S_2S_12April2013.pdf",dataID));

	/////////////////////// both 1S and 2S && RMS * pT / M
	infile_1S = new TFile(
			Form("DataFiles/%s/Psi1S/Fit/parameter/evaluateCtau/PR_2.5sigma/graphFrac_2.5sigma.root",dataID),"R");
	infile_2S = new TFile(
			Form("DataFiles/%s/Psi2S/Fit/parameter/evaluateCtau/PR_2sigma/graphFrac_2.0sigma.root",dataID),"R");

	if(!infile_1S){cout<<"no input file"<<endl; return;}
	if(!infile_2S){cout<<"no input file"<<endl; return;}
	TGraphAsymmErrors *graph_sigmaP_L_1S[2], *graph_sigmaP_L_2S[3];
	graph_sigmaP_L_1S[0] = (TGraphAsymmErrors*)infile_1S->Get("graph_sigmaP_L_0");
	graph_sigmaP_L_1S[1] = (TGraphAsymmErrors*)infile_1S->Get("graph_sigmaP_L_1");
	graph_sigmaP_L_2S[0] = (TGraphAsymmErrors*)infile_2S->Get("graph_sigmaP_L_0");
	graph_sigmaP_L_2S[1] = (TGraphAsymmErrors*)infile_2S->Get("graph_sigmaP_L_1");
	graph_sigmaP_L_2S[2] = (TGraphAsymmErrors*)infile_2S->Get("graph_sigmaP_L_2");

	//remove first pT bin for 2S
	//graph_sigmaP_L_2S[0] -> RemovePoint(0);
	//graph_sigmaP_L_2S[1] -> RemovePoint(0);
	//graph_sigmaP_L_2S[2] -> RemovePoint(0);



	graph_sigmaP_L_1S[0]->SetMarkerStyle(20);
	graph_sigmaP_L_1S[0]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_sigmaP_L_1S[0]->SetLineColor(onia::colour_rapForPTBins[2]);
	graph_sigmaP_L_1S[1]->SetMarkerStyle(24);
	graph_sigmaP_L_1S[1]->SetMarkerSize(1.2);
	graph_sigmaP_L_1S[1]->SetMarkerColor(onia::colour_rapForPTBins[2]);
	graph_sigmaP_L_1S[1]->SetLineColor(onia::colour_rapForPTBins[2]);

	graph_sigmaP_L_2S[0]->SetMarkerStyle(21);
	graph_sigmaP_L_2S[0]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_sigmaP_L_2S[0]->SetLineColor(onia::colour_rapForPTBins[3]);
	graph_sigmaP_L_2S[1]->SetMarkerStyle(25);
	graph_sigmaP_L_2S[1]->SetMarkerSize(1.2);
	graph_sigmaP_L_2S[1]->SetMarkerColor(onia::colour_rapForPTBins[3]);
	graph_sigmaP_L_2S[1]->SetLineColor(onia::colour_rapForPTBins[3]);
	graph_sigmaP_L_2S[2]->SetMarkerStyle(22);
	graph_sigmaP_L_2S[2]->SetMarkerSize(1.2);
	graph_sigmaP_L_2S[2]->SetMarkerColor(onia::colour_rapForPTBins[5]);
	graph_sigmaP_L_2S[2]->SetLineColor(onia::colour_rapForPTBins[5]);

	//fit
  TF1 *lineFunc=new TF1("lineFunc","[0]+[1]*x",10.,65.);
	lineFunc->SetParameters(0.073,0.0027);
	lineFunc->SetParNames("a","b");
	lineFunc->SetLineColor(kBlack);
	lineFunc->SetLineStyle(kDashed);
	graph_sigmaP_L_1S[1]->Fit("lineFunc","R+");
	double *par, *parErr;
	par    = lineFunc->GetParameters();
	parErr = lineFunc->GetParErrors();
	double chi2 = lineFunc->GetChisquare();
	int    ndf  = lineFunc->GetNDF();
	cout<<"a:    "<<par[0]<<" +/- "<<parErr[0]<<endl;
	cout<<"b:    "<<par[1]<<" +/- "<<parErr[1]<<endl;
	cout<<"chi2: "<<chi2<<endl;
	cout<<"ndf:  "<<ndf<<endl;
	//

	//blX = 0.18; blY = 0.63; trX = 0.51; trY = 0.84;
	blX = 0.12; blY = 0.65; trX = 0.4; trY = 0.93;
	TLegend* legend_sigmaP_L_1S_2S=new TLegend(blX,blY,trX,trY);
	legend_sigmaP_L_1S_2S->SetFillColor(kWhite);
	legend_sigmaP_L_1S_2S->SetTextFont(42);
	legend_sigmaP_L_1S_2S->SetTextSize(0.035);
	legend_sigmaP_L_1S_2S->SetBorderSize(0.);
	legend_sigmaP_L_1S_2S->AddEntry(graph_sigmaP_L_1S[0],"J/#psi |y| < 0.6","lp");
	legend_sigmaP_L_1S_2S->AddEntry(graph_sigmaP_L_1S[1],"J/#psi 0.6 < |y| < 1.2","lp");
	legend_sigmaP_L_1S_2S->AddEntry(graph_sigmaP_L_2S[0],"#psi(2S) |y| < 0.6","lp");
	legend_sigmaP_L_1S_2S->AddEntry(graph_sigmaP_L_2S[1],"#psi(2S) 0.6 < |y| < 1.2","lp");
	legend_sigmaP_L_1S_2S->AddEntry(graph_sigmaP_L_2S[2],"#psi(2S) 1.2 < |y| < 1.5","lp");
	legend_sigmaP_L_1S_2S->AddEntry(lineFunc,"Fit","lp");

	graph_sigmaP_L_1S[0]->GetYaxis()->SetTitle("L_{xy} r.m.s. (mm)");
	graph_sigmaP_L_1S[0]->Draw("AP");
	graph_sigmaP_L_1S[1]->Draw("P");
	graph_sigmaP_L_2S[0]->Draw("P");
	graph_sigmaP_L_2S[1]->Draw("P");
	graph_sigmaP_L_2S[2]->Draw("P");
	legend_sigmaP_L_1S_2S->Draw();

 

	double left=0.55, top=0.85, textSize=0.035;
	TLatex *latex=new TLatex();
	latex->SetTextFont(42);
	latex->SetNDC(kTRUE);
	latex->SetTextSize(textSize);
	double step=textSize*1.3;
	latex->DrawLatex(left,top,Form("#chi^{2} / ndf               %.4f / %d", chi2, ndf));
	top-=step;
	latex->DrawLatex(left,top,Form("a                %.3f #pm   %.4f (mm)", par[0], parErr[0]));
	top-=step;
	latex->DrawLatex(left,top,Form("b              %.4f #pm %.5f (mm/GeV)", par[1], parErr[1]));

	c1->SaveAs(Form("DataFiles/%s/RMS_L_1S_2S_12April2013.pdf",dataID));

	return ;
}
