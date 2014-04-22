#include "../interface/rootIncludes.inc"
#include "../interface/RooUtils.h"

const int rapBins_1s = 2;
const int rapBins_2s = 3;
const int ptBins_1s = 12;
const int ptBins_2s = 5;
double mean_1s[rapBins_1s][ptBins_1s], mean_2s[rapBins_2s][ptBins_2s];
double meanErr_1s[rapBins_1s][ptBins_1s], meanErr_2s[rapBins_2s][ptBins_2s];
double sigma_1s[rapBins_1s][ptBins_1s], sigma_2s[rapBins_2s][ptBins_2s];
double sigmaErr_1s[rapBins_1s][ptBins_1s], sigmaErr_2s[rapBins_2s][ptBins_2s];
char filePath[500];
char dataPath[500];

void fitPullMassLife(char *patch="SR"){

  sprintf(dataPath,"SetOfCuts11_ctauScen0_FracLSB-1_newMLfit_4Mar2013");

	double par[3];
	TF1* gauss = new TF1("gauss", "gaus", -5., 5.);

	TCanvas *c1=new TCanvas("c1","");
	c1->SetGridx();
	c1->SetGridy();
	gStyle->SetTitleFillColor(10);
	gStyle->SetOptFit(1);            
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
	gStyle->SetStatW(0.16);
	gStyle->SetStatH(0.15);
	gStyle->SetTitleX(0.15);
	gStyle->SetTitleY(0.96);

	// 1S 
	sprintf(filePath,"DataFiles/%s/Psi1S/Fit/pull",dataPath);
	gSystem->mkdir(Form("%s/Figures",filePath));
	for(int rap = 1; rap <= rapBins_1s; rap++){
		for(int pt = 1; pt <= ptBins_1s; pt++){
			cout<<"rap "<<rap<<" pt "<<pt<<endl;

			TFile *pullFile = new TFile(Form("%s/pull_rap%d_pt%d_%s.root",filePath,rap,pt,patch), "R");
			if(!pullFile){cout<<"no input file"<<endl; continue;}

			TH1F* pullHist = (TH1F*)pullFile->Get("pull");
			pullHist->SetXTitle("pull");
			pullHist->SetTitle("");

			double norm =pullHist->GetMaximum(), mean=0.,sigma=1.;

			for(int fit=0; fit<20; fit++){
				gauss -> SetParameters(norm,mean,sigma);
				pullHist->Fit(gauss,"R+");
				gauss->GetParameters(par);
				double ndof = pullHist->GetRMS()/pullHist->GetBinWidth(1) - 3;
				cout<<"chi2/ndf = "<<gauss->GetChisquare()/ndof<<endl;
				norm=par[0]*2.; mean=par[1]; sigma=par[2];
				if(gauss->GetChisquare()/ndof<5.) break;
			}
			mean_1s [rap-1][pt-1] = mean ;
			sigma_1s[rap-1][pt-1] = sigma ;

			meanErr_1s [rap-1][pt-1] = gauss->GetParError(1) ;
			sigmaErr_1s[rap-1][pt-1] = gauss->GetParError(2) ;

			pullHist->Draw();
			gauss->Draw("same");

			c1->SaveAs(Form("%s/Figures/pullFit_rap%d_pt%d_%s.pdf",filePath,rap,pt,patch));

		} //pt
	} //rap
	/////

	// 2S 
	sprintf(filePath,"DataFiles/%s/Psi2S/Fit/pull",dataPath);
	gSystem->mkdir(Form("%s/Figures",filePath));
	for(int rap = 1; rap <= rapBins_2s; rap++){
		for(int pt = 1; pt <= ptBins_2s; pt++){
			cout<<"rap "<<rap<<" pt "<<pt<<endl;

			TFile *pullFile = new TFile(Form("%s/pull_rap%d_pt%d_%s.root",filePath,rap,pt,patch), "R");
			if(!pullFile){cout<<"no input file"<<endl; continue;}

			TH1F* pullHist = (TH1F*)pullFile->Get("pull");
			pullHist->SetXTitle("pull");
			pullHist->SetTitle("");

			double norm =pullHist->GetMaximum(), mean=0.,sigma=1.;

			for(int fit=0; fit<20; fit++){
				gauss -> SetParameters(norm,mean,sigma);
				pullHist->Fit(gauss,"R+");
				gauss->GetParameters(par);
				double ndof = pullHist->GetRMS()/pullHist->GetBinWidth(1) - 3;
				cout<<"chi2/ndf = "<<gauss->GetChisquare()/ndof<<endl;
				norm=par[0]*2.; mean=par[1]; sigma=par[2];
				if(gauss->GetChisquare()/ndof<5.) break;
			}

			mean_2s [rap-1][pt-1] = mean ;
			sigma_2s[rap-1][pt-1] = sigma ;

			meanErr_2s [rap-1][pt-1] = gauss->GetParError(1) ;
			sigmaErr_2s[rap-1][pt-1] = gauss->GetParError(2) ;

			pullHist->Draw();
			gauss->Draw("same");

			c1->SaveAs(Form("%s/Figures/pullFit_rap%d_pt%d_%s.pdf",filePath,rap,pt,patch));

		} //pt
	} //rap
	/////


	////////////////////////////////////////////////////////////////////////

	double pT1s[2][12]={
		{11.0099, 12.9394, 14.9263, 16.9274, 18.9324, 20.9341, 23.3586, 27.1481, 32.1847, 37.2121, 44.0308, 56.8306},
		{10.9643, 12.9246, 14.9218, 16.9219, 18.9241, 20.9311, 23.3568, 27.1473, 32.1828, 37.2213, 44.0189, 56.8043}};
	double pT1sErr[2][12]={
		{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
		{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}};
	double pT2s[3][5]= {
		{11.8943, 15.7294, 19.7371, 25.1272, 36.1102},
		{11.7868, 15.7024, 19.7270, 25.1119, 36.1893},
		{11.7788, 15.7210, 19.7275, 25.1509, 36.1740}};
	double pT2sErr[3][5]= {
		{0., 0., 0., 0., 0.},
		{0., 0., 0., 0., 0.},
		{0., 0., 0., 0., 0.}};

	gStyle->SetPadBottomMargin(0.11);
	gStyle->SetPadLeftMargin(0.08); 
	gStyle->SetPadRightMargin(0.02);
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetOptFit(1);
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

	TCanvas *c3 = new TCanvas("c3","c3",1200,800);
	gPad->SetFillColor(kWhite);

	TGraphErrors *pull_mean_1s[rapBins_1s];
	TGraphErrors *pull_mean_2s[rapBins_1s];
	TGraphErrors *pull_sigma_1s[rapBins_1s];
	TGraphErrors *pull_sigma_2s[rapBins_1s];

	for(int rap=0; rap<rapBins_1s; rap++){
		pull_mean_1s[rap] = new TGraphErrors(ptBins_1s, pT1s[rap], mean_1s[rap], pT1sErr[rap], meanErr_1s[rap] );
		pull_mean_1s[rap]->SetTitle("");
		pull_mean_1s[rap]->GetXaxis()->SetTitle("p_{T} (GeV)");
		pull_mean_1s[rap]->GetYaxis()->SetTitle("mean of pull");
		pull_mean_1s[rap]->GetYaxis()->SetRangeUser(-1., 1.);
		pull_mean_1s[rap]->GetXaxis()->SetLimits(6., 72.);

		pull_sigma_1s[rap] = new TGraphErrors(ptBins_1s, pT1s[rap], sigma_1s[rap], pT1sErr[rap], sigmaErr_1s[rap] );
		pull_sigma_1s[rap]->SetTitle("");
		pull_sigma_1s[rap]->GetXaxis()->SetTitle("p_{T} (GeV)");
		pull_sigma_1s[rap]->GetYaxis()->SetTitle("sigma of pull");
		pull_sigma_1s[rap]->GetYaxis()->SetRangeUser(0., 2.);
		pull_sigma_1s[rap]->GetXaxis()->SetLimits(6., 72.);
	}

	for(int rap=0; rap<rapBins_2s; rap++){
		pull_mean_2s[rap] = new TGraphErrors(ptBins_2s, pT2s[rap], mean_2s[rap], pT2sErr[rap], meanErr_2s[rap] );
		pull_mean_2s[rap]->SetTitle("");
		pull_mean_2s[rap]->GetXaxis()->SetTitle("p_{T} (GeV)");
		pull_mean_2s[rap]->GetYaxis()->SetTitle("mean of pull");
		pull_mean_2s[rap]->GetYaxis()->SetRangeUser(-1., 1.);
		pull_mean_2s[rap]->GetXaxis()->SetLimits(6., 72.);

		pull_sigma_2s[rap] = new TGraphErrors(ptBins_2s, pT2s[rap], sigma_2s[rap], pT2sErr[rap], sigmaErr_2s[rap] );
		pull_sigma_2s[rap]->SetTitle("");
		pull_sigma_2s[rap]->GetXaxis()->SetTitle("p_{T} (GeV)");
		pull_sigma_2s[rap]->GetYaxis()->SetTitle("sigma of pull");
		pull_sigma_2s[rap]->GetYaxis()->SetRangeUser(0., 2.);
		pull_sigma_2s[rap]->GetXaxis()->SetLimits(6., 72.);
	}


	pull_mean_1s[0]->SetMarkerStyle(20);
	pull_mean_1s[0]->SetMarkerColor(4);
	pull_mean_1s[0]->SetLineColor(4);
	pull_mean_1s[1]->SetMarkerStyle(24);
	pull_mean_1s[1]->SetMarkerSize(1.2);
	pull_mean_1s[1]->SetMarkerColor(4);
	pull_mean_1s[1]->SetLineColor(4);
	pull_mean_2s[0]->SetMarkerStyle(21);
	pull_mean_2s[0]->SetMarkerColor(2);
	pull_mean_2s[0]->SetLineColor(2);
	pull_mean_2s[1]->SetMarkerStyle(25);
	pull_mean_2s[1]->SetMarkerSize(1.2);
	pull_mean_2s[1]->SetMarkerColor(2);
	pull_mean_2s[1]->SetLineColor(2);
	pull_mean_2s[2]->SetMarkerStyle(22);
	pull_mean_2s[2]->SetMarkerSize(1.2);
	pull_mean_2s[2]->SetMarkerColor(kMagenta+1);
	pull_mean_2s[2]->SetLineColor(kMagenta+1);

	pull_sigma_1s[0]->SetMarkerStyle(20);
	pull_sigma_1s[0]->SetMarkerColor(4);
	pull_sigma_1s[0]->SetLineColor(4);
	pull_sigma_1s[1]->SetMarkerStyle(24);
	pull_sigma_1s[1]->SetMarkerSize(1.2);
	pull_sigma_1s[1]->SetMarkerColor(4);
	pull_sigma_1s[1]->SetLineColor(4);
	pull_sigma_2s[0]->SetMarkerStyle(21);
	pull_sigma_2s[0]->SetMarkerColor(2);
	pull_sigma_2s[0]->SetLineColor(2);
	pull_sigma_2s[1]->SetMarkerStyle(25);
	pull_sigma_2s[1]->SetMarkerSize(1.2);
	pull_sigma_2s[1]->SetMarkerColor(2);
	pull_sigma_2s[1]->SetLineColor(2);
	pull_sigma_2s[2]->SetMarkerStyle(22);
	pull_sigma_2s[2]->SetMarkerSize(1.2);
	pull_sigma_2s[2]->SetMarkerColor(kMagenta+1);
	pull_sigma_2s[2]->SetLineColor(kMagenta+1);

	double blX = 0.7, blY = 0.7, trX = 0.95, trY = 0.93;
	TLegend* legend=new TLegend(blX,blY,trX,trY);
	legend->SetFillColor(kWhite);
	legend->SetTextFont(42);
	legend->SetTextSize(0.035);
	legend->SetBorderSize(0.);
	legend->AddEntry(pull_mean_1s[0],"J/#psi |y| < 0.6","lp");
	legend->AddEntry(pull_mean_1s[1],"J/#psi 0.6 < |y| < 1.2","lp");
	legend->AddEntry(pull_mean_2s[0],"#psi(2S) |y| < 0.6","lp");
	legend->AddEntry(pull_mean_2s[1],"#psi(2S) 0.6 < |y| < 1.2","lp");
	legend->AddEntry(pull_mean_2s[2],"#psi(2S) 1.2 < |y| < 1.5","lp");

	TLine *line_one = new TLine(6., 1., 72., 1.);
	line_one->SetLineWidth(1.5); line_one->SetLineColor(kBlack); line_one->SetLineStyle(7);

	TLine *line_zero = new TLine(6., 0., 72., 0.);
	line_zero->SetLineWidth(1.5); line_zero->SetLineColor(kBlack); line_zero->SetLineStyle(7);

	pull_mean_1s[0]->Draw("AP");
	pull_mean_1s[1]->Draw("P");
	pull_mean_2s[0]->Draw("P");
	pull_mean_2s[1]->Draw("P");
	pull_mean_2s[2]->Draw("P");
	line_zero->Draw("same");
	legend->Draw("same");
	c3->SaveAs(Form("DataFiles/%s/pull_mean_%s.pdf",dataPath,patch));

	pull_sigma_1s[0]->Draw("AP");
	pull_sigma_1s[1]->Draw("P");
	pull_sigma_2s[0]->Draw("P");
	pull_sigma_2s[1]->Draw("P");
	pull_sigma_2s[2]->Draw("P");
	line_one->Draw("same");
	legend->Draw("same");
	c3->SaveAs(Form("DataFiles/%s/pull_sigma_%s.pdf",dataPath,patch));

}
