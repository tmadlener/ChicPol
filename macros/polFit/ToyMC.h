namespace ToyMC{

double ScenarioSig [3][9]={{-1,-0.5,0,0.5,1,0,0,0.5,1},{0,0,0,0,0,0.5,-0.5,-0.75,-1},{0,0,0,0,0,0,0,0,0}};//lamth_Signal,lamph_Signal,lamthph_Signal

double ScenarioBkg [3][9]={{-1,-0.5,0,0.5,1,0,0,2,4},{0,0,0,0,0,0.5,0.8,0.4,-0.4},{0,0,0,0,0,0,0,0,0}};////lamth_Bkg,lamph_Bkg,lamthph_Bkg


/*int MarkerStyle[3] = {0, 20, 21};
int MarkerStyle2[3] = {0, 24, 25};
int MarkerStyle3[3] = {0, 26, 32};

int MarkerColor[3] = {0, 601, 632};
int MarkerColor2[3] = {0, 632, 632};
int MarkerColor3[3] = {0, 418, 632};
*/
//int MarkerStyle[6][4]={{0,0,0,0},{0,33,27,34},{0,20,24,29},{0,21,25,22},
//	{0,33,27,34},{0,20,24,29}}; // for each state, rapBin (1= closed, 2=open)
//int MarkerColor[6] = {0,1,1,1,1,1};//{0,600,632,418}; // for each frame
//double MarkerSize[6][4]={{0,0,0},{0,2.75,2.75,2.75},{0,1.65,1.65,1.65},{0,1.65,1.65,1.65},
//	{0,2.75,2.75,2.75},{0,1.65,1.65,1.65}};// for each state, rapBin

int MarkerStyle[8]={0,33,27,34,20,20,20,20}; // for each state, rapBin (1= closed, 2=open)
int MarkerColor[4] = {0,600,632,418}; // for each frame
double MarkerSize[8]={1.,1.,1.,1.,1.,1.,1.,1.};// for each state, rapBin



//Chic dummy values
const int nPtBins=5;
const int nRapBins=1;

//chic1
double fracBackground[nRapBins][nPtBins]={{0.3,0.3,0.3,0.3,0.3}};//0.001
double fracBackgrounderr[nRapBins][nPtBins]={{0.,0.,0.,0.,0.}};
double ptCentre[nRapBins][nPtBins]={{12.5, 17.5, 22.5, 27.5, 40.}};
int numEvents[nRapBins][nPtBins]={{20000, 20000, 10000, 5000, 4000}};
double meanRap[nRapBins][nPtBins]={{0.6, 0.6, 0.6, 0.6, 0.6}};

//chic2
//double fracBackground[nRapBins][nPtBins]={{0.4,0.4,0.4,0.4,0.4}};
//double fracBackgrounderr[nRapBins][nPtBins]={{0.,0.,0.,0.,0.}};
//double ptCentre[nRapBins][nPtBins]={{12.5, 17.5, 22.5, 27.5, 40.}};
//int numEvents[nRapBins][nPtBins]={{8000, 8000, 3500, 1500, 1100}};
//double meanRap[nRapBins][nPtBins]={{0.6, 0.6, 0.6, 0.6, 0.6}};


int binCosth[nRapBins][nPtBins]={{16, 16, 16, 16, 16}};
int binPhi  [nRapBins][nPtBins]={{16, 16, 16, 16, 16}};



double fracSignal[nRapBins][nPtBins]={{0, 0, 0, 0, 0}};

const int nEffs=3;
const int FidCuts=3;

}
