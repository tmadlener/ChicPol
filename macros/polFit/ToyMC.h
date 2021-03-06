namespace ToyMC{

double ScenarioSig [3][9]={{-1,-0.5,0,0.5,1,0,0,0.5,1},{0,0,0,0,0,0.5,-0.5,-0.75,-1},{0,0,0,0,0,0,0,0,0}};//lamth_Signal,lamph_Signal,lamthph_Signal

double ScenarioBkg [3][9]={{-1,-0.5,0,-0.25,1,0,0,2,4},{0,0,0,0.15,0,0.5,0.8,0.4,-0.4},{0,0,0,0,0,0,0,0,0}};////lamth_Bkg,lamph_Bkg,lamthph_Bkg

int MarkerStyle[6] = {25,26,27,28,30,32}; //Ilse
//int MarkerStyle[6][4]={{0,0,0,0},{0,33,27,34},{0,20,24,29},{0,21,25,22},
//	{0,33,27,34},{0,20,24,29}}; // for each state, rapBin (1= closed, 2=open)
int MarkerColor[6] = {0,1,2,3,1,1};//{0,600,632,418}; // for each frame
double MarkerSize[6] = {0,1.65,1.65,1.65,1.65,1.65}; //Ilse
//double MarkerSize[6][4]={{0,0,0},{0,2.75,2.75,2.75},{0,1.65,1.65,1.65},{0,1.65,1.65,1.65},
//	{0,2.75,2.75,2.75},{0,1.65,1.65,1.65}};// for each state, rapBin

const int nPtBins=5;
const int nRapBins=1;

// costh and phi bins from data BG histogram
int binCosth[nRapBins][nPtBins]={64, 48, 48, 48, 32};
int binPhi[nRapBins][nPtBins]={16, 16, 16, 16, 16};
//int binCosth[nRapBins][nPtBins]={{48, 34, 34, 32, 32, 32, 32, 32, 32, 32, 32, 16},{34, 29, 32, 32, 32, 32, 32, 32, 32, 32, 32, 16}};
//int binPhi  [nRapBins][nPtBins]={{16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16},{16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16}};

// ToyBackground
double fracBackground[nRapBins][nPtBins]={{0.01,0.01,0.01,0.01,0.01}};
double fracBackgrounderr[nRapBins][nPtBins]={{0.01,0.01,0.01,0.01,0.01}}; //Ilse
//// from real data
//double fracBackground[nRapBins][nPtBins]={{0.015, 0.015, 0.016, 0.018, 0.019, 0.020, 0.022, 0.025, 0.029, 0.033, 0.035, 0.046},{0.021, 0.021, 0.024, 0.025, 0.027, 0.029, 0.032, 0.034, 0.036, 0.043, 0.042, 0.060}};

double ptCentre[nRapBins][nPtBins]={12.5,17.5,22.5,27.,35.}; //random for now
int numEvents[nRapBins][nPtBins]={30000,15000,15000,8000,5000}; // random for now
double meanRap[nRapBins][nPtBins]={0.6,0.6,0.6,0.6,0.6}; //random for now
double fracSignal[nRapBins][nPtBins]={0,0,0,0,0};

//Psi(1):
//double ptCentre[nRapBins][nPtBins]={{11.0099,  12.9394, 14.9263, 16.9274, 18.9324, 20.9341, 23.3586, 27.1481, 32.1847, 37.2121, 44.0308, 56.8306},{10.9643, 12.9246, 14.9218, 16.9219, 18.9241, 20.9311, 23.3568, 27.1473, 32.1828, 37.2213, 44.0189, 56.8043}};
//int numEvents[nRapBins][nPtBins]={{436221, 327834, 200584, 118810, 72848, 45256, 39082, 30304, 12653, 5874, 4420, 1794},{555974, 357350, 205358, 117992, 69825, 42832, 37194, 28632, 11699, 5406, 3997, 1574}}; //Psi(1S)
//double meanRap[nRapBins][nPtBins]={{0.308202, 0.309691, 0.309141, 0.309356, 0.304804, 0.300372, 0.295868, 0.295074, 0.293911, 0.298604, 0, 0},{0.88879, 0.881283, 0.875584, 0.866698, 0.862725, 0.859332, 0.859927, 0.867519, 0.873316, 0.875993, 0, 0}};
//Total Number of Signal Events in safe region: 251987
//double fracSignal[nRapBins][nPtBins]={{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{4.996e-15, 4.29656e-14, 2.75335e-14, 2.45359e-14, 2.55351e-14, 4.62963e-14, 7.33857e-14, 2.13385e-13, 6.54921e-13, 6.82188e-12, 0, 0}};

const int nEffs=3;
const int FidCuts=3;

}
