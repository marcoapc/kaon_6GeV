#include <sys/stat.h>
#include "/media/marco/Data/Marco/Documents/Marco_CUA/codes/MLibraries/MallIncludes.h"
#include "/media/marco/Data/Marco/Documents/Marco_CUA/codes/MLibraries/MConfigGraphs.h"

struct DefRanges {
	Double_t Q2[2];
	Double_t W[2];
	Double_t WQ2cutW[4];// = {2.322, 2.422, 2.33, 2.267}; //x
	Double_t WQ2cutQ2[4];// = {2.02, 1.75, 2.218, 2.362}; //y
};

Double_t ReadScaler(TString filenameInit, Int_t runN);

int MDraw2DMM(TString name, TTree *EX, TCut mcut, Int_t nboxes, Double_t *xc, Double_t *xw, Double_t *yc, Double_t *yw, Double_t *range, TH2F *a, Int_t nbin=100, Double_t maxZ=0.0);

Double_t Fit2LinesRanges(Double_t *x, Double_t *par) {
	// x<=[8]  ->   [0]*([4]+[5]*x) + [1]*exp(-((x-[3])/[2])^2)
	// x>[8]   ->   [0]*([6]+[7]*x)
	Double_t value;
	if(x[0] <= par[8]) value = (par[0]*(par[4] + par[5]*x[0]));
	else value = (par[0]*(par[6] + par[7]*x[0]));
	return (value + par[1]*exp(-pow((x[0]-par[3])/par[2],2.0)));
}

Double_t FitPol3(Double_t *x, Double_t *par) {
	return (par[0]*(par[4]+par[5]*x[0]+par[6]*x[0]*x[0]+par[7]*x[0]*x[0]*x[0])+par[1]*exp(-pow((x[0]-par[3])/par[2],2)));
}

Double_t LorentzPeakPol3(Double_t *x, Double_t *par) {
	return (par[0]*(par[4]+par[5]*x[0]+par[6]*x[0]*x[0]+par[7]*x[0]*x[0]*x[0])+par[1]/(par[2]+pow(x[0]-par[3],2)/par[2]));
}

// use t_init and t_final only when autoRun=1
int MlookSIMC(Int_t runN=2, Int_t runNfinal=0, Int_t dummyRun=0, Int_t dummyRunfinal=0, Bool_t autoRun=0, Double_t t_init=0, Double_t t_final=0) {

	// FLAGS
	Int_t type = 2; //1-elastic e(p,e') on HMS only | 2-p(e,K+)L/S
	Int_t savePDF = 0;
	Int_t fast = 1; // only some plots
	Int_t writeLog = 0;
	Int_t SubtractDummy = 0;
	TString SIMCfolder = "Rosenbluth2015_006";
	Double_t t_range[2] = {0.0,1.0};
	//Double_t t_range[2] = {0.30,0.41};
	//Double_t t_range[2] = {0.41,0.49};
	//Double_t t_range[2] = {0.49,0.60};
	Bool_t scaleSIMC = 0;
	Double_t sSIMC = 1.0;//0.119072;//1308.74; 
	Int_t NotAnalyze[] = {47055, 47060, 47061, 47063, 47064, 47065, 47066, 47076, 47148, 47368, 47397, 47403, 47420, 47425, 47592, 47595, 47604, 47694, 47697, 47729, 47736}; // removing runs
	// For background analysis
	Bool_t analyzeBack = 1; // To run MM background fitting
	Int_t saveRes = 0; // save results from Phi cut on <filePhiRes>
	TString filePhiRes = "fPhi.csv"; // details of the analysis
	// phibin
	const int nphibins = 10;
	Double_t phiBinsLimits[nphibins+1] = {0.0, 36.0, 72.0, 108.0, 144.0, 180.0, 216.0, 252.0, 288.0, 324.0, 360.0};

	// IF RUNNING IN BATCH -> use script autoRunAll.C
	if(autoRun) {
		if(t_init >= t_final) { t_range[0]=0.0; t_range[1]=1.0; }
		else { t_range[0]=t_init; t_range[1]=t_final; }
	}

	// Constants and Variables
	Double_t eNorm = 1.0;
	if(eNorm!=1.0) {
		for(i=0; i<20; i++) cout << "*";
		cout << endl << "** ATTENTION!  the final yield by: " << eNorm << endl;
		for(i=0; i<20; i++) cout << "*";
	}
	TString SIMtype="l_kaon_h";
	Double_t LambdaMM = 1.115683;
	Int_t i,j,k;
	Long_t *id,*size,*flags,*mt;
	Double_t temp;
	TString cname;
	TString epstype;
	TString outfolderEnd = Form("t_gt_%.3f_lt_%.3f",t_range[0],t_range[1]);
	Int_t flagRanges;

	// Making special cases
	if(runN==0) { // elastic run to test
		runN = 47345; // elastic run to test
		type = 1;
	}
	if(runN==1) { //1st kinematics Kaon on FPI2
		runN = 47056;//47688;
		runNfinal = 47087;
		type = 2;
		epstype = "HighQ2_loweps";
		flagRanges=1;
	}
	if(runN==2) { //2nd kinematics Kaon on FPI2
		runN = 47358;//47688;
		runNfinal = 47371;
		type = 2;
		epstype = "HighQ2_higheps";
		flagRanges=1;
	}
	if(runN==3) { //2nd kinematics Kaon on FPI2
		runN = 47419;//47688;
		runNfinal = 47428;
		type = 2;
		epstype = "HighQ2_higheps";
		flagRanges=1;
	}
	if(runN==4) { //2nd kinematics Kaon on FPI2
		runN = 47393;//47688;
		runNfinal = 47405;
		type = 2;
		epstype = "HighQ2_higheps";
		flagRanges=1;
	}
	if(runN==5) { //1st kinematics Kaon on FPI2
		runN = 47120;//47688;
		runNfinal = 47154;
		type = 2;
		epstype = "HighQ2_loweps";
		flagRanges=1;
	}
	if(runN==6) { //3rd kinematics Kaon on FPI2 (low Q2, low eps)
		runN = 47587;
		runNfinal = 47597;
		type = 2;
		epstype = "LowQ2_loweps";
		flagRanges=2;
	}
	if(runN==7) { //3rd kinematics Kaon on FPI2 (low Q2, low eps)
		runN = 47600;
		runNfinal = 47613;
		type = 2;
		epstype = "LowQ2_loweps";
		flagRanges=2;
	}
	if(runN==8) { //4th kinematics Kaon on FPI2 (low Q2, high eps)
		runN = 47687;
		runNfinal = 47692;
		type = 2;
		epstype = "LowQ2_higheps";
		flagRanges=2;
	}
	if(runN==9) { //4th kinematics Kaon on FPI2 (low Q2, high eps)
		runN = 47695;
		runNfinal = 47701;
		type = 2;
		epstype = "LowQ2_higheps";
		flagRanges=2;
	}
	if(runN==10) { //4th kinematics Kaon on FPI2 (low Q2, high eps)
		runN = 47725;
		runNfinal = 47736;
		type = 2;
		epstype = "LowQ2_higheps";
		flagRanges=2;
	}
	if(runN==-1) { // load all Q2=2.45 and eps=0.54 together (no weight, be careful!) for Q2vsW region analysis
		runN=47358;
		runNfinal=47428;
		type=2;
		Int_t NotAnalyze[] = {47076, 47066, 47368, 47372, 47373, 47374, 47375, 47376, 47377, 47378, 47379, 47380, 47381, 47382, 47383, 47384, 47385, 47386, 47387, 47388, 47389, 47390, 47391, 47392, 47397, 47403, 47406, 47407, 47408, 47409, 47410, 47411, 47412, 47413, 47414, 47415, 47416, 47417, 47418, 47420, 47425}; // removing runs (comment the previous declared NotAnalyze variable)
		flagRanges=1;
	}
	if(runN==-2) { // load all Q2=2.45 and eps=0.27 together (no weight, be careful!) for Q2vsW region analysis
		runN=47056;
		runNfinal=47154;
		type=2;
		Int_t NotAnalyze[] = {47368, 47372, 47373, 47374, 47375, 47376, 47377, 47378, 47379, 47380, 47381, 47382, 47383, 47384, 47385, 47386, 47387, 47388, 47389, 47390, 47391, 47392, 47397, 47403, 47406, 47407, 47408, 47409, 47410, 47411, 47412, 47413, 47414, 47415, 47416, 47417, 47418, 47420, 47425, 47066, 47076, 47088, 47089, 47090, 47091, 47092, 47093, 47094, 47095, 47096, 47097, 47098, 47099, 47100, 47101, 47102, 47103, 47104, 47105, 47106, 47107, 47108, 47109, 47110, 47111, 47112, 47113, 47114, 47115, 47116, 47117, 47118}; // removing runs (comment the previous declared NotAnalyze variable)
		flagRanges=1;
	}

	// RANGES
	DefRanges ranges;
	if(flagRanges==1) { // FPI-2, high Q2
		// Q2
		ranges.Q2[0] = 1.45; // plot only
		ranges.Q2[1] = 2.8; // plot only
		// W
		ranges.W[0] = 2.15; // plot only
		ranges.W[1] = 2.6; // plot only
		// Defining vertices of diamont-cut W vs Q2
		// Q2=2.45, eps=0.54 - full high eps
		//ranges.WQ2cutW[] = {2.34, 2.57, 2.425, 2.21}; //x -- REAL CUT
		//ranges.WQ2cutQ2[] = {1.9, 1.525, 2.2575, 2.75}; //y -- REAL CUT
		// Q2=2.45, eps=0.27
		ranges.WQ2cutW[0]=2.322; ranges.WQ2cutW[1]=2.422; ranges.WQ2cutW[2]=2.33; ranges.WQ2cutW[3]=2.267; //x -- REAL CUT
		ranges.WQ2cutQ2[0]=2.02; ranges.WQ2cutQ2[1]=1.75; ranges.WQ2cutQ2[2]=2.218; ranges.WQ2cutQ2[3]=2.362; //y -- REAL CUT
	}
	else if(flagRanges==2) { // FPI-2, low Q2
		// Q2
		ranges.Q2[0] = 0.7; // plot only
		ranges.Q2[1] = 2.0; // plot only
		// W
		ranges.W[0] = 2.15; // plot only
		ranges.W[1] = 2.6; // plot only
		// Defining vertices of diamont-cut W vs Q2
		// Q2=2.45, eps=0.54 - full high eps
		//ranges.WQ2cutW[] = {2.34, 2.57, 2.425, 2.21}; //x -- REAL CUT
		//ranges.WQ2cutQ2[] = {1.9, 1.525, 2.2575, 2.75}; //y -- REAL CUT
		// Q2=2.45, eps=0.27
		ranges.WQ2cutW[0]=2.305; ranges.WQ2cutW[1]=2.37; ranges.WQ2cutW[2]=2.33; ranges.WQ2cutW[3]=2.265; //x -- REAL CUT
		ranges.WQ2cutQ2[0]=1.275; ranges.WQ2cutQ2[1]=1.15; ranges.WQ2cutQ2[2]=1.37; ranges.WQ2cutQ2[3]=1.59; //y -- REAL CUT
	}

	// Defining prefix_sufix for opening files
	TString EXtype;
	if(type==1) { //1-elastic
		SIMtype="eep_h";
		EXtype="hms";
	}
	else if(type==2) { //kaon
		EXtype="coin";
	}

	// Creating output folder, if it doesn't exist
	TString outFolder = "CrossSection/"+SIMCfolder+"/"+epstype+"/"+outfolderEnd;
	// MainFolder
	if(gSystem->GetPathInfo(outFolder,id,size,flags,mt) != 0) {
		if(gSystem->mkdir(outFolder.Data(),1)!=0) {
			cout << "ERROR: Could not access/create the folder: " << outFolder << endl;
			return 1;
		}
	}
	// Figures folder
	if(gSystem->GetPathInfo(outFolder+"/figures/"+Form("%d",runN),id,size,flags,mt) != 0) {
		if(gSystem->mkdir(outFolder+"/figures/"+Form("%d",runN),1)!=0) {
			cout << "ERROR: Could not access/create the folder: " << outFolder << "/figures/" << Form("%d",runN) << endl;
			return 1;
		}
	}
	// Figures of background fit folder
	if(gSystem->GetPathInfo(outFolder+"/figures/"+Form("%d",runN)+"/backfits",id,size,flags,mt) != 0) {
		if(gSystem->mkdir(outFolder+"/figures/"+Form("%d",runN)+"/backfits",1)!=0) {
			cout << "ERROR: Could not access/create the folder: " << outFolder << "/figures/" << Form("%d",runN) << "/backfits/" << endl;
			return 1;
		}
	}

	// Output logfile
	TString outLogName = "logAnalysis.csv";
	ofstream logFile;
	logFile.open(outLogName.Data(), std::ofstream::out | std::ofstream::app);
	Int_t logok = logFile.is_open();
	TDatime da;
	da.GetTime();
	if(logok && writeLog) {
		cout << "<Marco> Writing analysis log in file: " << outLogName << endl;
		logFile << endl << "->Run=" << runN << ",Now=" << da.AsSQLString() << ",";
		}
	else logok=0;
	if(!logok) {
		cout << "<Marco> Warning: analysis log NOT being saved." << endl;
		logok=0;
	}

////////////////////////
// Loading SIMC files //
////////////////////////
	// Normalization factor for runN (think about multiple files reading)
	cout << "SIMC: folder = " << SIMCfolder << endl;
	ifstream f_normfac("./files/simc/"+SIMCfolder+Form("/normfac_%d.dat",runN));
	Double_t normfac;
	f_normfac >> normfac;
	if(scaleSIMC) {
		cout << "********************" << endl
		     << "*** ATTENTION!!! ***  <- Scaling SIMC!" << endl
		     << "********************" << endl;
		normfac *= sSIMC;
	}
	cout << "SIMC: scaleSIMC*normfac = " << normfac << endl;
	if(logok) logFile << "SIMC.scale*normfac=" << normfac << ",";
	if(normfac==0.0) cout << "**************" << endl << "** WARNING! **  <- Did you remmember to run the files/simc/" << SIMCfolder << "/extractInfo.csh?" << endl << "*************" << endl;
	// Change
	ifstream f_simcQ("./files/simc/"+SIMCfolder+Form("/charge_%d.dat",runN));
	Double_t simcQ;
	f_simcQ >> simcQ;
	// Be careful, charge should be in miliC
	cout << "SIMC: charge = " << simcQ << " mC" << endl;
	if(logok) logFile << "SIMC.chargemC=" << simcQ << ",";
	// Ebeam
	ifstream f_simcEbeam("./files/simc/"+SIMCfolder+Form("/Ebeam_%d.dat",runN));
	Double_t simcEbeam;
	f_simcEbeam >> simcEbeam;
	cout << "SIMC: Ebeam = " << simcEbeam << endl;
	if(logok) logFile << "SIMC.Ebeam=" << simcEbeam << ",";
	// eTheta
	ifstream f_simceTh("./files/simc/"+SIMCfolder+Form("/eTheta_%d.dat",runN));
	Double_t simceTh;
	f_simceTh >> simceTh;
	cout << "SIMC: e_Theta = " << simceTh << endl;
	if(logok) logFile << "SIMC.eTheta=" << simceTh << ",";
	// eMomentum
	ifstream f_simceP("./files/simc/"+SIMCfolder+Form("/eP_%d.dat",runN));
	Double_t simceP;
	f_simceP >> simceP;
	cout << "SIMC: e_P = " << simceP << endl;
	if(logok) logFile << "SIMC.eP=" << simceP << ",";
	// pTheta
	ifstream f_simcpTh("./files/simc/"+SIMCfolder+Form("/pTheta_%d.dat",runN));
	Double_t simcpTh;
	f_simcpTh >> simcpTh;
	cout << "SIMC: p_Theta = " << simcpTh << endl;
	if(logok) logFile << "SIMC.pTheta=" << simcpTh << ",";
	// pMomentum
	ifstream f_simcpP("./files/simc/"+SIMCfolder+Form("/pP_%d.dat",runN));
	Double_t simcpP;
	f_simcpP >> simcpP;
	cout << "SIMC: p_P = " << simcpP << endl;
	if(logok) logFile << "SIMC.pP=" << simcpP << ",";

	TString inFolder = "./files/simc/"+SIMCfolder+"/";
	TString N, titles;
	if(runNfinal>runN) {
		N = outFolder+Form("/figures/%d/%d-%d",runN,runN,runNfinal); // name of output file
		Nfits = outFolder+Form("/figures/%d/backfits/%d-%d",runN,runN,runNfinal); // name of output file
		titles = Form("run%d-%d",runN,runNfinal);
	}
	else {
		N = outFolder+Form("/figures/%d/%d",runN,runN); // name of output file
		runNfinal=runN; //to avoid problems with TChain loading
		titles = Form("run%d",runN);
	}

	TChain *T = new TChain("h666");
	cout << "SIMC: Adding file " << runN << " - " << T->Add(inFolder + Form("marco_%d_",runN+i) + SIMtype + ".root") << endl;
	Double_t normSIMC = normfac/T->GetEntries()/simcQ;

///////////////////////
// Loading EXP files //
///////////////////////
	// For loading files
	TString EXinFolder = "./files/";//"./old_files_original_ntup_fpi2/";
	TString treename;
	if(EXtype.CompareTo("hms")==0) treename = "h9010";
	else if(EXtype.CompareTo("coin")==0) treename = "h9500";
	TChain *EX = new TChain(treename);
	Double_t *rawFileYields = new Double_t[runNfinal-runN+1];
	Int_t *FileNum = new Int_t[runNfinal-runN+1];
	Double_t *prescaler = new Double_t[runNfinal-runN+1];
	Double_t *track_eff = new Double_t[runNfinal-runN+1];
	Double_t *tr = new Double_t[runNfinal-runN+1];
	Double_t *ptr = new Double_t[runNfinal-runN+1];
	Double_t *hpre100 = new Double_t[runNfinal-runN+1];
	Double_t *hpre150 = new Double_t[runNfinal-runN+1];
	Double_t *helec = new Double_t[runNfinal-runN+1];
	Double_t *hcer_eEff = new Double_t[runNfinal-runN+1];
	Double_t *hcal_eEff = new Double_t[runNfinal-runN+1];
	Double_t *eff = new Double_t[runNfinal-runN+1];
	Double_t *norm = new Double_t[runNfinal-runN+1];
	Double_t *chargeA = new Double_t[runNfinal-runN+1];
	Double_t *chargeB = new Double_t[runNfinal-runN+1];
	Double_t *charge = new Double_t[runNfinal-runN+1];
	Double_t effCharge = 0.0;
	Int_t found;
	TFile *f;
	TString tempname;
	TTree *ttemp;
	cout << "EXP: Loading charges of runs (USING ONLY charge1):" << endl;
	Int_t nfiles=0;
	for(i=runN; i<=runNfinal; i++) {
		found=0;
		for(j=0;j<(sizeof(NotAnalyze)/sizeof(NotAnalyze[0])) && !found;j++)
			if(i==NotAnalyze[j]) found = 1;
		if(!found) {
			// Loading run files
			cout << endl << "----------" << endl;
			cout << "EXP: Adding run: " << i << endl;
			FileNum[nfiles] = i;
			rawFileYields[nfiles]=0.0;
			tempname = EXinFolder + EXtype + Form("%d",i) + ".root";
			if(!gSystem->AccessPathName(tempname)) {
				// to count raw number of events of individual files
				f = TFile::Open(tempname);
				ttemp = (TTree*)f->Get(treename);
				rawFileYields[nfiles] = rawFileYields[nfiles] + ttemp->GetEntries();
				// Chaining file
				cout << "\t- File: " << tempname << " - " << EX->Add(tempname) << " - " << ttemp->GetEntries() << " entries" << endl;
			}
			Int_t split=1;
			tempname = EXinFolder + EXtype + Form("%d.%d",i,split) + ".root";
			while(!gSystem->AccessPathName(tempname)) {
				// to count raw number of events of individual files
				f = TFile::Open(tempname);
				ttemp = (TTree*)f->Get(treename);
				rawFileYields[nfiles] = rawFileYields[nfiles] + ttemp->GetEntries();
				// Chaining file
				cout << "\t- File: " << tempname << " - " << EX->Add(tempname) << " - " << ttemp->GetEntries() << " entries" << endl;
				// Preparing next split
				tempname = EXinFolder + EXtype + Form("%d.%d",i,++split) + ".root";
			}
			// Reading efficiencies
			chargeA[nfiles] = ReadScaler("charge",i);
			chargeB[nfiles] = ReadScaler("charge2",i);
			// Attention! Charge should be in mC to match SIMC.
			chargeA[nfiles] /= 1000.0; // uC -> mC
			chargeB[nfiles] /= 1000.0; // uC -> mC
			charge[nfiles] = chargeA[nfiles]; // using 1st bcm only
			cout << Form("\t- charge1 = %.3f mC\t charge2 = %.3f mC\tCharge considered = %.3f milliCoulomb",i,chargeA[nfiles],chargeB[nfiles],charge[nfiles]) << endl;
			prescaler[nfiles] = ReadScaler("ps3",i);
			cout << "\t - prescaler: " << prescaler[nfiles] << endl;
			track_eff[nfiles] = ReadScaler("htrack",i);
			if(track_eff[nfiles]==0.0) {
				cout << "\t   * Warning: track_eff was ZERO! Setting it to one." << endl;
				track_eff[nfiles] = 1.0;
			}
			cout << "\t - track_eff: " << track_eff[nfiles] << endl;
			tr[nfiles] = ReadScaler("trig",i);
			cout << "\t - trigger eff: " << tr[nfiles] << endl;
			ptr[nfiles] = ReadScaler("ptrig",i);
			cout << "\t - pretrig_eff: " << ptr[nfiles] << endl;
			hpre100[nfiles] = ReadScaler("hpre100",i);
			cout << "\t - hpre100: " << hpre100[nfiles] << endl;
			hpre150[nfiles] = ReadScaler("hpre150",i);
			cout << "\t - hpre150: " << hpre150[nfiles] << endl;
			cout << "\t ---> Calculating:" << endl;
			helec[nfiles] = (60.0/50.0)*(hpre100[nfiles]-hpre150[nfiles])/hpre100[nfiles];
			helec[nfiles] = (60.0/50.0)*(hpre100[nfiles]-hpre150[nfiles])/hpre100[nfiles];
			cout << "\t\t - Electronic deadtime: " << helec[nfiles] << endl;
			eff[nfiles] = track_eff[nfiles]*(tr[nfiles]/ptr[nfiles])*(1.0-helec[nfiles]); // track_eff*(tr/ptr)*(1.0-helec)*(hcer_eEff*hcal_eEff)
			cout << "\t\t - Efficiencies: " << eff[nfiles] << endl;
			norm[nfiles] = (eff[nfiles]*charge[nfiles])/(prescaler[nfiles]);
			cout << "\t\t - Effective charge (Q*eff/prescaler): " << norm[nfiles] << endl;
			// Concluding... 
			rawFileYields[nfiles] /= charge[nfiles];
			cout << "\t - raw yield from this file = " << rawFileYields[nfiles] << " cnts/mC" << endl;
			effCharge += norm[nfiles];
			nfiles++;
		}
	}
	cout << endl << "EXP: Done reading files." << endl;
	// Normalization factor
	Double_t ExpNorm = (1.0/effCharge)*eNorm;
	cout << "EXP: eNorm = " << eNorm << " (manual normalization)" << endl;
	cout << "EXP: ExpNorm = eNorm*prescaler/(eff*charge) = " << ExpNorm << endl;
	// Checking consistency between rawFileYields
	Double_t meanRawYield = TMath::Mean(nfiles,rawFileYields);
	Double_t rmsRawYield = TMath::RMS(nfiles,rawFileYields);
	cout << "EXP: RawYield = " << meanRawYield << " +/- " << rmsRawYield << endl;
	for(Int_t i=0; i<nfiles; i++) {
		if(TMath::Abs(rawFileYields[i]-meanRawYield) > 2.0*rmsRawYield) {
			for(Int_t j=0; j<20; j++) cout << "*";
			cout << endl << "** WARNING!!! Run " << FileNum[i] << " with too different rawYield (greater than 2 sigmas)." << endl;
			cout << "** run " << FileNum[i] << ": RawYield = " << rawFileYields[i] << endl;
			cout << "** CAUTION HERE!" << endl;
			for(Int_t j=0; j<20; j++) cout << "*";
			cout << endl;
		}
	}

///////////////////
// Loading dummy //
///////////////////
	// Finding correspondent DUMMY, if not informed
	if(dummyRun==0 || dummyRunfinal<dummyRun) {
		ifstream f_dummyRun(Form("./files/flagsMarco/corrDummy_%d.dat",runN));
		f_dummyRun >> dummyRun;
		f_dummyRun >> dummyRunfinal;
	}
	if(dummyRunfinal==0) dummyRunfinal=dummyRun;
	if(dummyRun!=0) cout << Form("DUMMY: loading runs %d..%d as dummy",dummyRun,dummyRunfinal) << endl;
	else {
		cout << "ERROR: No dummy run defined. Does \"" << Form("./files/flagsMarco/corrDummy_%d.dat",runN) << "\" exists?" << endl;
		return 1;
	}
	if(logok) logFile << "DUM.run=" << dummyRun << ".." << dummyRunfinal << ",";
	TChain *DUM;
        if(EXtype.CompareTo("hms")==0) DUM = new TChain("h9010");
        else if(EXtype.CompareTo("coin")==0) DUM = new TChain("h9500");
	Double_t *DrawFileYields = new Double_t[dummyRunfinal-dummyRun+1];
	Int_t *DFileNum = new Int_t[dummyRunfinal-dummyRun+1];
	Double_t *Dprescaler = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t *Dtrack_eff = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t *Dtr = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t *Dptr = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t *Dhpre100 = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t *Dhpre150 = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t *Dhelec = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t *Dhcer_eEff = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t *Dhcal_eEff = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t *Deff = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t *Dnorm = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t *DchargeA = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t *DchargeB = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t *Dcharge = new Double_t[dummyRunfinal-dummyRun+1];
	Double_t DeffCharge = 0.0;
	TString tempname;
	cout << "DUMMY: Loading charges of runs (USING ONLY charge1):" << endl;
	Int_t Dnfiles=0;
	for(i=dummyRun; i<=dummyRunfinal; i++) {
		found=0;
		for(j=0;j<(sizeof(NotAnalyze)/sizeof(NotAnalyze[0])) && !found;j++)
			if(i==NotAnalyze[j]) found = 1;
		if(!found) {
			// Loading run files
			cout << endl << "----------" << endl;
			cout << "DUMMY: Adding run: " << i << endl;
			DFileNum[Dnfiles] = i;
			DrawFileYields[Dnfiles]=0.0;
			tempname = EXinFolder + EXtype + Form("%d",i) + ".root";
			if(!gSystem->AccessPathName(tempname)) {
				// to count raw number of events of individual files
				f = TFile::Open(tempname);
				ttemp = (TTree*)f->Get(treename);
				DrawFileYields[Dnfiles] = DrawFileYields[Dnfiles] + ttemp->GetEntries();
				// Chaining file
				cout << "\t- File: " << tempname << " - " << DUM->Add(tempname) << " - " << ttemp->GetEntries() << " entries" << endl;
			}
			Int_t split=1;
			tempname = EXinFolder + EXtype + Form("%d.%d",i,split) + ".root";
			while(!gSystem->AccessPathName(tempname)) {
				// to count raw number of events of individual files
				f = TFile::Open(tempname);
				ttemp = (TTree*)f->Get(treename);
				DrawFileYields[Dnfiles] = DrawFileYields[Dnfiles] + ttemp->GetEntries();
				// Chaining file
				cout << "\t- File: " << tempname << " - " << DUM->Add(tempname) << " - " << ttemp->GetEntries() << " entries" << endl;
				// Preparing next split
				tempname = EXinFolder + EXtype + Form("%d.%d",i,++split) + ".root";
			}
			// Reading efficiencies
			DchargeA[Dnfiles] = ReadScaler("charge",i);
			DchargeB[Dnfiles] = ReadScaler("charge2",i);
			// Attention! Charge should be in mC to match SIMC.
			DchargeA[Dnfiles] /= 1000.0; // uC -> mC
			DchargeB[Dnfiles] /= 1000.0; // uC -> mC
			Dcharge[Dnfiles] = DchargeA[Dnfiles]; // using 1st bcm only
			cout << Form("\t- charge1 = %.3f mC\t charge2 = %.3f mC\tCharge considered = %.3f milliCoulomb",i,DchargeA[Dnfiles],DchargeB[Dnfiles],Dcharge[Dnfiles]) << endl;
			Dprescaler[Dnfiles] = ReadScaler("ps3",i);
			cout << "\t - prescaler: " << Dprescaler[Dnfiles] << endl;
			Dtrack_eff[Dnfiles] = ReadScaler("htrack",i);
			if(Dtrack_eff[Dnfiles]==0.0) {
				cout << "\t   * Warning: track_eff was ZERO! Setting it to one." << endl;
				Dtrack_eff[Dnfiles] = 1.0;
			}
			cout << "\t - track_eff: " << Dtrack_eff[Dnfiles] << endl;
			Dtr[Dnfiles] = ReadScaler("trig",i);
			cout << "\t - trigger eff: " << Dtr[Dnfiles] << endl;
			Dptr[Dnfiles] = ReadScaler("ptrig",i);
			cout << "\t - pretrig_eff: " << Dptr[Dnfiles] << endl;
			Dhpre100[Dnfiles] = ReadScaler("hpre100",i);
			cout << "\t - hpre100: " << Dhpre100[Dnfiles] << endl;
			Dhpre150[Dnfiles] = ReadScaler("hpre150",i);
			cout << "\t - hpre150: " << Dhpre150[Dnfiles] << endl;
			cout << "\t ---> Calculating:" << endl;
			Dhelec[Dnfiles] = (60.0/50.0)*(Dhpre100[Dnfiles]-Dhpre150[Dnfiles])/Dhpre100[Dnfiles];
			Dhelec[Dnfiles] = (60.0/50.0)*(Dhpre100[Dnfiles]-Dhpre150[Dnfiles])/Dhpre100[Dnfiles];
			cout << "\t\t - Electronic deadtime: " << Dhelec[Dnfiles] << endl;
			Deff[Dnfiles] = Dtrack_eff[Dnfiles]*(Dtr[Dnfiles]/Dptr[Dnfiles])*(1.0-Dhelec[Dnfiles]); // Dtrack_eff*(Dtr/Dptr)*(1.0-Dhelec)*(Dhcer_eEff*Dhcal_eEff)
			cout << "\t\t - Efficiencies: " << Deff[Dnfiles] << endl;
			Dnorm[Dnfiles] = (Deff[Dnfiles]*Dcharge[Dnfiles])/(Dprescaler[Dnfiles]);
			cout << "\t\t - Effective charge (Q*eff/prescaler): " << Dnorm[Dnfiles] << endl;
			// Concluding... 
			DrawFileYields[Dnfiles] /= Dcharge[Dnfiles];
			cout << "\t - raw yield from this file = " << DrawFileYields[Dnfiles] << " cnts/mC" << endl;
			DeffCharge += Dnorm[Dnfiles];
			Dnfiles++;
		}
	}
	cout << endl << "DUM: Done reading files." << endl;
	// Normalization factor
	Double_t ExpNormD = (1.0/DeffCharge)*eNorm;
	cout << "DUMMY: eNorm = " << eNorm << " (manual normalization)" << endl;
	cout << "DUMMY: ExpNorm = eNorm*prescaler/(eff*charge) = " << ExpNormD << endl;
	// Checking consistency between DrawFileYields
	if(Dnfiles>1) {
		Double_t DmeanRawYield = TMath::Mean(Dnfiles,DrawFileYields);
		Double_t DrmsRawYield = TMath::RMS(Dnfiles,DrawFileYields);
		cout << "DUMMY: RawYield = " << DmeanRawYield << " +/- " << DrmsRawYield << endl;
		for(Int_t i=0; i<Dnfiles; i++) {
			if(TMath::Abs(DrawFileYields[i]-DmeanRawYield) > 2.0*DrmsRawYield) {
				for(Int_t j=0; j<20; j++) cout << "*";
				cout << endl << "** WARNING!!! DUMMY run " << FileNum[i] << " with too different rawYield (greater than 2 sigmas)." << endl;
				cout << "** DUMMY run " << DFileNum[i] << ": RawYield = " << DrawFileYields[i] << endl;
				cout << "** CAUTION HERE!" << endl;
				for(Int_t j=0; j<20; j++) cout << "*";
				cout << endl;
			}
		}
	}
return 2;

///////////////////
// Defining CUTs //
///////////////////
	//T->StartViewer();

	//////////////////
	// Elastic CUTS //
	//////////////////
	// HMS
	TCut chcer_el = "hcer_npe>0.3"; // gas cerenkov
	TCut chcal_el = "hsshsum>0.5"; // calorimeter
//	TCut chcal_max = "hsshsum<1.15"; // calorimeter
//	TCut chaero = ; // 
//	TCut chbeta = "abs(hsbeta-0.9585)<0.1"; // p/E
//	TCut chaero = haero_ne;
	//TCut cW_el = "abs(W-0.938272)<(3.0*1.09186e-02)"; //<- ATTENTION HERE!!!!
	TCut cW_el = "abs(W-0.938272)<0.1"; //<- ATTENTION HERE!!!!
	//TCut cW_el = "";//"abs(W-0.938272)<0.1"; //<- ATTENTION HERE!!!!
	TCut cWmax_el = "W < 1.07";
	TCut cdelta_el = "abs(hsdelta)<8.5";
	TCut cxptar_el = "abs(hsxptar)<0.09";
	TCut cyptar_el = "abs(hsyptar)<0.055";
	TCut cWbackFit_el = "abs(W-0.96)>(6.0*1.09186e-02)";
	Double_t limitInt[2]; // for background fit
	limitInt[0] = 0.938272 - 3.0*1.09186e-02;
	limitInt[1] = 0.938272 + 3.0*1.09186e-02;

	///////////////
	// Kaon CUTS //
	///////////////
	// SOS
//	TCut cssp = // momentum
	TCut cssshsum = "ssshsum > 0.7";// calorimeter
	TCut cscer = "scer_npe>0.5";// gas cerenkov
	TCut cssdelta = "abs(ssdelta)<20.0";

	// HMS
	TCut chcer_k = "hcer_npe<2.0"; // gas cerenkov
	// Threshold momentum
	// pi+:  3840.90 MeV/c
	// K+:  13585.77 MeV/c
	// p+:  25820.83 MeV/c
	TCut chcal_k = "";//"hsshsum>0.5"; // calorimeter
//	TCut chcal_max_k = "hsshsum<1.15"; // calorimeter
	TCut chaero_k = "haero_su > 1.5"; 
	// Thresholds (M/sqrt(1.03^2-1)):
	// - pi+:  565.566 MeV/c
	// - K+:  2000.480 MeV/c
	// - p+:  3802.070 MeV/c
//	TCut chbeta_k = "abs(hsbeta-0.9585)<0.1"; // p/E
	//TCut cW_k = "abs(W-0.938272)<(3.0*1.09186e-02)"; //<- ATTENTION HERE!!!!
//	TCut cW_k = "abs(W-0.938272)<0.1"; //<- ATTENTION HERE!!!!
	//TCut cW_k = "";//"abs(W-0.938272)<0.1"; //<- ATTENTION HERE!!!!
//	TCut cWmax_k = "W < 1.07";
	TCut cdelta_k = "abs(hsdelta)<8.5";
	TCut cxptar_k = "abs(hsxptar)<0.09";
	TCut cyptar_k = "abs(hsyptar)<0.055";

	// t
	TCut tCut = Form("t>%f && t<%f",t_range[0],t_range[1]);

	// ExtraCUT
	TCut extraCut = "";

	// W vs Q2 cut
	Double_t p0_WQ2cut[4], p1_WQ2cut[4];
	for(Int_t i=0; i<4; i++) {
		p1_WQ2cut[i] = (ranges.WQ2cutQ2[(i+1)%4]-ranges.WQ2cutQ2[i])/(ranges.WQ2cutW[(i+1)%4]-ranges.WQ2cutW[i]);
		p0_WQ2cut[i] = ranges.WQ2cutQ2[i] - p1_WQ2cut[i]*ranges.WQ2cutW[i];
	}
	TCut cWQ2_k = Form("(Q2>(%f*W+%f)) && (Q2<(%f*W+%f)) && (Q2<(%f*W+%f)) && (Q2>(%f*W+%f))",p1_WQ2cut[0],p0_WQ2cut[0],p1_WQ2cut[1],p0_WQ2cut[1],p1_WQ2cut[2],p0_WQ2cut[2],p1_WQ2cut[3],p0_WQ2cut[3]);

	// REACTION
	//missmass
	Double_t cmissmass1_lim[2]; // Kaons!
	cmissmass1_lim[0] = 1.115683-0.020;
	cmissmass1_lim[1] = 1.115683+0.032;
	Double_t cmissmass1_center = (cmissmass1_lim[1]+cmissmass1_lim[0])/2.0; //1.115683; //Lambda
	Double_t cmissmass1_width  = (cmissmass1_lim[1]-cmissmass1_lim[0])/2.0; //2.0*0.0125
	TCut cmissmass1 = Form("abs(missmass-%f)<%f",cmissmass1_center,cmissmass1_width);

	Double_t cmissmass1err_center[2];
	cmissmass1err_center[0] = cmissmass1_center-1.5*cmissmass1_width;
	cmissmass1err_center[1] = cmissmass1_center+1.5*cmissmass1_width;
	Double_t cmissmass1err_width[2];
	cmissmass1err_width[0] = 1.5*cmissmass1_width;
	cmissmass1err_width[1] = 1.5*cmissmass1_width;
	TCut cmissmass1err = Form("(abs(missmass-%f)<%f || abs(missmass-%f)<%f)",cmissmass1err_center[0],cmissmass1err_width[0],cmissmass1err_center[1],cmissmass1err_width[1]);
//	TCut cQ2_k = Q2
	// cointime
	Double_t coin1_lim[2] = {-15.0,15.0}; // coarse time cut (acceptance)
	Double_t coin1_center = (coin1_lim[1]+coin1_lim[0])/2.0;
	Double_t coin1_width =  (coin1_lim[1]-coin1_lim[0])/2.0;
	TCut ccointime1 = Form("abs(cointime-%f)<%f",coin1_center,coin1_width);
	Double_t coin2_lim[2] = {-0.6,0.9}; // Kaons!
	// if change here, make sure to match the background subtraction binning to have the same range
	Double_t coin2_center = (coin2_lim[1]+coin2_lim[0])/2.0;
	Double_t coin2_width =  (coin2_lim[1]-coin2_lim[0])/2.0;
	TCut ccointime2 = Form("abs(cointime-%f)<%f",coin2_center,coin2_width);

	TCut coinback = Form("abs(cointime-%f+3.0*2.0)<%f || abs(cointime-%f+4.0*2.0)<%f || abs(cointime-%f+5.0*2.0)<%f",coin2_center,coin2_width,coin2_center,coin2_width,coin2_center,coin2_width);

	/////////
	// ALL //
	/////////
	TCut call1, call1backFit, SIMCcall1, call_acept, call2back;
	// Elastic
	if(type==1) {
		call1 = cdelta_el && cxptar_el && cyptar_el && chcer_el && chcal_el && cW_el && cWmax_el; 
		call_acept = cdelta_el && cxptar_el && cyptar_el; // cutting only acceptance
		call1backFit = cdelta_el && cxptar_el && cyptar_el && chcer_el && chcal_el && cWbackFit_el; // to be used for W background fit
		SIMCcall1 = cdelta_el && cxptar_el && cyptar_el && cW_el && cWmax_el;
	}
	else if(type==2) { // KAON
		call1 = cWQ2_k && cssdelta && cssshsum && cscer && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k && cmissmass1 && ccointime2 && chaero_k;// && chsmass2;
		call2 = cWQ2_k && cssdelta && cssshsum && cscer && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && cmissmass1 && ccointime2 && tCut && extraCut; // kaons only
		call2tmp = cWQ2_k && cssdelta && cssshsum && cscer && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && cmissmass1 && ccointime2; // kaons only with NO extraCut
		call2back = cWQ2_k && cssdelta && cssshsum && cscer && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && cmissmass1 && coinback;
		call_acept = cWQ2_k && cssdelta && cssshsum && cscer && cdelta_k && cxptar_k && cyptar_k && ccointime1 && chaero_k;// && chsmass2;
		call1backFit = cWQ2_k && cssdelta && cssshsum && cscer && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k && cmissmass1 && ccointime2; // to be used for W background fit
		call1err = cWQ2_k && cssdelta && cssshsum && cscer && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k && cmissmass1err && ccointime2; // to be used for background analysis around call1
		SIMCcall1 = cWQ2_k && cssdelta && cdelta_k && cxptar_k && cyptar_k && cmissmass1 && tCut && extraCut;
	}

	if(logok) logFile << "EXP_DUM.CUTS=" << call1.GetTitle() << ",";
	if(logok) logFile << "SIMC.CUTS=" << SIMCcall1.GetTitle() << ",";

//********************************
//********************************
//********************************

if(1==0) {
	cout << "c0" << endl;
	// EXP: hcer_npe vs calorimeter
	TCanvas *c0 = new TCanvas("c0","c0",600,600);
	MConfigCanvas(c0,0,0);
	//-> EXP.
	//TH1F *phcernpe_all = new TH1F("phcernpe_all","Delta",100,-15.0,10.0);
	EX->Draw("hcer_npe:hsshsum",call_acept,"goff");
	TGraph *gr0_all = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	MConfigPoints(gr0_all,4,0.5,8);
	EX->Draw("hcer_npe:hsshsum",call1,"goff");
	TGraph *gr0_cut = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	MConfigPoints(gr0_cut,3,0.5,8);
	gr0_all->GetXaxis()->SetTitle("NPE Gas Cherenkov (hcer_npe)");
	gr0_all->GetYaxis()->SetTitle("Calorimeter (hsshsum)");
	gr0_all->SetTitle(titles + " - EXP: Gas Cherenkov vs. Calorimeter");
	gr0_all->Draw("ap");
	gr0_cut->Draw("psame");
	//// Other things...
	Double_t c0xlim[2], c0ylim[2];
	c0->Update();
	c0xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c0xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c0ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c0ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l0;
	l0.SetTextSize(0.04);
	l0.DrawLatex(c0xlim[0]+0.7*(c0xlim[1]-c0xlim[0]),c0ylim[0]+0.7*(c0ylim[1]-c0ylim[0]),"#color[4]{All data}");
	l0.DrawLatex(c0xlim[0]+0.7*(c0xlim[1]-c0xlim[0]),c0ylim[0]+0.65*(c0ylim[1]-c0ylim[0]),"#color[3]{After cut}");
	l0.Draw();
	pad2png(c0,N+"_hs.png");
	if(savePDF) c0->Print(N+".pdf(");
}

//********************************
//********************************
//********************************

if(1==0) {
	cout << "c1" << endl;
	// hsdelta
	TCanvas *c1 = new TCanvas("c1","c1",600,600);
	MConfigCanvas(c1,0,0);
	//-> EXP.
	TH1F *phsdelta_all = new TH1F("phsdelta_all","Delta",90,-20.0,20.0);
	MConfigHist(phsdelta_all,"Momentum Delta","Delta","Weighted counts",41,0);
	EX->Draw("hsdelta>>+phsdelta_all",call2 * Form("%g",ExpNorm));//call_acept
	phsdelta_all->GetYaxis()->SetRangeUser(0.0,0.12);
	phsdelta_all->GetXaxis()->SetRangeUser(-20.0,20.0);
	phsdelta_all->SetDirectory(0);
	TH1F *phsdelta = new TH1F("phsdelta","",90,-20.0,20.0);
	MConfigHist(phsdelta,"","","",3,0);
	EX->Draw("hsdelta>>+phsdelta",call2 * Form("%g",ExpNorm),"goff");
	phsdelta->SetDirectory(0);
	TH1F *phsdelta_noWeight = new TH1F("phsdelta_noWeight","",90,-20.0,20.0);
	EX->Draw("hsdelta>>+phsdelta_noWeight",call2,"goff");
	//->Background from shift on cointime
	TH1F *phsdeltaBack = new TH1F("phsdeltaBack","phsdeltaBack",90,-20.0,20.0);
	EX->Draw("hsdelta>>+phsdeltaBack",call2back * Form("%g*%g",1.0/3.0,ExpNorm),"same");
	MConfigHist(phsdeltaBack,"","","",6,0);
	phsdeltaBack->SetDirectory(0);
	//-> DUMMY
	//TH1F *phsdeltaD_all = new TH1F("phsdeltaD_all","",90,-20.0,20.0);
	//MConfigHist(phsdeltaD_all,"","","",1,0);
	//DUM->Draw("hsdelta>>+phsdeltaD_all",Form("%g",call_acept * ExpNormD));
	//phsdeltaD_all->SetDirectory(0);
	TH1F *phsdeltaD = new TH1F("phsdeltaD","DeltaCut",90,-20.0,20.0);
	MConfigHist(phsdeltaD,"","","",8,0);
	DUM->Draw("hsdelta>>+phsdeltaD",call2 * Form("%g",ExpNormD),"goff");
	phsdeltaD->SetDirectory(0);
	TH1F *phsdeltaD_noWeight = new TH1F("phsdeltaD_noWeight","",90,-20.0,20.0);
	DUM->Draw("hsdelta>>+phsdeltaD_noWeight",call2,"goff");
	//-> EXP - DUMMY
	TH1F *edhsdelta = phsdelta->Clone();
	edhsdelta->SetName("edhsdelta");
	edhsdelta->Add(phsdeltaD,-1*SubtractDummy);
	MConfigHist(edhsdelta,"","","",1,0);
	edhsdelta->Draw("same");
	//-> SIMC
	TH1F *mhsdelta_all = new TH1F("mhsdelta_all","Delta",90,-20.0,20.0);
	MConfigHist(mhsdelta_all,titles+"","","",30,0);
	T->Draw("hsdelta>>+mhsdelta_all",call_acept * Form("Weight*%g",normSIMC),"goff");//,"","");
	mhsdelta_all->SetDirectory(0);
	TH1F *mhsdelta = new TH1F("mhsdelta","DeltaCut",90,-20.0,20.0);
	MConfigHist(mhsdelta,"","","",2,0);
	//mhcer->SetFillStyle(3004)
	T->Draw("hsdelta>>+mhsdelta",SIMCcall1 * Form("(Weight*%g)",normSIMC),"same");
	mhsdelta->SetDirectory(0);
	//// Other things...
	Double_t c1xlim[2], c1ylim[2];
	c1->Update();
	c1xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c1xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c1ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c1ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l1;
	l1.SetTextSize(0.04);
	l1.DrawLatex(c1xlim[0]+0.05*(c1xlim[1]-c1xlim[0]),c1ylim[0]+0.95*(c1ylim[1]-c1ylim[0]),"#color[2]{SIMC (weighted + cut)}");
	l1.DrawLatex(c1xlim[0]+0.05*(c1xlim[1]-c1xlim[0]),c1ylim[0]+0.9*(c1ylim[1]-c1ylim[0]),"#color[1]{Experimental - Dummy (weighted + cut)}");
	l1.DrawLatex(c1xlim[0]+0.05*(c1xlim[1]-c1xlim[0]),c1ylim[0]+0.85*(c1ylim[1]-c1ylim[0]),"#color[6]{Exp. - Background from shift on cointime}");
	//l1.DrawLatex(c1xlim[0]+0.05*(c1xlim[1]-c1xlim[0]),c1ylim[0]+0.8*(c1ylim[1]-c1ylim[0]),"#color[30]{SIMC (weighted)}");
	//l1.DrawLatex(c1xlim[0]+0.05*(c1xlim[1]-c1xlim[0]),c1ylim[0]+0.75*(c1ylim[1]-c1ylim[0]),"#color[41]{EXP (weighted)}");
	//l1.DrawLatex(c1xlim[0]+0.05*(c1xlim[1]-c1xlim[0]),c1ylim[0]+0.7*(c1ylim[1]-c1ylim[0]),"#color[8]{DUMMY (weighted + cut)}");
	cout << "EXP:  Selected/Total entries: " << EX->GetSelectedRows() << "/" << EX->GetEntries() << Form(" (%.1f\%)",100.0*EX->GetSelectedRows()/EX->GetEntries()) << endl
	     << "DUMMY:  Selected/Total entries: " << DUM->GetSelectedRows() << "/" << DUM->GetEntries() << Form(" (%.1f\%)",100.0*DUM->GetSelectedRows()/DUM->GetEntries()) << endl
	     << "SIMC: Selected/Total entries:" << T->GetSelectedRows() << "/" <<  T->GetEntries() << Form(" (%.1f\%)",100.0*T->GetSelectedRows()/T->GetEntries()) << endl;
	//cout << Form("%.3f\t%.3f",c1ylim[0],c1ylim[1]) << endl;
	///////////////////////
	// CALCULATING RATIO //
	///////////////////////
	Double_t limitDeltaInt[] = {-20.0,20.0};//{-5.0,5.0};
	// Integrating EXP-DUMMY
        Double_t integral_sub = MHistIntegrate(edhsdelta,limitDeltaInt[0],limitDeltaInt[1]);
	if(logok) logFile << "EXP.yield=" << MHistIntegrate(phsdelta,limitDeltaInt[0],limitDeltaInt[1]) - MHistIntegrate(phsdeltaD,limitDeltaInt[0],limitDeltaInt[1]);
        cout	<< "For ratio analysis:";
	if(!SubtractDummy) cout << "<Marco> ATTENTION! Dummy NOT subtracted.";
	cout << endl << "    EXP-DUMMY integral: " << integral_sub << " +/- ";
	// Uncertainties EXP-DUMMY
	Double_t integral_Exp_noWeight = MHistIntegrate(phsdelta_noWeight,limitDeltaInt[0],limitDeltaInt[1]);
	Double_t integral_Dum_noWeight = MHistIntegrate(phsdeltaD_noWeight,limitDeltaInt[0],limitDeltaInt[1]);
	Double_t sig_sub = sqrt(pow(ExpNorm,2.0)*integral_Exp_noWeight+pow(ExpNormD,2.0)*integral_Dum_noWeight);
	cout << sig_sub << Form("  -> (%.2f\%)",100.0*sig_sub/integral_sub) << endl;
	// Integrating SIMC
        Double_t integral_SIMC = MHistIntegrate(mhsdelta,limitDeltaInt[0],limitDeltaInt[1]);
	Double_t sigma_Ratio = (sig_sub/integral_sub)*integral_sub/integral_SIMC;
	cout	<< "    SIMC integral: " << integral_SIMC << " +/- ???" << endl
		<< "    RATIO: " << integral_sub/integral_SIMC << " +/- " << sigma_Ratio << endl
		<< "-----------------------------------" << endl;
	// Saving in log file
	if(logok) logFile << "EXP.yield=" << MHistIntegrate(phsdelta,limitDeltaInt[0],limitDeltaInt[1]) << ",DUM.yield=" << MHistIntegrate(phsdeltaD,limitDeltaInt[0],limitDeltaInt[1]) << ",EXP-DUM=" << integral_sub << ",s(EXP-DUM)=" << sig_sub << ",SIMC.yield=" << integral_SIMC << ",RATIO=" << integral_sub/integral_SIMC << ",s(RATIO)=" << sigma_Ratio << ",";
	// Printing
	pad2png(c1,N+"_hsdelta.png");
	if(savePDF) c1->Print(N+".pdf");
}

//********************************
//********************************
//********************************

if(1==1) {
	cout << "c1a - HMS delta plot" << endl;
	// hsdelta_allPlots
	TCanvas *c1a = new TCanvas("c1a","c1a",1800,600);
	//c1a->Divide(3,1);
	TPad *p1a1 = c1a->cd(1);
	p1a1->SetGrid();
	p1a1->GetFrame()->SetBorderSize(10);
	//-> EXP. - no hsdelta cut
	TH1F *phsdelta_ano = new TH1F("phsdelta_ano","Delta",150,-20.0,20.0);
	MConfigHist(phsdelta_ano,"HMS Delta","hsdelta","Yield",14,0);
	TCut c1a_cutEXP = cWQ2_k && cssdelta && cssshsum && cscer && cxptar_k && cyptar_k && chcer_k && chcal_k && chaero_k && cmissmass1 && ccointime2 && tCut && extraCut;
	EX->Draw("hsdelta>>+phsdelta_ano", c1a_cutEXP * Form("%g",ExpNorm));
	phsdelta_ano->SetMaximum(1.4*phsdelta_ano->GetBinContent(phsdelta_ano->GetMaximumBin()));
	// phsdelta_ano->GetYaxis()->SetRangeUser(0.0,14000.0);
	phsdelta_ano->SetDirectory(0);
	//-> EXP.
	TH1F *phsdelta_a = new TH1F("phsdelta_a","Delta",150,-20.0,20.0);
	MConfigHist(phsdelta_a,"","","",4,0);
	EX->Draw("hsdelta>>+phsdelta_a",call2 * Form("%g",ExpNorm),"same");
	// phsdelta_a->GetYaxis()->SetRangeUser(0.0,14000.0);
	phsdelta_a->SetDirectory(0);
	//-> DUMMY
	//TPad *p1a2 = c1a->cd(2);
	//p1a2->SetGrid();
	//p1a2->GetFrame()->SetBorderSize(10);
	TH1F *phsdeltaD_a = new TH1F("phsdeltaD_a","Delta",150,-20.0,20.0);
	MConfigHist(phsdeltaD_a,"Delta - Dummy","Delta","Counts",1,0);
	DUM->Draw("hsdelta>>+phsdeltaD_a",call2 * Form("%g",ExpNormD),"same");
	phsdeltaD_a->SetDirectory(0);
	//-> SIMC
	//TPad *p1a3 = c1a->cd(3);
	//p1a3->SetGrid();
	//p1a3->GetFrame()->SetBorderSize(10);
	TH1F *mhsdelta_a = new TH1F("mhsdelta_a","Delta",150,-20.0,20.0);
	MConfigHist(mhsdelta_a,"Delta - SIMC","Delta","Counts",3,0);
	T->Draw("hsdelta>>+mhsdelta_a",SIMCcall1 * Form("Weight*%g",normSIMC),"same");
	mhsdelta_a->SetDirectory(0);
	//// Other things...
	c1a->Update();
	Double_t c1axlim[2], c1aylim[2];
	c1axlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c1axlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c1aylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c1aylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l1a;
	l1a.SetTextSize(0.03);
	l1a.DrawLatex(c1axlim[0]+0.05*(c1axlim[1]-c1axlim[0]),c1aylim[0]+0.95*(c1aylim[1]-c1aylim[0]),"#color[14]{EXP - All cuts, but hsdelta}");
	l1a.DrawLatex(c1axlim[0]+0.05*(c1axlim[1]-c1axlim[0]),c1aylim[0]+0.90*(c1aylim[1]-c1aylim[0]),"#color[4]{EXP - All cuts}");
	l1a.DrawLatex(c1axlim[0]+0.05*(c1axlim[1]-c1axlim[0]),c1aylim[0]+0.85*(c1aylim[1]-c1aylim[0]),"#color[1]{Dummy}");
	l1a.DrawLatex(c1axlim[0]+0.05*(c1axlim[1]-c1axlim[0]),c1aylim[0]+0.80*(c1aylim[1]-c1aylim[0]),"#color[3]{SIMC}");
	l1a.Draw();
	// Printing
	pad2png(c1a,N+"_hsdelta_steps.png");
	if(savePDF) c1a->Print(N+".pdf");
	//if(fast) return 1;
}

//********************************
//********************************
//********************************

if(1==1) { // SOS delta
	cout << "c1b - SOS delta plot" << endl;
	// hsdelta_allPlots
	TCanvas *c1b = new TCanvas("c1b","c1b",1800,600);
	//c1b->Divide(3,1);
	TPad *p1b1 = c1b->cd(1);
	p1b1->SetGrid();
	p1b1->GetFrame()->SetBorderSize(10);
	//-> EXP. - no ssdelta cut
	TH1F *pssdelta_ano = new TH1F("pssdelta_ano","SOS Delta",150,-20.0,20.0);
	MConfigHist(pssdelta_ano,"SOS delta","ssdelta","Yield",14,0);
	TCut c1b_cutEXP = cWQ2_k && cdelta_k && cssshsum && cscer && cxptar_k && cyptar_k && chcer_k && chcal_k && chaero_k && cmissmass1 && ccointime2 && tCut && extraCut;
	EX->Draw("ssdelta>>+pssdelta_ano", c1b_cutEXP * Form("%g",ExpNorm));
	pssdelta_ano->SetMaximum(1.4*pssdelta_ano->GetBinContent(pssdelta_ano->GetMaximumBin()));
	// pssdelta_ano->GetYaxis()->SetRangeUser(0.0,14000.0);
	pssdelta_ano->SetDirectory(0);
	//-> EXP.
	TH1F *pssdelta_a = new TH1F("pssdelta_a","Delta",150,-20.0,20.0);
	MConfigHist(pssdelta_a,"","","",4,0);
	EX->Draw("ssdelta>>+pssdelta_a",call2 * Form("%g",ExpNorm),"same");
	// pssdelta_a->GetYaxis()->SetRangeUser(0.0,14000.0);
	pssdelta_a->SetDirectory(0);
	//-> DUMMY
	//TPad *p1b2 = c1b->cd(2);
	//p1b2->SetGrid();
	//p1b2->GetFrame()->SetBorderSize(10);
	TH1F *pssdeltaD_a = new TH1F("pssdeltaD_a","Delta",150,-20.0,20.0);
	MConfigHist(pssdeltaD_a,"Delta - Dummy","Delta","Counts",1,0);
	DUM->Draw("ssdelta>>+pssdeltaD_a",call2 * Form("%g",ExpNormD),"same");
	pssdeltaD_a->SetDirectory(0);
	//-> SIMC
	//TPad *p1b3 = c1b->cd(3);
	//p1b3->SetGrid();
	//p1b3->GetFrame()->SetBorderSize(10);
	TH1F *mssdelta_a = new TH1F("mssdelta_a","Delta",150,-20.0,20.0);
	MConfigHist(mssdelta_a,"Delta - SIMC","Delta","Counts",3,0);
	T->Draw("ssdelta>>+mssdelta_a",SIMCcall1 * Form("Weight*%g",normSIMC),"same");
	mssdelta_a->SetDirectory(0);
	//// Other things...
	c1b->Update();
	Double_t c1bxlim[2], c1bylim[2];
	c1bxlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c1bxlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c1bylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c1bylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l1b;
	l1b.SetTextSize(0.03);
	l1b.DrawLatex(c1bxlim[0]+0.05*(c1bxlim[1]-c1bxlim[0]),c1bylim[0]+0.95*(c1bylim[1]-c1bylim[0]),"#color[14]{EXP - All cuts, but ssdelta}");
	l1b.DrawLatex(c1bxlim[0]+0.05*(c1bxlim[1]-c1bxlim[0]),c1bylim[0]+0.90*(c1bylim[1]-c1bylim[0]),"#color[4]{EXP - All cuts}");
	l1b.DrawLatex(c1bxlim[0]+0.05*(c1bxlim[1]-c1bxlim[0]),c1bylim[0]+0.85*(c1bylim[1]-c1bylim[0]),"#color[1]{Dummy}");
	l1b.DrawLatex(c1bxlim[0]+0.05*(c1bxlim[1]-c1bxlim[0]),c1bylim[0]+0.80*(c1bylim[1]-c1bylim[0]),"#color[3]{SIMC}");
	l1b.Draw();
	// Printing
	pad2png(c1b,N+"_ssdelta_steps.png");
	if(savePDF) c1b->Print(N+".pdf");
	//if(fast) return 1;
}

//********************************
//********************************
//********************************

	// cointime vs beta
if(type==2 && 1==0) { // coincidence
	cout << "c2 - 3D cointime vs beta plot" << endl;
	TCanvas *c2 = new TCanvas("c2","c2",600,600);
	MConfigCanvas(c2,0,0);
	EX->Draw("cointime:hsbeta",call_acept,"goff");
	TGraph *gr2_all = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	MConfigPoints(gr2_all,4,0.4,8);
	gr2_all->GetXaxis()->SetRangeUser(-30.0,10.0);
	gr2_all->GetYaxis()->SetRangeUser(-0.1,1.4);
	gr2_all->SetTitle("Beam structure");
	gr2_all->GetXaxis()->SetTitle("Coincidence Time / ns");
	gr2_all->GetYaxis()->SetTitle("#beta");
	EX->Draw("cointime:hsbeta",call1,"goff");
	TGraph *gr2 = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	MConfigPoints(gr2,3,0.4,8);
	gr2_all->Draw("ap");
	gr2->Draw("psame");
	c2->Update();
	Double_t c2xlim[2], c2ylim[2];
	c2xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c2xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c2ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c2ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l2;
	l2.SetTextSize(0.04);
	//l2.DrawLatex(c2xlim[0]+0.5*(c2xlim[1]-c2xlim[0]),c2ylim[0]+0.3*(c2ylim[1]-c2ylim[0]),"#color[4]{All data}");
	//l2.DrawLatex(c2xlim[0]+0.5*(c2xlim[1]-c2xlim[0]),c2ylim[0]+0.25*(c2ylim[1]-c2ylim[0]),"#color[3]{Cut center structure}");
	l2.Draw();
	c2->Update();
	pad2png(c2,N+"_hsbeta_cointime.png");
	c2->Print(N+".pdf");
}

//********************************
//********************************
//********************************
	// hsbeta vs cointime 2D-histogram ZOOM
if(type==2 && 1==0) {
	TCanvas *c2aa = new TCanvas("c2aa","c2aa",1200,600);
	MConfigCanvas(c2aa,0,0);
	TH2F *h2aa = new TH2F("h2aa","h2aa",75,-1.0,2.0,75,0.7,1.3);
	gStyle->SetNumberContours(250);
	h2aa->GetZaxis()->SetRangeUser(0.0, 200.0);
	gStyle->SetOptStat(0);
	EX->Draw("hsbeta:cointime>>h2aa",call_acept,"COLZ");
	h2aa->GetYaxis()->SetTitle("HMS #beta");
	h2aa->GetXaxis()->SetTitle("Coincidence time / ns");
	h2aa->SetTitle("");
	c2aa->Update();
	pad2png(c2aa,N+"_hsbeta_cointime_hist.png");
	c2aa->Print(N+"_hsbeta_cointime_hist.eps");
}

//********************************
//********************************
//********************************
	// W vs Q2 2D-histogram
if(type==2) {
	cout << "c2aab - W vs Q2 phase-space plot" << endl;
	TCanvas *c2aab = new TCanvas("c2aab","c2aab",600,600);
	MConfigCanvas(c2aab,0,0);
	TH2F *h2aab = new TH2F("h2aab","h2aab",150,ranges.Q2[0],ranges.Q2[1],150,ranges.W[0],ranges.W[1]);
	//gStyle->SetNumberContours(250);
	//h2aa->GetZaxis()->SetRangeUser(0.0, 200.0);
	gStyle->SetOptStat(0);
	EX->Draw("W:Q2>>h2aab",call2,"COLZ");
	//T->Draw("Q2:W>>h2aab",SIMCcall1,"COLZ"); // <- SIMC, only to check
	cout << "  + Number of events: " << EX->GetSelectedRows() << "/" << EX->GetEntries() << endl;
	if(EX->GetSelectedRows()==0) {
		cout << "Stopping execution! Number of elements selected with this cut equals zero. Nothing to do." << endl;
		return 1;
	}
	Int_t numevsTotal[2];
	numevsTotal[0] = EX->GetSelectedRows();
	numevsTotal[1] = EX->GetEntries();
	h2aab->SetTitle(Form("run %d-%d;Q^{2} (GeV^{2});W (GeV)",runN,runNfinal));
	// Drawing lines
	TLine *l2aab[4];
	for(Int_t i=0;i<4;i++) {
		l2aab[i] = new TLine(ranges.WQ2cutQ2[i],ranges.WQ2cutW[i],ranges.WQ2cutQ2[(i+1)%4],ranges.WQ2cutW[(i+1)%4]); 
		l2aab[i]->SetLineColor(1);
		l2aab[i]->SetLineWidth(2);
		l2aab[i]->Draw();
	}
	c2aab->Update();
	pad2png(c2aab,N+"_W_Q2.png");
}

//********************************
//********************************
//********************************
	// cointime vs W
if(type==2 && 1==0) { // coincidence
	cout << "c2a - cointime vs W" << endl;
	TCanvas *c2a = new TCanvas("c2a","c2a",600,600);
	MConfigCanvas(c2a,0,0);
	EX->Draw("cointime:W",call_acept,"goff");
	TGraph *gr2a_all = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	MConfigPoints(gr2a_all,4,0.4,8);
	gr2a_all->GetXaxis()->SetRangeUser(-30.0,10.0);
	gr2a_all->GetYaxis()->SetRangeUser(1.5,3.0);
	gr2a_all->SetTitle("Beam structure");
	gr2a_all->GetXaxis()->SetTitle("Coincidence Time / ns");
	gr2a_all->GetYaxis()->SetTitle("W");
	EX->Draw("cointime:W",call1,"goff");
	TGraph *gr2a = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	MConfigPoints(gr2a,3,0.4,8);
	gr2a_all->Draw("ap");
	gr2a->Draw("psame");
	c2a->Update();
	Double_t c2axlim[2], c2aylim[2];
	c2axlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c2axlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c2aylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c2aylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l2a;
	l2a.SetTextSize(0.04);
	//l2a.DrawLatex(c2axlim[0]+0.5*(c2axlim[1]-c2axlim[0]),c2aylim[0]+0.3*(c2aylim[1]-c2aylim[0]),"#color[4]{All data}");
	//l2a.DrawLatex(c2axlim[0]+0.5*(c2axlim[1]-c2axlim[0]),c2aylim[0]+0.2a5*(c2aylim[1]-c2aylim[0]),"#color[3]{Cut center structure}");
	l2a.Draw();
	c2a->Update();
	pad2png(c2a,N+"_W_cointime.png");
	c2a->Print(N+".pdf");
}

//********************************
//********************************
//********************************
	// W, Q2, t, and eps full dist
if(1==1) {
	// Preparing to save mean W, Q2, t, eps on file
	ofstream meanF;
	meanF.open(outFolder+Form("/summary_%d.csv",runN),std::ofstream::trunc);
	meanF << "# Now = " << da.AsSQLString() << endl;
	meanF << "# runN = " << runN << " .. " << runNfinal << endl;
	meanF << "# Analyzed events (within cuts) = " << numevsTotal[0] << "/" << numevsTotal[1] << endl;
	meanF << "# CUTS:" << endl << "# EXP: " << call2.GetTitle() << endl << "# SIMC: " << SIMCcall1.GetTitle() << endl;
	TString var3;
	TCanvas *c3 = new TCanvas("c3","c3",970,600);
	cout << "c3" << endl;
	MConfigCanvas(c3,0,0);
	gStyle->SetOptStat(0);
	TH1F *Wmmiss, *Wmmiss2, *Wmmiss3;
	Double_t c3xlim[2], c3ylim[2];
	TLatex l3;
	l3.SetTextSize(0.03);
	for(Int_t i=0; i<4; i++) {
		if(i==0) var3 = "Q2";
		else if(i==1) var3 = "W";
		else if(i==2) var3 = "t";
		else if(i==3) var3 = "epsilon";
		TH1F *Wmmiss = new TH1F("Wmmiss","Wmmiss",300,-0.1,3.0);
		EX->Draw(var3+">>Wmmiss",call2 * Form("%g",ExpNorm));//,chcer && ccointime1);
		cout << "  Mean("+var3+") = " << Wmmiss->GetMean() << endl;
		meanF << "# " << var3 << " \\pm s(" << var3 << ")" << endl << Wmmiss->GetMean() << "," << Wmmiss->GetRMS() << endl;
		MConfigHist(Wmmiss,var3,var3,"Yield",4,0);
		if(var3.EqualTo("t")) Wmmiss->GetXaxis()->SetRangeUser(0.2,0.6);
		else if(var3.EqualTo("W")) Wmmiss->GetXaxis()->SetRangeUser(ranges.W[0],ranges.W[1]);
		else if(var3.EqualTo("Q2")) Wmmiss->GetXaxis()->SetRangeUser(ranges.Q2[0],ranges.Q2[1]);
		else if(var3.EqualTo("epsilon")) Wmmiss->GetXaxis()->SetRangeUser(0.2,0.75);
		//Wmmiss->GetYaxis()->SetRangeUser(0.0,650.0);
//		// Fitting elastic peak
//		TF1 *f3 = new TF1("f3","[0]*exp(-0.5*((x-[1])/[2])^2)",0.8,1.1);
//		f3->SetLineColor(2);
//		f3->SetParameters(10000.0,0.938272,0.01);
//		f3->FixParameter(1,0.938272);
//		mmiss->Fit("f3","0");
		// Applying cut
		//Wmmiss2 = new TH1F("Wmmiss2","Wmmiss2",300,-0.1,3.0);
		//EX->Draw(var3+">>Wmmiss2",Form("%g",ExpNorm) * call2,"same");
		//MConfigHist(Wmmiss2,"","","",4,0);
		Wmmiss3 = new TH1F("Wmmiss3","Wmmiss3",300,-0.1,3.0);
		T->Draw(var3+">>Wmmiss3",SIMCcall1 * Form("Weight*%g",normSIMC),"same");
		MConfigHist(Wmmiss3,"","","",3,0);
		//f3->Draw("lsame");
		c3->Update();
		c3xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
		c3xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
		c3ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
		c3ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
		//l3.DrawLatex(c3xlim[0]+0.1*(c3xlim[1]-c3xlim[0]),pow(10.0,c3ylim[0]+0.8*(c3ylim[1]-c3ylim[0])),"#color[14]{No Q2, W or t cut}");
		l3.DrawLatex(c3xlim[0]+0.1*(c3xlim[1]-c3xlim[0]),c3ylim[0]+0.85*(c3ylim[1]-c3ylim[0]),"#color[4]{Experiment}");
		l3.DrawLatex(c3xlim[0]+0.1*(c3xlim[1]-c3xlim[0]),c3ylim[0]+0.8*(c3ylim[1]-c3ylim[0]),"#color[3]{SIMC}");
		//l3.DrawLatex(c3xlim[0]+0.1*(c3xlim[1]-c3xlim[0]),pow(10.0,c3ylim[0]+0.8*(c3ylim[1]-c3ylim[0])),"#color[2]{SIMC}");
		l3.Draw();
		pad2png(c3,N+"_EXP_"+var3+".png");
		if(savePDF) c3->Print(N+".pdf");
		Wmmiss->Delete();
		Wmmiss3->Delete();
	}
	meanF.close();
}
//return 0;

//********************************
//********************************
//********************************
	// W, Q2 or t, bin-by-bin dists (in phi)
if(0==1) {
//phiBinsLimits[nphibins+1]
	TString var3p = var3; //"Q2";
	TCanvas *c3p = new TCanvas("c3p","c3p",800,800);
	Int_t nlin_c3p = ((int)TMath::Sqrt(nphibins));
	Int_t ncol_c3p = ((int)(nphibins/nlin_c3p+(nphibins%nlin_c3p!=0)));
	c3p->Divide(nlin_c3p,ncol_c3p,0,0);
	TVirtualPad *p3p[nphibins];
	cout << "c3p" << endl;
	//MConfigCanvas(c3p,0,0);
	TH1F *Wmmissp[nphibins], *mmiss2p[nphibins], *mmiss3p[nphibins];
	TLatex l3p[nphibins];
	Double_t c3pxlim[2], c3pylim[2];
	for(Int_t ii=0;ii<nphibins;ii++) {
		p3p[ii] = c3p->cd(ii+1);
		TCut phi3p = Form("(180.0*phi_pq/3.1415)>=%f && (180.0*phi_pq/3.1415)<%f",phiBinsLimits[ii],phiBinsLimits[ii+1]);
		TCut phi3pSIMC = Form("(180.0*phipq/3.1415)>=%f && (180.0*phipq/3.1415)<%f",phiBinsLimits[ii],phiBinsLimits[ii+1]);
		Wmmissp[ii] = new TH1F(Form("Wmmissp%d",ii),"Wmmiss",200,-0.1,3.0);
		EX->Draw(var3p+Form(">>Wmmissp%d",ii),(call2 && phi3p) * Form("%g",ExpNorm));//,chcer && ccointime1);
		MConfigHist(Wmmissp[ii],"",var3p,"Yield",14,0);
		//Wmmissp[ii]->GetXaxis()->SetRangeUser(1.5,3.0);
		if(var3p.EqualTo("t")) Wmmissp[ii]->GetXaxis()->SetRangeUser(0.0,1.0);
		else if(var3p.EqualTo("W")) Wmmissp[ii]->GetXaxis()->SetRangeUser(2.0,3.0);
		else if(var3p.EqualTo("Q2")) Wmmissp[ii]->GetXaxis()->SetRangeUser(1.5,2.5);
		//Wmmissp[ii]->GetYaxis()->SetRangeUser(0.0,650.0);
		Wmmissp[ii]->GetXaxis()->SetLabelSize(0.08);
		Wmmissp[ii]->GetYaxis()->SetLabelSize(0.08);
		//Wmmissp[ii]->GetXaxis()->SetTitleSize(0.08);
		//Wmmissp[ii]->GetYaxis()->SetTitleSize(0.08);
		// Applying cut
		mmiss2p[ii] = new TH1F(Form("Wmmiss2p%d",ii),"Wmmiss2",200,-0.1,3.0);
		EX->Draw(var3p+Form(">>+Wmmiss2p%d",ii),(call2 && phi3p) * Form("%g",ExpNorm),"same");
		MConfigHist(mmiss2p[ii],"","","",3,0);
		// SIMC
		mmiss3p[ii] = new TH1F(Form("Wmmiss3p%d",ii),"Wmmiss3",200,-0.1,3.0);
		T->Draw(var3p+Form(">>+Wmmiss3p%d",ii),(SIMCcall1 && phi3pSIMC) * Form("Weight*%g",normSIMC),"same");
		MConfigHist(mmiss3p[ii],"","","",2,0);
		// Final config
		Wmmissp[ii]->SetDirectory(0);
		mmiss2p[ii]->SetDirectory(0);
		mmiss3p[ii]->SetDirectory(0);
		c3p->Update();
		c3pxlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
		c3pxlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
		c3pylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
		c3pylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
		l3p[ii].SetTextSize(0.1);
		//l3[ii].DrawLatex(c3pxlim[0]+0.1*(c3pxlim[1]-c3pxlim[0]),pow(10.0,c3pylim[0]+0.8*(c3pylim[1]-c3pylim[0])),"#color[14]{No Q2, W or t cut}");
		l3p[ii].DrawLatex(c3pxlim[0]+0.1*(c3pxlim[1]-c3pxlim[0]),c3pylim[0]+0.9*(c3pylim[1]-c3pylim[0]),Form("#color[1]{%.1f^{o} <= #phi < %.1f^{o}}",phiBinsLimits[ii],phiBinsLimits[ii+1]));
		if(ii==0) {
			l3p[ii].DrawLatex(c3pxlim[0]+0.1*(c3pxlim[1]-c3pxlim[0]),c3pylim[0]+0.8*(c3pylim[1]-c3pylim[0]),"#color[3]{Experiment}");
			l3p[ii].DrawLatex(c3pxlim[0]+0.1*(c3pxlim[1]-c3pxlim[0]),c3pylim[0]+0.7*(c3pylim[1]-c3pylim[0]),"#color[2]{SIMC}");
			//l3p[ii].DrawLatex(c3pxlim[0]+0.1*(c3pxlim[1]-c3pxlim[0]),pow(10.0,c3pylim[0]+0.8*(c3pylim[1]-c3pylim[0])),"#color[2]{SIMC}");
		}
		if(ii==(ncol_c3p-1)) l3p[ii].DrawLatex(c3pxlim[0]+0.8*(c3pxlim[1]-c3pxlim[0]),c3pylim[0]+0.8*(c3pylim[1]-c3pylim[0]),"#color[1]{"+var3p+"}");
		l3p[ii].Draw();
	}
	pad2png(c3p,N+"_EXP_"+var3+"_BinByBin.png");
	if(savePDF) c3p->Print(N+".pdf");
}
//return 0;
//********************************
//********************************
//********************************
	// missmass only cut events
if(type==2) {
	TCanvas *c3a = new TCanvas("c3a","c3a",600,600);
	MConfigCanvas(c3a,0,0);
	c3a->SetLeftMargin(0.13);
	TH1F *mmissa = new TH1F("mmissa","mmissa",1200,-0.5,2.0);
	TCut c3a_cutEXPall = cWQ2_k && cssdelta && cssshsum && cscer && cxptar_k && cyptar_k && chcer_k && chcal_k && chaero_k && ccointime2 && tCut && extraCut;
	EX->Draw("missmass>>+mmissa",c3a_cutEXPall);
	MConfigHist(mmissa,"","Missing mass (GeV/c^{2})","Counts",12,0);
	mmissa->GetYaxis()->SetTitleOffset(1.4);
	//mmissa->GetXaxis()->SetRangeUser(0.8,1.25);
	mmissa->GetXaxis()->SetRangeUser(1.05,1.25);
	//mmissa->GetYaxis()->SetRangeUser(0.0,7500.0);
	TH1F *mmisscut = new TH1F("mmisscut","mmisscut",1200,-0.5,2.0);
	EX->Draw("missmass>>+mmisscut",call2,"same");
	MConfigHist(mmisscut,"","","",2);
	// Sigma (just ilustration for now)
	TH1F *mmissSig = new TH1F("mmissSig","mmissSig",1200,-0.5,2.0);
	TCut c3a_cutEXPsig = call2 && "missmass>1.1825 && missmass<1.2025";
	EX->Draw("missmass>>+mmissSig",c3a_cutEXPsig,"same");
	MConfigHist(mmissSig,"","","",3);
	mmissa->SetDirectory(0);
	//mmisscut->SetDirectory(0);
	c3a->Update();
	Double_t c3axlim[2], c3aylim[3];
	c3axlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c3axlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c3aylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c3aylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l3a;
	l3a.SetTextSize(0.04);
	//l3a.DrawLatex(c3axlim[0]+0.1*(c3axlim[1]-c3axlim[0]),c3aylim[0]+0.8*(c3aylim[1]-c3aylim[0]),"#color[3]{Kaons?}");
	//l3a.Draw();
	pad2png(c3a,N+"_missmass2.png");
	c3a->Print(N+"_missmass2.eps");
}
//return 0;

//********************************
//********************************
//********************************
	// missmass vs cointime
if(type==2 && 1==1) {
	TCanvas *c4 = new TCanvas("c4","c4",600,600);
	MConfigCanvas(c4,0,0);
	EX->Draw("missmass:cointime",call_acept,"goff");
	TGraph *gr4_all = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	EX->Draw("missmass:cointime",call1,"goff");
	TGraph *gr4a = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	EX->Draw("missmass:cointime",call1err,"goff");
	TGraph *gr4err = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	gr4_all->GetXaxis()->SetRangeUser(-0.5,1.5);
	//gr4->GetYaxis()->SetRangeUser();
	MConfigPoints(gr4_all,4,0.5,8);
	MConfigPoints(gr4a,3,0.5,8);
	MConfigPoints(gr4err,2,0.5,8);
	gr4_all->GetXaxis()->SetTitle("Missing mass / GeV");
	gr4_all->GetYaxis()->SetTitle("Coincidence time / ns");
	gr4_all->SetTitle("Missing mass cut");
	gr4_all->Draw("ap");
	gr4a->Draw("psame");
	gr4err->Draw("psame");
	c4->Update();
	Double_t c4xlim[2], c4ylim[2];
	c4xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c4xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c4ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c4ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l4;
	l4.SetTextSize(0.04);
	l4.DrawLatex(c4xlim[0]+0.1*(c4xlim[1]-c4xlim[0]),c4ylim[0]+0.8*(c4ylim[1]-c4ylim[0]),"#color[4]{All data}");
	l4.DrawLatex(c4xlim[0]+0.1*(c4xlim[1]-c4xlim[0]),c4ylim[0]+0.75*(c4ylim[1]-c4ylim[0]),"#color[3]{Cut missmass and cointime}");
	l4.DrawLatex(c4xlim[0]+0.1*(c4xlim[1]-c4xlim[0]),c4ylim[0]+0.7*(c4ylim[1]-c4ylim[0]),"#color[2]{For background subtraction}");
	l4.Draw();
	// making line around cut
	TLine *li4[4];
	li4[0] = new TLine(cmissmass1_lim[0],coin2_lim[0],cmissmass1_lim[1],coin2_lim[0]);
	li4[1] = new TLine(cmissmass1_lim[1],coin2_lim[0],cmissmass1_lim[1],coin2_lim[1]);
	li4[2] = new TLine(cmissmass1_lim[1],coin2_lim[1],cmissmass1_lim[0],coin2_lim[1]);
	li4[3] = new TLine(cmissmass1_lim[0],coin2_lim[1],cmissmass1_lim[0],coin2_lim[0]);
	for(i=0;i<4;i++) {
		li4[i]->SetLineColor(3);
		li4[i]->Draw();
	}
	TLine *li4e[8];
	Double_t eps8 = 0.005;
	li4e[0] = new TLine(cmissmass1err_center[0]-cmissmass1err_width[0],coin2_lim[0],cmissmass1err_center[0]+cmissmass1err_width[0]-eps8,coin2_lim[0]);
	li4e[1] = new TLine(cmissmass1err_center[0]+cmissmass1err_width[0]-eps8,coin2_lim[0],cmissmass1err_center[0]+cmissmass1err_width[0]-eps8,coin2_lim[1]);
	li4e[2] = new TLine(cmissmass1err_center[0]+cmissmass1err_width[0]-eps8,coin2_lim[1],cmissmass1err_center[0]-cmissmass1err_width[0],coin2_lim[1]);
	li4e[3] = new TLine(cmissmass1err_center[0]-cmissmass1err_width[0],coin2_lim[1],cmissmass1err_center[0]-cmissmass1err_width[0],coin2_lim[0]);
	li4e[4] = new TLine(cmissmass1err_center[1]-cmissmass1err_width[1]+eps8,coin2_lim[0],cmissmass1err_center[1]+cmissmass1err_width[1],coin2_lim[0]);
	li4e[5] = new TLine(cmissmass1err_center[1]+cmissmass1err_width[1],coin2_lim[0],cmissmass1err_center[1]+cmissmass1err_width[1],coin2_lim[1]);
	li4e[6] = new TLine(cmissmass1err_center[1]+cmissmass1err_width[1],coin2_lim[1],cmissmass1err_center[1]-cmissmass1err_width[1]+eps8,coin2_lim[1]);
	li4e[7] = new TLine(cmissmass1err_center[1]-cmissmass1err_width[1]+eps8,coin2_lim[1],cmissmass1err_center[1]-cmissmass1err_width[1]+eps8,coin2_lim[0]);
	for(i=0;i<8;i++) {
		li4e[i]->SetLineColor(2);
		li4e[i]->Draw();
	}
	pad2png(c4,N+"_missmass_cointime.png");
	c4->Print(N+"_missmass_cointime.eps");
}
//return 0;

//********************************
//********************************
//********************************
	// beta-betap vs cointime
if(1==2) { // I dont want this plot now
	TCanvas *c4b = new TCanvas("c4b","c4b",600,600);
	MConfigCanvas(c4b,0,0);
	EX->Draw("(hsbeta-hsbeta_p):cointime",call_acept,"goff");
	TGraph *gr4b_all = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	EX->Draw("(hsbeta-hsbeta_p):cointime",call1,"goff");
	TGraph *gr4b = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	gr4b_all->GetXaxis()->SetRangeUser(-10.0,10.0);
	//gr4->GetYaxis()->SetRangeUser();
	MConfigPoints(gr4b_all,4,0.5,8);
	MConfigPoints(gr4b,3,0.5,8);
	gr4b_all->GetXaxis()->SetTitle("Missing mass / GeV");
	gr4b_all->GetYaxis()->SetTitle("(#beta - #beta_{p})");
	gr4b_all->SetTitle("");
	gr4b_all->Draw("ap");
	gr4b->Draw("psame");
	c4->Update();
	Double_t c4bxlim[2], c4bylim[2];
	c4bxlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c4bxlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c4bylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c4bylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l4b;
	l4b.SetTextSize(0.04);
	l4b.DrawLatex(c4bxlim[0]+0.1*(c4bxlim[1]-c4bxlim[0]),c4bylim[0]+0.8*(c4bylim[1]-c4bylim[0]),"#color[4]{All data}");
	l4b.DrawLatex(c4bxlim[0]+0.1*(c4bxlim[1]-c4bxlim[0]),c4bylim[0]+0.75*(c4bylim[1]-c4bylim[0]),"#color[3]{Cut missmass and cointime}");
	l4b.DrawLatex(c4bxlim[0]+0.1*(c4bxlim[1]-c4bxlim[0]),c4bylim[0]+0.7*(c4bylim[1]-c4bylim[0]),"#color[2]{For background subtraction}");
	l4b.Draw();
	pad2png(c4b,N+"_betabetap_cointime.png");
	c4b->Print(N+"_betabetap_cointime.eps");
}
//return 0;

/////////////////////////////
/////////////////////////////
/////////////////////////////
/////////////////////////////
/////////////////////////////
//    FITTING BACKGROUND   //
// (cointime and phi bins) //
/////////////////////////////
/////////////////////////////
/////////////////////////////
/////////////////////////////
/////////////////////////////

//********************************
//********************************
//********************************
	// missmass vs cointime 2D-histogram
if(type==2) {
	cname = "c4c";
	Int_t nboxes = 2;
	Double_t xc[] = {coin2_center,coin2_center};
	Double_t xw[] = {coin2_width,coin2_width};
	Double_t yc[] = {cmissmass1_center,1.18937};
	Double_t yw[] = {cmissmass1_width,cmissmass1_width};
	Double_t range[] = {-4.0,0.8,4.0,1.4};//{-4.0,0.5,4.0,1.5}; //{-1.0,1.06,2.0,1.16};
	TH2F *hh1a;
	MDraw2DMM(cname, EX, call_acept, nboxes, xc, xw, yc, yw, range, &hh1a, 150,170.0);
	gPad->Print(N+"_missmass_cointime_hist.png");
}
//return 0;

//********************************
//********************************
//********************************
	// Whole pion peak with Lambda missing mass
Double_t cx_1 = 1.2;
Double_t cw_1 = 0.4;
if(type==2 && analyzeBack) {
	cname = "c4back1";
	Int_t nboxes = 1;
	Double_t xc[] = {cx_1};
	Double_t xw[] = {cw_1};
	Double_t yc[] = {cmissmass1_center};
	Double_t yw[] = {cmissmass1_width*2.5};
	Double_t range[] = {-1.0,1.06,2.0,1.16};
	TH2F *hh1a;
	MDraw2DMM(cname, EX, call_acept, nboxes, xc, xw, yc, yw, range, &hh1a, 100, 0.0);
	gPad->Print(N+"_"+cname+"_missmass_cointime_hist.png");
}

//********************************
//********************************
//********************************

/////////////////////////////
// PREPARING PHI BIN HERE: //
////////////////////////////
Double_t Phicut[2];
Double_t NumSubtractedBack[nphibins], NumTotalEv[nphibins], MeanPhi[nphibins];

for(k=0;k<nphibins;k++) {
	NumSubtractedBack[k]=0.0;
	NumTotalEv[k]=0.0;
	MeanPhi[k]=0.0;
}

///////////////////////////////
// DEFINE COINTIME BIN HERE: //
///////////////////////////////
// if you change the here, make sure ccointime2 matches is still consistent. Otherwise the analyzed background and the total yield will be different
// 20140911: using ccointime2: [-0.6,0.9]
const int nbincointime = 1;//8;
Double_t cx_2[nbincointime] = {coin2_center};//{-0.475,-0.225,0.025,0.275,0.525,0.775,1.025,1.275}; // center of cointime bins (ns)
// IMPORTANT:
Double_t coinCut = 0.775+0.125; //ns - will calculate kaon yield below this num
Double_t cwx_2[nbincointime] = {coin2_width};//{0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125}; // width of cointime bins (ns)
Double_t cy_2[nbincointime], cwy_2[nbincointime];
for(i=0;i<nbincointime;i++) {
	cy_2[i] = cmissmass1_center;
	cwy_2[i] = cmissmass1_width*3.0;
	//cout << "  cx_2[" << i << "] = " << cx_2[i] << " \tcwx_2[" << i << "] = " << cwx_2[i] << endl;
}	

//********************************
//********************************
//********************************
// Read background if it is not to be fitted
if(!analyzeBack || type!=2) { // read background from filename: <num>+"_backPhi.csv"
	// Reading fitted background from file
	if(gSystem->GetPathInfo(Form("./fittedParameters/%d-%d_backFit.csv",runN,runNfinal),id,size,flags,mt) != 0) {// file does not exist
		cout << "<Marco> Error: could NOT read file with background: " << Form("fittedParameters/%d-%d",runN,runNfinal) << "_backFit.csv" << endl << endl;
		return 1;
	}
	ifstream fBackR(Form("./fittedParameters/%d-%d_backFit.csv",runN,runNfinal));//, std::ios_base::in);
	Double_t temp2, temp2a;
	fBackR >> temp2;
	for(k=0;k<((int)(temp2));k++) {
	    	fBackR >> temp2a >> temp2a;
		MeanPhi[k] = temp2a;
	    	fBackR >> temp2a;
		NumTotalEv[k] = temp2a;
	    	fBackR >> temp2a;
		NumSubtractedBack[k] = temp2a;
	}
	fBackR.close();
}

//********************************
//********************************
//********************************

// Fit background
else { 
for(k=0;k<nphibins;k++) {
	Phicut[0] =  phiBinsLimits[k];
	Phicut[1] =  phiBinsLimits[k+1];
	//Double_t Phicut[] = {180.0,210.0}; // 0,120,150,180,210,240,360
	TCut phi4back = Form("(180.0*phi_pq/3.1415)>=%f && (180.0*phi_pq/3.1415)<%f",Phicut[0],Phicut[1]);
	//... continue all the analyses below until end of background analysis with this loop opened ...


//********************************
//********************************
//********************************

// First fit of whole peak at Lambda missing mass
if(type==2 && analyzeBack) {
	if(saveRes) {
		Int_t fileExists = gSystem->GetPathInfo(filePhiRes,id,size,flags,mt);
		ofstream fPhi;
		fPhi.open(filePhiRes,std::ofstream::out | std::ofstream::app);
		if(fileExists != 0) // file doesn't exist, so print header
			fPhi << "# backmodels: 1-twoLines+Gaus, 2-pol3+Gaus, 3-pol3+Lorentz" << endl
			     << "# runN,Phicut[0],Phicut[1],meanPhi,centerMM,widthMM,centerCointimeBin,backmodel,KaonPeak,subtractedBack" << endl;
		fPhi << "# Now = " << da.AsSQLString() << endl;
	}
	TCanvas *c4back2 = new TCanvas("c4back2","c4back2",600,600);
	MConfigCanvas(c4back2,0,0);
	TH1F *m4back2all = new TH1F("m4back2all","m4back2all",300,0.96,1.26);
	EX->Draw("missmass>>+m4back2all",call_acept+phi4back+Form("abs(cointime-%f)<%f",cx_1,cw_1));
	MConfigHist(m4back2all,N+Form(" - %.1f<=phi<%.1f",Phicut[0],Phicut[1]),"Missing mass / GeV","counts",16,0);
	//m4back2all->GetXaxis()->SetRangeUser(0.25,1.5);
	//m4back2all->GetYaxis()->SetRangeUser(0.0,7500.0);
	m4back2all->SetDirectory(0);
	// Fits
	Double_t div4b2[] = {1.05,1.17,1.27};
	TF1 *f4back2a = new TF1("f4back2a","pol1",div4b2[0],div4b2[1]);
	TF1 *f4back2b = new TF1("f4back2b","pol1",div4b2[1],div4b2[2]);
	TF1 *f4back2c = new TF1("f4back2c","pol3",div4b2[0],div4b2[2]);
	f4back2a->SetLineColor(2);
	f4back2b->SetLineColor(2);
	f4back2c->SetLineColor(3);
	m4back2all->Fit("f4back2a","R");
	m4back2all->Fit("f4back2b","R+");
	m4back2all->Fit("f4back2c","R+");
	// Showing future cut
	TH1F *m4back2 = new TH1F("m4back2","m4back2",300,0.96,1.26);
	EX->Draw("missmass>>+m4back2",call_acept+phi4back+Form("abs(cointime-%f)<%f && abs(missmass-%f)<%f",cx_1,cw_1,cmissmass1_center,cmissmass1_width),"same");
	MConfigHist(m4back2,"","","",4,0);
	m4back2->SetDirectory(0);
	c4back2->Update();
	// Writing on plot
	Double_t c4b2xlim[2], c4b2ylim[3];
	c4b2xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c4b2xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c4b2ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c4b2ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l4b2;
	l4b2.SetTextSize(0.04);
	//l4b2.DrawLatex(c3axlim[0]+0.1*(c3axlim[1]-c3axlim[0]),c3aylim[0]+0.8*(c3aylim[1]-c3aylim[0]),"#color[3]{Kaons?}");
	gPad->Print(Nfits+Form("_c4back2_missmass_cointime_hist_phi_%.0f_%.0f.png",Phicut[0],Phicut[1]));
}

//********************************
//********************************
//********************************
	// Slicing cointime
if(type==2 && analyzeBack) {
	cout << "c4back3" << endl;
	cname = "c4back3";
	nboxes = 1;//8; // <- NUMBER OF SLICES
	Double_t range[] = {-1.0,1.06,2.0,1.16};
	TH2F *hh4b3;
	MDraw2DMM(cname, EX, call_acept, nboxes, cx_2, cwx_2, cy_2, cwy_2, range, &hh4b3, 100, 0.0);
	gPad->Print(N+"_"+cname+"_missmass_cointime_hist.png");
}

//********************************
//********************************
//********************************
	// Slicing cointime - histograms and scaling fits
if(type==2 && analyzeBack) {
	cout << "c4back4" << endl;
	Bool_t plotScaled = 1;
	TCanvas *c4back4 = new TCanvas("c4back4","c4back4",1300,600);
	//c4back4->Divide(4,2,0,0);
	const int nbin4back4 = 1;//8;
	TVirtualPad *p4back4[nbin4back4];
	TH1F *cmmKnoM4back4[nbin4back4], *cmmK4back4[nbin4back4];
	Double_t center4back4, width4back4;
	TLatex *l4back4[nbin4back4];
	Double_t meanPhi[nbin4back4];
	Int_t nii = 0; // to count number of cointime bins used to take mean phi
	// For fit with linear background
	TF1 *f4back4[nbin4back4];
	const int npar4back4 = 9;
	Double_t par4back4[npar4back4];
	par4back4[0] = 1.0; // scale factor for background
	par4back4[1] = 200.0; // Peak area
	par4back4[2] = 0.007; // peak width
	par4back4[3] = LambdaMM;//cmissmass1_center;
	par4back4[4] = f4back2a->GetParameter(0);
	par4back4[5] = f4back2a->GetParameter(1);
	par4back4[6] = f4back2b->GetParameter(0);
	par4back4[7] = f4back2b->GetParameter(1);
	par4back4[8] = div4b2[1];
	Double_t fittedScale1[nbin4back4];
	Double_t fittedArea1[nbin4back4];
	Double_t fittedBack1[nbin4back4];
	Double_t fittedScaleError1[nbin4back4];
	Double_t fittedAreaError1[nbin4back4];
	Double_t fittedBackError1[nbin4back4];
	// For fit with pol3 background
	const int npar4back4pol3 = 8;
	TF1 *f4back4pol3[nbin4back4];
	Double_t par4back4pol3[npar4back4pol3];
	par4back4pol3[0] = 1.0; // scale factor for background
	par4back4pol3[1] = 200.0; // Peak area
	par4back4pol3[2] = 0.007; // Peak width
	par4back4pol3[3] = LambdaMM;//cmissmass1_center;
	par4back4pol3[4] = f4back2c->GetParameter(0);
	par4back4pol3[5] = f4back2c->GetParameter(1);
	par4back4pol3[6] = f4back2c->GetParameter(2);
	par4back4pol3[7] = f4back2c->GetParameter(3);
	Double_t fittedScale2[nbin4back4];
	Double_t fittedArea2[nbin4back4];
	Double_t fittedBack2[nbin4back4];
	Double_t fittedScaleError2[nbin4back4];
	Double_t fittedAreaError2[nbin4back4];
	Double_t fittedBackError2[nbin4back4];
	// For fit with pol3 background and Lorentz peak (Breit-Wigner) 
	const int npar4back4pol3lorentz = 8;
	TF1 *f4back4pol3lorentz[nbin4back4];
	TF1 *f4back4pol3lorentzPeak[nbin4back4];
	Double_t par4back4pol3lorentz[npar4back4pol3lorentz];
	par4back4pol3lorentz[0] = 1.0; // scale factor for background
	par4back4pol3lorentz[1] = 200.0; // Peak area
	par4back4pol3lorentz[2] = 0.001; // Peak width
	par4back4pol3lorentz[3] = LambdaMM;//cmissmass1_center;
	par4back4pol3lorentz[4] = f4back2c->GetParameter(0);
	par4back4pol3lorentz[5] = f4back2c->GetParameter(1);
	par4back4pol3lorentz[6] = f4back2c->GetParameter(2);
	par4back4pol3lorentz[7] = f4back2c->GetParameter(3);
	Double_t fittedScale3[nbin4back4];
	Double_t fittedArea3[nbin4back4];
	Double_t fittedBack3[nbin4back4];
	Double_t fittedScaleError3[nbin4back4];
	Double_t fittedAreaError3[nbin4back4];
	Double_t fittedBackError3[nbin4back4];
	// loop
	for(Int_t ii=0;ii<nbin4back4;ii++) {
		p4back4[ii] = c4back4->cd(ii+1);
		center4back4 = cx_2[ii];//coin2_center+coin2_width*(-1.0+1.0/nbin4back4+((double)ii)*2.0/nbin4back4);
		width4back4 = cwx_2[ii];//coin2_width/nbin4back4;
		//No missmass cut
		cmmKnoM4back4[ii] = new TH1F(Form("cmmKnoM4back4%d",ii),Form("cmmKnoM4back4%d",ii),100,1.0,1.4);
		EX->Draw(Form("missmass>>+cmmKnoM4back4%d",ii),(cWQ2_k && cssdelta && cssshsum && cscer && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && Form("abs(cointime-(%f))<%f",center4back4, width4back4) && phi4back && tCut && extraCut));
		cmmKnoM4back4[ii]->GetXaxis()->SetRangeUser(1.05,1.2);//(1.05,1.25);
		if(plotScaled) cmmKnoM4back4[ii]->GetYaxis()->SetRangeUser(0.0,70.0);//600.0);
		//cmmKnoM4back4[ii]->GetYaxis()->SetRangeUser(-70+f4back4pol3lorentz[ii]->Eval(1.15),f4back4pol3lorentz[ii]->Eval(1.15)+70); // for peak tail zoom
		cmmKnoM4back4[ii]->GetXaxis()->SetLabelSize(0.06);
		if(plotScaled) cmmKnoM4back4[ii]->GetYaxis()->SetLabelSize(0.06);
		else cmmKnoM4back4[ii]->GetYaxis()->SetLabelSize(0.0);
		gPad->Update();
		MConfigHist(cmmKnoM4back4[ii],"","","",14,0);
		// Get weighed value of phi_pq
		EX->Draw("(180.0*phi_pq/3.1415)",(cWQ2_k && cssdelta && cssshsum && cscer && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && Form("abs(cointime-(%f))<%f",center4back4, width4back4) && phi4back && tCut && extraCut),"goff");
		meanPhi[ii] = TMath::Mean(EX->GetSelectedRows(),EX->GetV1());
		//With missmass cut, for peak area
		cmmK4back4[ii] = new TH1F(Form("cmmK4back4%d",ii),Form("cmmK4back4%d",ii),100,1.0,1.4);
		MConfigHist(cmmK4back4[ii],"","","",6,0);
		EX->Draw(Form("missmass>>+cmmK4back4%d",ii),(cWQ2_k && cssdelta && cssshsum && cscer && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && Form("abs(cointime-(%f))<%f",center4back4, width4back4) && cmissmass1 && phi4back && tCut && extraCut),"goff");
		//cmmK4back4[ii]->SetFillColor(6);
		// FITS
		// two pol1 background fit + gaussian
		cout << endl << "  - two pol1 background fit + gaussian:" << endl;
		f4back4[ii] = new TF1(Form("f4back4%d",ii),Fit2LinesRanges,div4b2[0],div4b2[2],npar4back4);
		f4back4[ii]->SetParameters(par4back4);
		for(i=3;i<=8;i++) f4back4[ii]->FixParameter(i,par4back4[i]);
		f4back4[ii]->SetParLimits(0,0.0,4.0);
		f4back4[ii]->SetParLimits(1,0.0,1.0e3);
		f4back4[ii]->SetParLimits(2,0.004,0.008);
		cmmKnoM4back4[ii]->Fit(Form("f4back4%d",ii),"MR0");

		f4back4[ii]->Update();
		f4back4Peak[ii] = new TF1(Form("f4back4Peak%d",ii),Fit2LinesRanges,div4b2[0],div4b2[2],npar4back4);
		for(Int_t ij=0; ij<npar4back4; ij++) f4back4Peak[ii]->SetParameter(ij,f4back4[ii]->GetParameter(ij));
		f4back4Peak[ii]->SetLineColor(2);
		f4back4Peak[ii]->Draw("lsame");

		f4back4[ii]->FixParameter(1,0.0); // Peak area
		f4back4[ii]->Draw("lsame");
		fittedScale1[ii]=f4back4[ii]->GetParameter(0);
		fittedArea1[ii]=cmmK4back4[ii]->GetEntries() - f4back4[ii]->Integral(cmissmass1_center-cmissmass1_width,cmissmass1_center+cmissmass1_width)/cmmK4back4[ii]->GetBinWidth(1);
		fittedBack1[ii]=f4back4[ii]->Integral(cmissmass1_center-cmissmass1_width,cmissmass1_center+cmissmass1_width)/cmmK4back4[ii]->GetBinWidth(1);
		fittedScaleError1[ii]=f4back4[ii]->GetParError(0);
		fittedAreaError1[ii]=0.0;//fittedArea1[ii]*sqrt(pow(f4back4[ii]->GetParError(1)/f4back4[ii]->GetParameter(1),2)+pow(f4back4[ii]->GetParError(2)/f4back4[ii]->GetParameter(2),2));
		fittedBackError1[ii]=0.0;//fittedWidth1[ii]*f4back4[ii]->GetParError(2)/f4back4[ii]->GetParameter(2);
		if(saveRes) fPhi << Form("%d,%g,%g,%g,%g,%g,%g,%d,%g,%g",runN,Phicut[0],Phicut[1],meanPhi[ii],cmissmass1_center,cmissmass1_width,center4back4,1,fittedArea1[ii],fittedBack1[ii]) << endl;
		// pol3 background fit + gaussian
		cout << endl << "  - pol3 background fit + gaussian:" << endl;
		f4back4pol3[ii] = new TF1(Form("f4back4pol3%d",ii),FitPol3,div4b2[0],div4b2[2],npar4back4pol3);
		f4back4pol3[ii]->SetLineColor(3);
		f4back4pol3[ii]->SetParameters(par4back4pol3);
		for(i=3;i<=7;i++) f4back4pol3[ii]->FixParameter(i,par4back4pol3[i]);
		f4back4pol3[ii]->SetParLimits(0,0.0,4.0);
		f4back4pol3[ii]->SetParLimits(1,0.0,1.0e3);
		f4back4pol3[ii]->SetParLimits(2,0.004,0.010);
		cmmKnoM4back4[ii]->Fit(Form("f4back4pol3%d",ii),"MR0+");

		f4back4pol3[ii]->Update();
		f4back4pol3Peak[ii] = new TF1(Form("f4back4pol3Peak%d",ii),FitPol3,div4b2[0],div4b2[2],npar4back4pol3);
		for(Int_t ij=0; ij<npar4back4pol3; ij++) f4back4pol3Peak[ii]->SetParameter(ij,f4back4pol3[ii]->GetParameter(ij));
		f4back4pol3Peak[ii]->SetLineColor(3);
		f4back4pol3Peak[ii]->Draw("lsame");

		f4back4pol3[ii]->FixParameter(1,0.0); // Peak area
		f4back4pol3[ii]->Draw("lsame"); // Peak area
		fittedScale2[ii]=f4back4pol3[ii]->GetParameter(0);
		fittedArea2[ii]=cmmK4back4[ii]->GetEntries() - f4back4pol3[ii]->Integral(cmissmass1_center-cmissmass1_width,cmissmass1_center+cmissmass1_width)/cmmK4back4[ii]->GetBinWidth(1);
		fittedBack2[ii]=f4back4pol3[ii]->Integral(cmissmass1_center-cmissmass1_width,cmissmass1_center+cmissmass1_width)/cmmK4back4[ii]->GetBinWidth(1);
		fittedScaleError2[ii]=f4back4pol3[ii]->GetParError(0);
		fittedAreaError2[ii]=0.0;
		fittedBackError2[ii]=0.0;
		if(saveRes) fPhi << Form("%d,%g,%g,%g,%g,%g,%g,%d,%g,%g",runN,Phicut[0],Phicut[1],meanPhi[ii],cmissmass1_center,cmissmass1_width,center4back4,2,fittedArea2[ii],fittedBack2[ii]) << endl;
		// pol3 background fit + Lorentz
		cout << endl << "  - pol3 background fit + Lorentz:" << endl;
		f4back4pol3lorentz[ii] = new TF1(Form("f4back4pol3lorentz%d",ii),LorentzPeakPol3,div4b2[0],div4b2[2],npar4back4pol3lorentz);
		f4back4pol3lorentz[ii]->SetLineColor(4);
		f4back4pol3lorentz[ii]->SetParameters(par4back4pol3lorentz);
		for(i=3;i<=7;i++) f4back4pol3lorentz[ii]->FixParameter(i,par4back4pol3lorentz[i]);
		f4back4pol3lorentz[ii]->SetParLimits(0,0.0,4.0);
		f4back4pol3lorentz[ii]->SetParLimits(1,0.0,10.0);
		f4back4pol3lorentz[ii]->SetParLimits(2,0.0001,0.010);
		cmmKnoM4back4[ii]->Fit(Form("f4back4pol3lorentz%d",ii),"RM0+");

		f4back4pol3lorentz[ii]->Update();
		f4back4pol3lorentzPeak[ii] = new TF1(Form("f4back4pol3lorentzPeak%d",ii),LorentzPeakPol3,div4b2[0],div4b2[2],npar4back4pol3lorentz);
		for(Int_t ij=0; ij<npar4back4pol3lorentz; ij++) f4back4pol3lorentzPeak[ii]->SetParameter(ij,f4back4pol3lorentz[ii]->GetParameter(ij));
		f4back4pol3lorentzPeak[ii]->SetLineColor(4);
		f4back4pol3lorentzPeak[ii]->Draw("lsame");

		f4back4pol3lorentz[ii]->FixParameter(1,0.0); // Peak area
		f4back4pol3lorentz[ii]->Draw("lsame"); // Peak area
		fittedScale3[ii]=f4back4pol3lorentz[ii]->GetParameter(0);
		fittedArea3[ii]=cmmK4back4[ii]->GetEntries() - f4back4pol3lorentz[ii]->Integral(cmissmass1_center-cmissmass1_width,cmissmass1_center+cmissmass1_width)/cmmK4back4[ii]->GetBinWidth(1);
		fittedBack3[ii]=f4back4pol3lorentz[ii]->Integral(cmissmass1_center-cmissmass1_width,cmissmass1_center+cmissmass1_width)/cmmK4back4[ii]->GetBinWidth(1);
		fittedScaleError3[ii]=f4back4pol3lorentz[ii]->GetParError(0);
		fittedAreaError3[ii]=0.0;
		fittedBackError3[ii]=0.0;
		if(saveRes) fPhi << Form("%d,%g,%g,%g,%g,%g,%g,%d,%g,%g",runN,Phicut[0],Phicut[1],meanPhi[ii],cmissmass1_center,cmissmass1_width,center4back4,3,fittedArea3[ii],fittedBack3[ii]) << endl;
		// Hist with cut region for peak area
		cmmK4back4[ii]->Draw("same");
		// Latex legend
		l4back4[ii] = new TLatex();
		l4back4[ii]->SetTextSize(0.05);
		l4back4[ii]->DrawLatex(1.05,0.9*gPad->GetUymax(),Form("(%d) #rightarrow |cointime - (%.5f)| < %.5f",ii,center4back4,width4back4));
		l4back4[ii]->SetTextSize(0.04);
		l4back4[ii]->DrawLatex(1.15,0.82*gPad->GetUymax(),Form("#color[2]{%f <= phi < %f}",Phicut[0],Phicut[1]));
		l4back4[ii]->DrawLatex(1.15,0.74*gPad->GetUymax(),Form("#color[2]{meanPhi = %.2f}",meanPhi[ii]));
		// Making it on root
		cmmK4back4[ii]->SetDirectory(0);
		cmmKnoM4back4[ii]->SetDirectory(0);
		// Saving results for Phi dependent background fit
		// ii - cointime cut bin
		// k - phi bin
		if(center4back4<=coinCut) { // APPLYING COINTIME CUT
			NumSubtractedBack[k] += ((fittedBack1[ii]+fittedBack2[ii]+fittedBack3[ii])/3.0) / (Phicut[1]-Phicut[0]); // per unit of angle!
			NumTotalEv[k] += (cmmK4back4[ii]->GetEntries()) / (Phicut[1]-Phicut[0]); // per unit of angle
			MeanPhi[k] += (meanPhi[ii]);
			nii++;
		}
	}
	if(nii>0) {
		MeanPhi[k] /= nii; // to take the mean value of phi, for all cointime bins in this phi bin.
	}
	c4back4->Update();
	pad2png(c4back4,Nfits+Form("_c4back4_cointime_background_sliceCointime_phi_%.0f_%.0f.png",Phicut[0],Phicut[1]));
}

//********************************
//********************************
//********************************
	// Plotting results from fits
if(type==2 && analyzeBack) {
	cout << "c4back5" << endl;
	TCanvas *c4back5 = new TCanvas("c4back5","c4back5",1300,600);
	c4back5->Divide(3,1);
	TVirtualPad *p4back5[3];
	// Scaling factor background
	p4back5[0] = c4back5->cd(1);
	TGraphErrors *gr4back5a1 = new TGraphErrors(nbin4back4,cx_2,fittedScale1,0,fittedScaleError1);
	TGraphErrors *gr4back5a2 = new TGraphErrors(nbin4back4,cx_2,fittedScale2,0,fittedScaleError2);
	TGraphErrors *gr4back5a3 = new TGraphErrors(nbin4back4,cx_2,fittedScale3,0,fittedScaleError3);
	gr4back5a1->SetMarkerSize(1.0);
	gr4back5a1->SetMarkerStyle(8);
	gr4back5a1->SetMarkerColor(2);
	gr4back5a2->SetMarkerSize(1.0);
	gr4back5a2->SetMarkerStyle(8);
	gr4back5a2->SetMarkerColor(3);
	gr4back5a3->SetMarkerSize(1.0);
	gr4back5a3->SetMarkerStyle(8);
	gr4back5a3->SetMarkerColor(4);
	gStyle->SetEndErrorSize(3);
	MConfigAxis(gr4back5a1,"Background scaling factor","Coincidence time bin (ns)","Scaling factor");
	gr4back5a1->Draw("ap");
	gr4back5a2->Draw("psame");
	gr4back5a3->Draw("psame");
	c4back5->Update();
	// Latex legend
	TLatex *l4back5 = new TLatex();
	l4back5->SetTextSize(0.025);
	l4back5->DrawLatex(gPad->GetUxmin()+0.1*(gPad->GetUxmax()-gPad->GetUxmin()),0.85*gPad->GetUymax(),Form("#color[2]{%f <= phi < %f}",Phicut[0],Phicut[1]));
	// Gaussian area
	p4back5[1] = c4back5->cd(2);
	TGraphErrors *gr4back5b1 = new TGraphErrors(nbin4back4,cx_2,fittedArea1,0,fittedAreaError1);
	TGraphErrors *gr4back5b2 = new TGraphErrors(nbin4back4,cx_2,fittedArea2,0,fittedAreaError2);
	TGraphErrors *gr4back5b3 = new TGraphErrors(nbin4back4,cx_2,fittedArea3,0,fittedAreaError3);
	gr4back5b1->SetMarkerSize(1.0);
	gr4back5b1->SetMarkerStyle(8);
	gr4back5b1->SetMarkerColor(2);
	gr4back5b2->SetMarkerSize(1.0);
	gr4back5b2->SetMarkerStyle(8);
	gr4back5b2->SetMarkerColor(3);
	gr4back5b3->SetMarkerSize(1.0);
	gr4back5b3->SetMarkerStyle(8);
	gr4back5b3->SetMarkerColor(4);
	gStyle->SetEndErrorSize(3);
	MConfigAxis(gr4back5b1,"Kaon peak","Coincidence time bin (ns)","counts");
	//gr4back5b1->GetYaxis()->SetRangeUser(0.0,10.0);
	gr4back5b1->Draw("ap");
	gr4back5b2->Draw("psame");
	gr4back5b3->Draw("psame");
	TF1 *f4back5area1 = new TF1("f4back5area1","gaus",-1.5,1.5);
	TF1 *f4back5area2 = new TF1("f4back5area2","gaus",-1.5,1.5);
	TF1 *f4back5area3 = new TF1("f4back5area3","gaus",-1.5,1.5);
	f4back5area1->SetLineColor(2);
	f4back5area2->SetLineColor(3);
	f4back5area3->SetLineColor(4);
	f4back5area1->SetLineWidth(0.5);
	f4back5area2->SetLineWidth(0.5);
	f4back5area3->SetLineWidth(0.5);
	gr4back5b1->Fit("f4back5area1","R");
	gr4back5b2->Fit("f4back5area2","R+");
	gr4back5b3->Fit("f4back5area3","R+");
	// Background counts
	p4back5[2] = c4back5->cd(3);
	TGraphErrors *gr4back5c1 = new TGraphErrors(nbin4back4,cx_2,fittedBack1,0,fittedBackError1);
	TGraphErrors *gr4back5c2 = new TGraphErrors(nbin4back4,cx_2,fittedBack2,0,fittedBackError2);
	TGraphErrors *gr4back5c3 = new TGraphErrors(nbin4back4,cx_2,fittedBack3,0,fittedBackError3);
	gr4back5c1->SetMarkerSize(1.0);
	gr4back5c1->SetMarkerStyle(8);
	gr4back5c1->SetMarkerColor(2);
	gr4back5c2->SetMarkerSize(1.0);
	gr4back5c2->SetMarkerStyle(8);
	gr4back5c2->SetMarkerColor(3);
	gr4back5c3->SetMarkerSize(1.0);
	gr4back5c3->SetMarkerStyle(8);
	gr4back5c3->SetMarkerColor(4);
	gStyle->SetEndErrorSize(3);
	MConfigAxis(gr4back5c1,"Subtracted background","Coincidence time bin (ns)","Counts");
	gr4back5c1->Draw("ap");
	gr4back5c2->Draw("psame");
	gr4back5c3->Draw("psame");
	//gr4back5b->SetDirectory(0);
	pad2png(c4back5,Nfits+Form("_c4back5_cointime_background_sliceCointime_%.0f_%.0f.png",Phicut[0],Phicut[1]));
}

////////////////////////////
} // DO NOT ERASE
////////////////////////////

// Saving fitted background to file
Int_t fileExists2 = gSystem->GetPathInfo(Form("fittedParameters/%d-%d_backFit.csv",runN,runNfinal),id,size,flags,mt);
ofstream fBack;
fBack.open(Form("fittedParameters/%d-%d_backFit.csv",runN,runNfinal),std::ofstream::trunc);//,std::ofstream::out | std::ofstream::app);
fBack << nphibins << endl;
for(k=0;k<nphibins;k++) fBack << k << " " << MeanPhi[k] << " " << NumTotalEv[k] << " " << NumSubtractedBack[k] << endl;
fBack << "# First line: number of lines with data" << endl;
fBack << "# Loop to read: k MeanPhi[k] NumTotalEv[k] NumSubtractedBack[k]" << endl;
fBack << "# MeanPhi: mean value of Phi (degrees) in the analyzed Phi bin" << endl;
fBack << "# NumTotalEv: total number of events (accept., cointime and aero cut) in that Phi bin divided by the size of the bin (events per degree)" << endl;
fBack << "# NumSubtractedBack: integrated background (accept., cointime and aero cut) in that Phi bin divided by the size of the bin (events per degree)" << endl;
fBack.close();

////////////////////////////
} // end of cycle over each background phi bin
////////////////////////////

//********************************
//********************************
//********************************
// Fitting background for rebinning and future subtraction from yield - not really using this (double check)
if(type==2) {
	cout << "c4back6" << endl;
	// Accounting for background
	TCanvas *c4back6 = new TCanvas("c4back6","c4back6",600,600);
	c4back6->SetTicks(1,1);
	// Percentual of background
	Double_t back6percent[nphibins], TotalBackEv6[nphibins];
	cout << "k\tMeanPhi[k]\tNumTotalEv[k]\tNumSubtractedBack[k]" << endl;
	for(k=0;k<nphibins;k++) {
		back6percent[k] = 100.0*NumSubtractedBack[k]/NumTotalEv[k];
		TotalBackEv6[k] = NumSubtractedBack[k]*(phiBinsLimits[k+1]-phiBinsLimits[k]);
		cout << k << "\t" << MeanPhi[k] << "\t" << NumTotalEv[k] << "\t" << NumSubtractedBack[k] << endl;
	}
	/*
	// Drawing
	TGraph *grb6 = new TGraph(nphibins,MeanPhi,back6percent);//NumSubtractedBack);//TotalBackEv6);//back6percent);//NumSubtractedBack);
	grb6->SetMarkerSize(1.0);
	grb6->SetMarkerStyle(20);
	grb6->SetMarkerColor(1);
	grb6->GetXaxis()->SetLimits(0.0,360.0);
	grb6->GetYaxis()->SetRangeUser(0.0,60.0);
	MConfigAxis(grb6,"Fitting background","Weighed phi (deg)","Background area / Total number of events (%)");
	//MConfigAxis(grb6,"Fitting background","Weighed phi (deg)","Background area per unit phi (events/deg)");
	grb6->Draw("ap");
	// Fitting
	TLegend *leg4b6 = new TLegend(0.12,0.7,0.3,0.88);
	// simple gaus
	TF1 *fbackPhi6 = new TF1("fbackPhi6","gaus",0.0,360.0);
	grb6->Fit("fbackPhi6","MR+");
	leg4b6->AddEntry(fbackPhi6,"Gaus.","l");
	// gaus + const
	cout << "Fit gaus + pol1" << endl;
	TF1 *fbackPhi6a = new TF1("fbackPhi6a","[0]+[1]*x+[2]*TMath::Exp(-0.5*TMath::Power((x-[3])/[4],2.0))",0.0,360.0);
	fbackPhi6a->SetParLimits(0,0.0,20.0);
	fbackPhi6a->SetParLimits(1,0.0,0.002);
	fbackPhi6a->SetParLimits(2,0.01*fbackPhi6->GetParameter(0),20.0*fbackPhi6->GetParameter(0));
	fbackPhi6a->SetParLimits(3,0.25*fbackPhi6->GetParameter(1),2.0*fbackPhi6->GetParameter(1));
	fbackPhi6a->SetParLimits(4,0.5*fbackPhi6->GetParameter(2),1.5*fbackPhi6->GetParameter(2));
	fbackPhi6a->SetLineColor(3);
	grb6->Fit("fbackPhi6a","MR+");
	leg4b6->AddEntry(fbackPhi6a,"Gaus. + const.","l");
	// A + cos(B*x)
	TF1 *fbackPhi6c = new TF1("fbackPhi6c","[0]+[1]*cos(x*(2*TMath::Pi()/360.0))",0.0,360.0);
	fbackPhi6c->SetLineColor(6);
	grb6->Fit("fbackPhi6c","MR+");
	leg4b6->AddEntry(fbackPhi6c,"[0] + [1]*cos(x)","l");
	// Legend
	leg4b6->Draw();
	// Save image
	pad2png(c4back6,Nfits+Form("_c4back6_Fit_background_phiDependence.png"));
	*/
}

//********************************
//********************************
//********************************
	// Plot all to finalize background analysis
if(type==2) {
	const int nbinFinalPhi = nphibins;//16;
	// meanT, meanW, meanQ2
	Double_t meanT[nbinFinalPhi], meanW2[nbinFinalPhi], meanQ2[nbinFinalPhi], meanTheta[nbinFinalPhi], meanEPS[nbinFinalPhi], meanS[nbinFinalPhi];
	Double_t phi0, phi1;
	for(i=0; i<nbinFinalPhi; i++) {
		phi0 = (((double)i)/(nbinFinalPhi))*360.0;
		phi1 = phi0+360.0/nbinFinalPhi;
		TCut phi4h1a = Form("(180.0*phi_pq/3.1415)>=%f && (180.0*phi_pq/3.1415)<%f",phi0,phi1);
		EX->Draw("t:W2:Q2:th_pq",call2 && phi4h1a,"goff");
		meanT[i] = TMath::Mean(EX->GetSelectedRows(),EX->GetV1());
		meanW2[i] = TMath::Mean(EX->GetSelectedRows(),EX->GetV2());
		meanQ2[i] = TMath::Mean(EX->GetSelectedRows(),EX->GetV3());
		meanTheta[i] = TMath::Mean(EX->GetSelectedRows(),EX->GetV4());
		EX->Draw("epsilon:W2",call2 && phi4h1a,"goff");
		meanEPS[i] = TMath::Mean(EX->GetSelectedRows(),EX->GetV1());
		meanS[i] = TMath::Mean(EX->GetSelectedRows(),EX->GetV2());
	}
	// finalize back analysis
	TCanvas *c4h1 = new TCanvas("c4h1","c4h1",700,600);
	MConfigCanvas(c4h1,0,0);
	c4h1->SetRightMargin(0.15);
	c4h1->SetLeftMargin(0.15);
	TH1F *phipq1 = new TH1F("phipq1","phipq1",nbinFinalPhi,0.0,360.0);
	EX->Draw("(phi_pq*180.0/3.14159265)>>+phipq1",call2 * Form("%g",ExpNorm));
	// consistency
	cout << endl << "CONSISTENCY TEST:" << endl << "Total number of events selected by \"call2\" cut: " << EX->GetSelectedRows() << endl;
	MConfigHist(phipq1,"","Phi_pq",Form("Yield (cts/mC) * %g",1.0/eNorm),0,0);//4,0);
	c4h1->Update();
	//phipq1->GetXaxis()->SetRangeUser(0.0,360.0);
//	phipq1->GetYaxis()->SetRangeUser(0.0,0.00015);//0.025);//gPad->GetUymax());
	phipq1->GetYaxis()->SetTitleOffset(1.4);
	phipq1->GetYaxis()->SetTitleSize(0.05);
	phipq1->GetXaxis()->SetTitleSize(0.05);
	phipq1->GetYaxis()->SetLabelSize(0.05);
	phipq1->GetXaxis()->SetLabelSize(0.05);
	TH1F *phipqBack = new TH1F("phipqBack","phipqBack",nbinFinalPhi,0.0,360.0);
	EX->Draw("(phi_pq*180.0/3.14159265)>>+phipqBack",call2back * Form("%g*%g",1.0/3.0,ExpNorm),"goff"); // Background from shift /// "goff" 
	MConfigHist(phipqBack,"","","",6,0);
	// New right axis for counts
	c4h1->Update();
	TGaxis *ax4h1 = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),gPad->GetUymin()/ExpNorm,gPad->GetUymax()/ExpNorm,phipq1->GetYaxis()->GetNdivisions(),"-=L");
	ax4h1->SetTitle("Counts (experimental)");
	ax4h1->SetTitleSize(0.04);
	ax4h1->SetTitleOffset(1.9);
	ax4h1->SetLabelOffset(-0.015);
	ax4h1->SetLabelSize(0.04);
	ax4h1->Draw();
	// Checking consistency
	Double_t temp4h1 = 0.0;
	for(k=0;k<nphibins;k++) temp4h1 += NumTotalEv[k]*(phiBinsLimits[k+1]-phiBinsLimits[k]);
	cout << "Total	number of events from background fit file: " << temp4h1 << endl << endl;
	// SIMC
	TH1F *phiSIMC4 = new TH1F("phiSIMC4","phiSIMC4",nbinFinalPhi,0.0,360.0);
	T->Draw("(phipq*180.0/3.14159265)>>+phiSIMC4",SIMCcall1 * Form("Weight*%g",normSIMC),"same");
	MConfigHist(phiSIMC4,"","","",2,0);
	// DUMMY
	TH1F *phiDUM4 = new TH1F("phiDUM4","phiDUM4",nbinFinalPhi,0.0,360.0);
	DUM->Draw("(phi_pq*180.0/3.14159265)>>+phiDUM4",call2back * Form("%g",ExpNorm),"same");
	MConfigHist(phiDUM4,"","","",7,0);
	// Graph Background from fbackPhi6, scaled to the new phi binning
	Double_t TotalBackEv6weigh[nphibins];
	//Double_t nbinsTemp;
	for(k=0;k<nphibins;k++) {
		//nbinsTemp = (phiBinsLimits[k+1]-phiBinsLimits[k])/phipq1->GetBinWidth(1);
		TotalBackEv6weigh[k] = NumSubtractedBack[k]*ExpNorm*phipq1->GetBinWidth(1);
	}
	TGraph *gr4h1B = new TGraph(nphibins,MeanPhi,TotalBackEv6weigh);//back6percent);//NumSubtractedBack);
	gr4h1B->SetMarkerStyle(20);
	gr4h1B->SetMarkerSize(1.0);
	//gr4h1B->Draw("psame"); // do not plot
	// Subtracting background
	TH1F *phiSubB = new TH1F("phiSubB","phiSubB",nbinFinalPhi,0.0,360.0);
	MConfigHist(phiSubB,"","Phi_pq","Yield",3,0);
	Double_t centerPhi4h1;
	Int_t backInd4h1;
	cout << "Background subtraction: bin by bin." << endl; //using simple gaus. fit (fbackPhi6)" << endl;
	Double_t sub;
	for(i=1;i<=nbinFinalPhi;i++) {
		centerPhi4h1 = phiSubB->GetBinCenter(i);
		//for(k=0;k<=nphibins;k++) if(centerPhi4h1>=phiBinsLimits[k] && centerPhi4h1<phiBinsLimits[k+1]) backInd4h1=k;
		for(k=0;k<=nphibins;k++) if(centerPhi4h1>=phiBinsLimits[k] && centerPhi4h1<phiBinsLimits[k+1]) backInd4h1=k; //{ backInd4h1=k; cout << "lim[0]=" << phiBinsLimits[k] << "\tcenter=" << centerPhi4h1 << "\tlim[1]=" << phiBinsLimits[k+1] << endl; }
		//if(bin4h1>0.0) 
			sub = TMath::Max(0.0,phipq1->GetBinContent(i)-TMath::Max(0.0,TotalBackEv6weigh[backInd4h1])); // REMOVING NEGATIVE!!! BE VEEEERRYYY CAREFULLL HERE!!!!!!!!!
			//sub = phipq1->GetBinContent(i)*(1.0-fbackPhi6->Eval(centerPhi4h1)/100.0); // if want to use fit, mainly if number of background bins is different of number of phi_pq bins
			phiSubB->SetBinContent(i,sub); 
			cout << "i=" << i << "\tcenterPhi=" << centerPhi4h1 << "\tEXP-BACK = " << phiSubB->GetBinContent(i) << endl;
		//else phiSubB->SetBinContent(i,0.0); // removing negative bins
	}
	phiSubB->Draw("same");

	// Conclude plots
	phipq1->SetDirectory(0);
	phipqBack->SetDirectory(0);
	phiSIMC4->SetDirectory(0);
	c4h1->Update();
	Double_t c4h1xlim[2], c4h1ylim[3];
	c4h1xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c4h1xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c4h1ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c4h1ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l4h1;
	l4h1.SetTextSize(0.03);
	//l4h1.DrawLatex(c4h1xlim[0]+0.05*(c4h1xlim[1]-c4h1xlim[0]),c4h1ylim[0]+0.9*(c4h1ylim[1]-c4h1ylim[0]),"#color[4]{Data}");
	//l4h1.DrawLatex(c4h1xlim[0]+0.05*(c4h1xlim[1]-c4h1xlim[0]),c4h1ylim[0]+0.85*(c4h1ylim[1]-c4h1ylim[0]),"#color[6]{Background from shift}");
	l4h1.DrawLatex(c4h1xlim[0]+0.05*(c4h1xlim[1]-c4h1xlim[0]),c4h1ylim[0]+0.85*(c4h1ylim[1]-c4h1ylim[0]),"#color[7]{Dummy}");
	l4h1.DrawLatex(c4h1xlim[0]+0.05*(c4h1xlim[1]-c4h1xlim[0]),c4h1ylim[0]+0.8*(c4h1ylim[1]-c4h1ylim[0]),"#color[2]{SIMC}");
	l4h1.DrawLatex(c4h1xlim[0]+0.05*(c4h1xlim[1]-c4h1xlim[0]),c4h1ylim[0]+0.75*(c4h1ylim[1]-c4h1ylim[0]),"#color[3]{Data - background subtracted}");
	// Box of cut
	if(0==1) { // To draw all cuts on phi
		TBox b4h2[nphibins];
		for(i=0;i<nphibins;i++) {
			b4h2[i].SetLineWidth(2);
			b4h2[i].SetLineColor(1);
			b4h2[i].SetLineStyle(2);
			b4h2[i].SetFillStyle(0);
			b4h2[i].DrawBox(phiBinsLimits[i],0.0,phiBinsLimits[i+1],0.999*gPad->GetUymax());
		}
	}
	// Fitting Rosenbluth form
/*	TF1 *f4h2data = new TF1("f4h2data","[0]+[1]*TMath::Cos(x*TMath::DegToRad())+[2]*TMath::Cos(2.0*x*TMath::DegToRad())",0.0,360.0);
	f4h2data->SetLineColor(2);
	phiSubB->Fit("f4h2data","RM0+");
	TF1 *f4h2SIMC = new TF1("f4h2SIMC","[0]+[1]*TMath::Cos(x*TMath::DegToRad())+[2]*TMath::Cos(2.0*x*TMath::DegToRad())",0.0,360.0);
	f4h2SIMC->SetLineColor(6);
	phiSIMC4->Fit("f4h2SIMC","RM0+");
*/
	// Saving image
	pad2png(c4h1,N+"_phipq.png");
	//c4h1->Print(N+"_phipq.eps");

	// Final counts for scaling calculation
	cout << endl << "Integrals:" << endl << "SIMC:\t\t" << phiSIMC4->Integral() << endl << "EXP-Back:\t\t" << phiSubB->Integral() << endl << "(EXP-Back)/SIMC:\t" << phiSubB->Integral()/phiSIMC4->Integral() << endl << endl;
	// Output phi bin scales for sig_factorized original
	//ofstream outScale;
	//outScale.open(Form("run%d_scaleOutput.csv",runN),std::ofstream::trunc);
	//outScale << "# runN, centerPhi, meanT, meanW, meanQ2, meanP, NevtsEXP, YieldEXP, YieldSIMC, ratioEXPSIMC" << endl;
	//for(i=0;i<nbinFinalPhi;i++) {
	//	outScale << runN << "," << phiSubB->GetBinCenter(i+1) << ",";
	//	if(phiSubB->GetBinContent(i+1)>0) outScale << meanT[i] << "," << meanW[i] << "," << meanQ2[i] << "," << meanP[i] << "," << (phiSubB->GetBinContent(i+1)/ExpNorm) << "," << phiSubB->GetBinContent(i+1) << "," << phiSIMC4->GetBinContent(i+1) << "," << phiSubB->GetBinContent(i+1)/phiSIMC4->GetBinContent(i+1) << endl;
	//	else outScale << "0.0,0.0,0.0,0.0,0.0,0.0," << phiSIMC4->GetBinContent(i+1) << ",0.0" << endl;
	//	//else phiSubB->SetBinContent(i,0.0); // removing negative bins
	//}
	//outScale.close();
	// Output phi bin scales for 
	ofstream outScale;
	outScale.open(outFolder+Form("/run%d_scaleOutput.csv",runN),std::ofstream::trunc);
	outScale << "# runN, centerPhi, meanTheta, meanT, meanQ2, meanS, meanW2, meanEPS, mtar, NevtsEXP, YieldEXP, YieldSIMC, ratioEXPSIMC" << endl;
	for(i=0;i<nbinFinalPhi;i++) {
		if(phiSubB->GetBinContent(i+1)>0) {
			outScale << runN << "," << phiSubB->GetBinCenter(i+1) << "," << meanTheta[i] << ",";
			outScale << meanT[i] << "," << meanQ2[i] << "," << meanS[i] << "," << meanW2[i] << "," << meanEPS[i] << "," << LambdaMM << "," << (phiSubB->GetBinContent(i+1)/ExpNorm) << "," << phiSubB->GetBinContent(i+1) << "," << phiSIMC4->GetBinContent(i+1) << "," << phiSubB->GetBinContent(i+1)/phiSIMC4->GetBinContent(i+1) << endl;
		}
		else {
			outScale << "# " << runN << "," << phiSubB->GetBinCenter(i+1) << "," << meanTheta[i] << ",";
			outScale << "0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0," << phiSIMC4->GetBinContent(i+1) << ",0.0" << endl;
		}
		//else phiSubB->SetBinContent(i,0.0); // removing negative bins
	}
	outScale.close();
	
}

///////////////////////////////
///////////////////////////////
///////////////////////////////
// END OF BACKGROUND FITTING //
///////////////////////////////
///////////////////////////////
///////////////////////////////

//return 0;

//********************************
//********************************
//********************************
	// missmass vs cointime 2D-histogram - background
if(type==2 && 1==0) {
	TCanvas *c4cB = new TCanvas("c4cB","c4cB",1200,600);
	MConfigCanvas(c4c,0,0);
	TH2F *h4c1B = new TH2F("h4c1B","h4c1B",2*150,-4.0,12.0,150,0.5,1.5);
	//h4c1B->SetContour(80);
	gStyle->SetNumberContours(250);
	h4c1B->GetZaxis()->SetRangeUser(0.0, 200.0);
	gStyle->SetOptStat(0);
	EX->Draw("missmass:cointime>>h4c1B",call_acept,"COLZ");
	h4c1B->GetYaxis()->SetTitle("Missing mass / GeV");
	h4c1B->GetXaxis()->SetTitle("Coincidence time / ns");
	h4c1B->SetTitle("Background");
	c4cB->Update();
	//TLatex l4cB;
	//l4cB.SetTextSize(0.04);
	//l4cB.DrawLatex(c4xlim[0]+0.1*(c4xlim[1]-c4xlim[0]),c4ylim[0]+0.8*(c4ylim[1]-c4ylim[0]),"#color[4]{All data}");
	const int n4celBorig = 5;
	Double_t x4crecBorig[n4celBorig], y4crecBorig[n4celBorig];
	y4crecBorig[0] = cmissmass1_center - cmissmass1_width;
	y4crecBorig[1] = cmissmass1_center - cmissmass1_width;
	y4crecBorig[2] = cmissmass1_center + cmissmass1_width;
	y4crecBorig[3] = cmissmass1_center + cmissmass1_width;
	y4crecBorig[4] = cmissmass1_center - cmissmass1_width;
	x4crecBorig[0] = coin2_center - coin2_width; 
	x4crecBorig[1] = coin2_center + coin2_width; 
	x4crecBorig[2] = coin2_center + coin2_width; 
	x4crecBorig[3] = coin2_center - coin2_width; 
	x4crecBorig[4] = coin2_center - coin2_width; 
	TGraph *gr4crecBorig = new TGraph(n4celBorig, x4crecBorig, y4crecBorig);
	MConfigGraphLines(gr4crecBorig,1,2,1);
	gr4crecBorig->Draw("lsame");
	const int n4celB1 = 5;
	Double_t x4crecB1[n4celB1], y4crecB1[n4celB1];
	y4crecB1[0] = cmissmass1_center - cmissmass1_width;
	y4crecB1[1] = cmissmass1_center - cmissmass1_width;
	y4crecB1[2] = cmissmass1_center + cmissmass1_width;
	y4crecB1[3] = cmissmass1_center + cmissmass1_width;
	y4crecB1[4] = cmissmass1_center - cmissmass1_width;
	x4crecB1[0] = coin2_center - coin2_width + 3.0*2.0; 
	x4crecB1[1] = coin2_center + coin2_width + 3.0*2.0; 
	x4crecB1[2] = coin2_center + coin2_width + 3.0*2.0; 
	x4crecB1[3] = coin2_center - coin2_width + 3.0*2.0; 
	x4crecB1[4] = coin2_center - coin2_width + 3.0*2.0; 
	TGraph *gr4crecB1 = new TGraph(n4celB1, x4crecB1, y4crecB1);
	MConfigGraphLines(gr4crecB1,1,2,1);
	gr4crecB1->Draw("lsame");
	const int n4celB2 = 5;
	Double_t x4crecB2[n4celB2], y4crecB2[n4celB2];
	y4crecB2[0] = cmissmass1_center - cmissmass1_width;
	y4crecB2[1] = cmissmass1_center - cmissmass1_width;
	y4crecB2[2] = cmissmass1_center + cmissmass1_width;
	y4crecB2[3] = cmissmass1_center + cmissmass1_width;
	y4crecB2[4] = cmissmass1_center - cmissmass1_width;
	x4crecB2[0] = coin2_center - coin2_width + 4.0*2.0; 
	x4crecB2[1] = coin2_center + coin2_width + 4.0*2.0; 
	x4crecB2[2] = coin2_center + coin2_width + 4.0*2.0; 
	x4crecB2[3] = coin2_center - coin2_width + 4.0*2.0; 
	x4crecB2[4] = coin2_center - coin2_width + 4.0*2.0; 
	TGraph *gr4crecB2 = new TGraph(n4celB2, x4crecB2, y4crecB2);
	MConfigGraphLines(gr4crecB2,1,2,1);
	gr4crecB2->Draw("lsame");
	const int n4celB3 = 5;
	Double_t x4crecB3[n4celB3], y4crecB3[n4celB3];
	y4crecB3[0] = cmissmass1_center - cmissmass1_width;
	y4crecB3[1] = cmissmass1_center - cmissmass1_width;
	y4crecB3[2] = cmissmass1_center + cmissmass1_width;
	y4crecB3[3] = cmissmass1_center + cmissmass1_width;
	y4crecB3[4] = cmissmass1_center - cmissmass1_width;
	x4crecB3[0] = coin2_center - coin2_width + 5.0*2.0; 
	x4crecB3[1] = coin2_center + coin2_width + 5.0*2.0; 
	x4crecB3[2] = coin2_center + coin2_width + 5.0*2.0; 
	x4crecB3[3] = coin2_center - coin2_width + 5.0*2.0; 
	x4crecB3[4] = coin2_center - coin2_width + 5.0*2.0; 
	TGraph *gr4crecB3 = new TGraph(n4celB3, x4crecB3, y4crecB3);
	MConfigGraphLines(gr4crecB3,1,2,1);
	gr4crecB3->Draw("lsame");
	pad2png(c4cB,N+"_missmass_cointime_histBackground.png");
	//c4cB->Print(N+Form"_missmass_cointime_histBackground.eps");
}


//********************************
//********************************
//********************************
	// missmass vs cointime 2D-histogram ZOOM
if(type==2 && 1==0) {
	TCanvas *c4d = new TCanvas("c4d","c4d",1200,600);
	MConfigCanvas(c4d,0,0);
	TH2F *h4d1 = new TH2F("h4d1","h4d1",100,-1.0,2.0,100,1.06,1.16);
	//h4c1->SetContour(80);
	gStyle->SetNumberContours(250);
	//h4d1->GetZaxis()->SetRangeUser(0.0, 200.0);
	gStyle->SetOptStat(0);
	EX->Draw("missmass:cointime>>h4d1",call_acept,"COLZ"); //call_acept
	h4d1->GetYaxis()->SetTitle("Missing mass / GeV");
	h4d1->GetXaxis()->SetTitle("Coincidence time / ns");
	h4d1->SetTitle("");
	c4d->Update();
	// Drawing ellipse
	//TEllipse el4d;
	//el4d.SetFillStyle(0);
	//el4d.SetFillColor(0);
	//el4d.SetLineWidth(2);
	//el4d.SetLineColor(1);
	//el4d.DrawEllipse(0.25,1.115,0.45,0.0125,0,360,0);
	//el4d.Draw("same");
	// Drawing rectangle
	const int n4el = 5;
	Double_t x4rec[n4el], y4rec[n4el];
	y4rec[0] = cmissmass1_center - cmissmass1_width;
	y4rec[1] = cmissmass1_center - cmissmass1_width;
	y4rec[2] = cmissmass1_center + cmissmass1_width;
	y4rec[3] = cmissmass1_center + cmissmass1_width;
	y4rec[4] = cmissmass1_center - cmissmass1_width;
	x4rec[0] = coin2_center - coin2_width; 
	x4rec[1] = coin2_center + coin2_width; 
	x4rec[2] = coin2_center + coin2_width; 
	x4rec[3] = coin2_center - coin2_width; 
	x4rec[4] = coin2_center - coin2_width; 
	TGraph *gr4rec = new TGraph(n4el, x4rec, y4rec);
	MConfigGraphLines(gr4rec,1,2,1);
	gr4rec->Draw("lsame");
	pad2png(c4d,N+"_missmass_cointime_histZOOM.png");
	c4d->Print(N+"_missmass_cointime_histZOOM.eps");

	//TCanvas *c4d2 = new TCanvas("c4d2","c4d2",1200,600);
	//MConfigCanvas(c4d2,0,0);
	//h4d1->ProjectionX()->Draw();
	//h4d1->GetXaxis()->SetTitle("Coincidence time / ns");
	//pad2png(c4d2,N+"_missmass_cointime_histZOOMprojection.png");
	//c4d2->Print(N+"_missmass_cointime_histZOOMprojection.eps");
}

//********************************
//********************************
//********************************
	// missmass vs cointime 2D-histogram ZOOM LEGO
if(1==2 && 1==0) {
	TCanvas *c4e = new TCanvas("c4e","c4e",1200,600);
	MConfigCanvas(c4e,0,0);
	TH2F *h4e1 = new TH2F("h4e1","h4e1",75,-1.0,2.0,75,1.06,1.16);
	//h4c1->SetContour(80);
	gStyle->SetNumberContours(20);
	//h4d1->GetZaxis()->SetRangeUser(0.0, 200.0);
	gStyle->SetOptStat(0);
	EX->Draw("missmass:cointime>>h4e1",call_acept,"LEGO2Z");
	h4e1->GetYaxis()->SetTitle("Missing mass / GeV");
	h4e1->GetXaxis()->SetTitle("Coincidence time / ns");
	h4e1->SetTitle("");
	c4e->Update();
	pad2png(c4e,N+"_missmass_cointime_histZOOMlego.png");
	c4e->Print(N+"_missmass_cointime_histZOOMlego.eps");
}

//********************************
//********************************
//********************************
	// missmass vs cointime 2D-histogram - shiffted
if(type==2 && 1==0) {
	TCanvas *c4f = new TCanvas("c4f","c4f",1200,600);
	MConfigCanvas(c4f,0,0);
	TH2F *h4f1 = new TH2F("h4f1","h4f1",150,2.0,10.0,150,0.5,1.5);
	//h4f1->SetContour(80);
	gStyle->SetNumberContours(250);
	h4f1->GetZaxis()->SetRangeUser(0.0, 200.0);
	gStyle->SetOptStat(0);
	EX->Draw("missmass:cointime>>h4f1",call_acept,"COLZ");
	h4f1->GetYaxis()->SetTitle("Missing mass / GeV");
	h4f1->GetXaxis()->SetTitle("Coincidence time / ns");
	h4f1->SetTitle("");
	c4f->Update();
	TLatex l4f;
	l4f.SetTextSize(0.04);
	//l4f.DrawLatex(c4xlim[0]+0.1*(c4xlim[1]-c4xlim[0]),c4ylim[0]+0.8*(c4ylim[1]-c4ylim[0]),"#color[4]{All data}");
	pad2png(c4f,N+"_missmass_cointime_hist_shiftT.png");
	c4f->Print(N+"_missmass_cointime_hist_shiftT.eps");
}

//********************************
//********************************
//********************************
	// missmass vs cointime 2D-histogram ZOOM - shiffted
if(type==2 && 1==0) {
	TCanvas *c4g = new TCanvas("c4g","c4g",1200,600);
	MConfigCanvas(c4g,0,0);
	TH2F *h4g1 = new TH2F("h4g1","h4g1",100,5.0,8.0,100,1.06,1.16);
	//h4c1->SetContour(80);
	gStyle->SetNumberContours(250);
	//h4g1->GetZaxis()->SetRangeUser(0.0, 200.0);
	gStyle->SetOptStat(0);
	EX->Draw("missmass:cointime>>h4g1",call_acept,"COLZ");
	h4g1->GetYaxis()->SetTitle("Missing mass / GeV");
	h4g1->GetXaxis()->SetTitle("Coincidence time / ns");
	h4g1->SetTitle("");
	c4g->Update();
	// Drawing ellipse
	TEllipse el4g;
	el4g.SetFillStyle(0);
	el4g.SetFillColor(0);
	el4g.SetLineWidth(2);
	el4g.SetLineColor(1);
	el4g.DrawEllipse(6.25,1.115,0.4,0.01,0,360,0);
	el4g.Draw("same");
	pad2png(c4g,N+"_missmass_cointime_histZOOM_shiftT.png");
	c4g->Print(N+"_missmass_cointime_histZOOM_shiftT.eps");

	TCanvas *c4g2 = new TCanvas("c4g2","c4g2",1200,600);
	MConfigCanvas(c4g2,0,0);
	h4g1->ProjectionX()->Draw();
	h4g1->GetXaxis()->SetTitle("Coincidence time / ns");
	pad2png(c4g2,N+"_missmass_cointime_histZOOMprojection_shiftT.png");
	c4g2->Print(N+"_missmass_cointime_histZOOMprojection_shiftT.eps");
}


//********************************
//********************************
//********************************
	// missmass with kaon and pion mass - pion background studies
if(type==2 && 1==0) {
	TCanvas *c4i = new TCanvas("c4i","c4i",600,600);
	MConfigCanvas(c4i,0,0);
	//No missmass cut
	TH1F *mmissKnoMcut = new TH1F("mmissKnoMcut","mmissKnoMcut",1000,-0.5,2.0);
	EX->Draw("missmass>>+mmissKnoMcut",(cWQ2_k && cssshsum && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && ccointime2));
	mmissKnoMcut->GetXaxis()->SetRangeUser(1.0,1.4);//(1.05,1.25);
	//mmissKnoMcut->GetYaxis()->SetRangeUser(0.0,7500.0);
	MConfigHist(mmissKnoMcut,"Missing mass","Missing mass / GeV","counts",14,0);
	//No missmass cut - Pion Mass
	TH1F *mmissPinoMcut = new TH1F("mmissPinoMcut","mmissPinoMcut",1000,-0.5,2.0);
	EX->Draw("M_misspi>>+mmissPinoMcut",(cWQ2_k && cssshsum && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && ccointime2),"same");//"same");
	MConfigHist(mmissPinoMcut,"","","",8,0);
	//With missmass cut
	TH1F *mmissK = new TH1F("mmissK","mmissK",1000,-0.5,2.0);
	EX->Draw("missmass>>+mmissK",call2,"same");
	MConfigHist(mmissK,"","","",4,0);
	TH1F *mmissmPi = new TH1F("mmissmPi","mmissmPi",1000,-0.5,2.0);
	EX->Draw("M_misspi>>+mmissmPi",call2,"same");//"same");
	MConfigHist(mmissmPi,"","","",3,0);
	//mmissmPi->Draw("same");
	TH1F *mmissShift = new TH1F("mmissShift","mmissShift",1000,-0.5,2.0);
	EX->Draw("missmass>>+mmissShift",call2back * Form("%g",1.0/3.0),"same");
	MConfigHist(mmissShift,"","","",6,0);
	mmissK->SetDirectory(0);
	mmissmPi->SetDirectory(0);
	mmissShift->SetDirectory(0);
	mmissKnoMcut->SetDirectory(0);
	c4i->Update();
	Double_t c4ixlim[2], c4iylim[3];
	c4ixlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c4ixlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c4iylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c4iylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l4i;
	l4i.SetTextSize(0.04);
	l4i.DrawLatex(c4ixlim[0]+0.45*(c4ixlim[1]-c4ixlim[0]),c4iylim[0]+0.8*(c4iylim[1]-c4iylim[0]),"#color[4]{Kaons mass}");
	l4i.DrawLatex(c4ixlim[0]+0.45*(c4ixlim[1]-c4ixlim[0]),c4iylim[0]+0.75*(c4iylim[1]-c4iylim[0]),"#color[3]{Pion mass}");
	l4i.DrawLatex(c4ixlim[0]+0.45*(c4ixlim[1]-c4ixlim[0]),c4iylim[0]+0.7*(c4iylim[1]-c4iylim[0]),"#color[6]{From shifted cointime}");
	l4i.DrawLatex(c4ixlim[0]+0.45*(c4ixlim[1]-c4ixlim[0]),c4iylim[0]+0.65*(c4iylim[1]-c4iylim[0]),"#color[14]{K mass - no missmass cut}");
	l4i.DrawLatex(c4ixlim[0]+0.45*(c4ixlim[1]-c4ixlim[0]),c4iylim[0]+0.6*(c4iylim[1]-c4iylim[0]),"#color[8]{Pi mass - no missmass cut}");
	l4i.Draw();
	pad2png(c4i,N+"_missmass_mK_mPi.png");
	c4i->Print(N+"_missmass_mK_mPi.eps");
}

//********************************
//********************************
//********************************
	// cointime kaon mass - background studies
if(type==2 && 1==0) {
	TCanvas *c4j = new TCanvas("c4j","c4j",1200,600);
	MConfigCanvas(c4j,0,0);
	//No missmass cut
	TH1F *ctimeKnoTcut = new TH1F("ctimeKnoTcut","ctimeKnoTcut",100,-1.0,2.5);
	EX->Draw("cointime>>+ctimeKnoTcut",(cWQ2_k && cssshsum && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && cmissmass1));
	ctimeKnoTcut->GetXaxis()->SetRangeUser(-0.75,2.0);//(1.05,1.25);
	//ctimeKnoTcut->GetYaxis()->SetRangeUser(0.0,7500.0);
	MConfigHist(ctimeKnoTcut,"Coincidence Time","cointime / ns","counts",14,0);
	//With missmass cut
	TH1F *ctimeK = new TH1F("ctimeK","ctimeK",100,-1.0,2.5);
	EX->Draw("cointime>>+ctimeK",call2,"same");
	MConfigHist(ctimeK,"","","",4,0);
	TH1F *ctimeShift = new TH1F("ctimeShift","ctimeShift",100,-1.0,2.5);
	EX->Draw("cointime>>+ctimeShift",call2back * Form("%g",1.0/3.0),"same");
	MConfigHist(ctimeShift,"","","",6,0);
	ctimeK->SetDirectory(0);
	ctimeShift->SetDirectory(0);
	ctimeKnoTcut->SetDirectory(0);
	c4j->Update();
	Double_t c4jxlim[2], c4jylim[2];
	c4jxlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c4jxlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c4jylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c4jylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l4j;
	l4j.SetTextSize(0.04);
	l4j.DrawLatex(c4jxlim[0]+0.1*(c4jxlim[1]-c4jxlim[0]),c4jylim[0]+0.8*(c4jylim[1]-c4jylim[0]),"#color[4]{Kaons mass}");
	l4j.DrawLatex(c4jxlim[0]+0.1*(c4jxlim[1]-c4jxlim[0]),c4jylim[0]+0.7*(c4jylim[1]-c4jylim[0]),"#color[6]{From shifted cointime}");
	l4j.DrawLatex(c4jxlim[0]+0.1*(c4jxlim[1]-c4jxlim[0]),c4jylim[0]+0.65*(c4jylim[1]-c4jylim[0]),"#color[14]{K mass - no cointime cut}");
	l4j.Draw();
	pad2png(c4j,N+"_cointime_background.png");
}
//return 0;

//********************************
//********************************
//********************************
	// slice cointime kaon mass - background studies
if(type==2 && 1==0) {
	cout << "c4j1" << endl;
	TCanvas *c4j1 = new TCanvas("c4j1","c4j1",1200,600);
	//MConfigCanvas(c4j1,0,0);
	c4j1->Divide(3,2,0,0);
	const int nbin4j1 = 6;
	TVirtualPad *p4j1[nbin4j1];
	TH1F *ctimeKnoTcut1[nbin4j1], *ctimeK1[nbin4j1];
	Double_t center4j1, width4j1;
	TLatex *l4j1[nbin4j1];
	//TText *l4j1[nbin4j1];
	cout << "Orig: " << endl << "   Center = " << cmissmass1_center << " \tWidth = " << cmissmass1_width << endl << "Sub-bin:" << endl;
	for(Int_t ii=0;ii<nbin4j1;ii++) {
		p4j1[ii] = c4j1->cd(ii+1);
		center4j1 = cmissmass1_center+cmissmass1_width*(-1.0+1.0/nbin4j1+((double)ii)*2.0/nbin4j1);
		width4j1 = cmissmass1_width/nbin4j1;
		cout << "  ii=" << ii << " \tCenter = " << center4j1 << " \tWidth = " << width4j1 << endl;
		//No cointime cut
		ctimeKnoTcut1[ii] = new TH1F(Form("ctimeKnoTcut1%d",ii),Form("ctimeKnoTcut1%d",ii),100,-1.0,2.5);
		EX->Draw(Form("cointime>>+ctimeKnoTcut1%d",ii),(cWQ2_k && cssshsum && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && Form("abs(missmass-%f)<%f",center4j1, width4j1)));
		ctimeKnoTcut1[ii]->GetXaxis()->SetRangeUser(-0.75,2.0);//(1.05,1.25);
		ctimeKnoTcut1[ii]->GetYaxis()->SetRangeUser(0.0,75.0);
		MConfigHist(ctimeKnoTcut1[ii],"","","",14,0);
		//With cointime cut
		ctimeK1[ii] = new TH1F(Form("ctimeK1%d",ii),Form("ctimeK1%d",ii),100,-1.0,2.5);
		EX->Draw(Form("cointime>>+ctimeK1%d",ii),(cWQ2_k && cssshsum && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && Form("abs(missmass-%f)<%f",center4j1, width4j1) && ccointime2),"same");
		MConfigHist(ctimeK1[ii],"","","",4,0);
		// Latex legend
		l4j1[ii] = new TLatex(-0.1,70.0,Form("|MM - %.5f| < %.5f",center4j1,width4j1));
		l4j1[ii]->SetTextSize(0.05);
		l4j1[ii]->Draw();
		ctimeK1[ii]->SetDirectory(0);
		ctimeKnoTcut1[ii]->SetDirectory(0);
	}
	c4j1->Update();
	pad2png(c4j1,N+"_cointime_background_sliceMissmass.png");
}

//********************************
//********************************
//********************************
	// slice missing mass kaon mass - background studies
if(type==2 && 1==0) {
	cout << "c4j2" << endl;
	TCanvas *c4j2 = new TCanvas("c4j2","c4j2",1200,600);
	c4j2->Divide(3,2,0,0);
	const int nbin4j2 = 6;
	TVirtualPad *p4j2[nbin4j2];
	TH1F *cmmKnoMcut1[nbin4j2], *cmmK1[nbin4j2];
	Double_t center4j2, width4j2;
	TLatex *l4j2[nbin4j2];
	cout << "Orig: " << endl << "   Center = " << coin2_center << " \tWidth = " << coin2_width << endl << "Sub-bin:" << endl;
	for(Int_t ii=0;ii<nbin4j2;ii++) {
		p4j2[ii] = c4j2->cd(ii+1);
		center4j2 = coin2_center+coin2_width*(-1.0+1.0/nbin4j2+((double)ii)*2.0/nbin4j2);
		width4j2 = coin2_width/nbin4j2;
		cout << "  ii=" << ii << " \tCenter = " << center4j2 << " \tWidth = " << width4j2 << endl;
		//No cointime cut
		cmmKnoMcut1[ii] = new TH1F(Form("cmmKnoMcut1%d",ii),Form("cmmKnoMcut1%d",ii),100,1.0,1.4);
		EX->Draw(Form("missmass>>+cmmKnoMcut1%d",ii),(cWQ2_k && cssshsum && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && Form("abs(cointime-(%f))<%f",center4j2, width4j2)));
		cmmKnoMcut1[ii]->GetXaxis()->SetRangeUser(1.0,1.4);//(1.05,1.25);
		cmmKnoMcut1[ii]->GetYaxis()->SetRangeUser(0.0,200.0);
		MConfigHist(cmmKnoMcut1[ii],"","","",14,0);
		//With cointime cut
		cmmK1[ii] = new TH1F(Form("cmmK1%d",ii),Form("cmmK1%d",ii),100,1.0,1.4);
		EX->Draw(Form("missmass>>+cmmK1%d",ii),(cWQ2_k && cssshsum && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && Form("abs(cointime-(%f))<%f",center4j2, width4j2) && cmissmass1),"same");
		MConfigHist(cmmK1[ii],"","","",4,0);
		// Latex legend
		l4j2[ii] = new TLatex(1.1,185.0,Form("|cointime - (%.5f)| < %.5f",center4j2,width4j2));
		l4j2[ii]->SetTextSize(0.05);
		l4j2[ii]->Draw();
		cmmK1[ii]->SetDirectory(0);
		cmmKnoMcut1[ii]->SetDirectory(0);
	}
	c4j2->Update();
	pad2png(c4j2,N+"_cointime_background_sliceCointime.png");
}

//********************************
//********************************
//********************************
	// slice missing mass kaon mass - FITTING
if(type==2 && 1==0) {
	cout << "c4j3" << endl;
	TCanvas *c4j3 = new TCanvas("c4j3","c4j3",1200,600);
	c4j3->Divide(4,2,0,0);
	const int nbin4j3 = 8; // 2 extremes are outside cut 
	TVirtualPad *p4j3[nbin4j3];
	TH1F *cmmKnoMcut3[nbin4j3], *cmmK3[nbin4j3];
	Double_t center4j3, width4j3;
	TLatex *l4j3[nbin4j3];
	cout << "Orig: " << endl << "   Center = " << coin2_center << " \tWidth = " << coin2_width << endl << "Sub-bin:" << endl;
	for(Int_t ii=0;ii<nbin4j3;ii++) {
		p4j3[ii] = c4j3->cd(ii+1);
		center4j3 = coin2_center+coin2_width*(-1.0+1.0/(nbin4j3-2)+((double)(ii-1))*2.0/(nbin4j3-2));
		width4j3 = coin2_width/(nbin4j3-2);
		cout << "  ii=" << ii << " \tCenter = " << center4j3 << " \tWidth = " << width4j3 << endl;
		//No cointime cut
		cmmKnoMcut3[ii] = new TH1F(Form("cmmKnoMcut3%d",ii),Form("cmmKnoMcut3%d",ii),100,1.0,1.4);
		EX->Draw(Form("missmass>>+cmmKnoMcut3%d",ii),(cWQ2_k && cssshsum && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && Form("abs(cointime-(%f))<%f",center4j3, width4j3)));
		cmmKnoMcut3[ii]->GetXaxis()->SetRangeUser(1.0,1.4);//(1.05,1.25);
		cmmKnoMcut3[ii]->GetYaxis()->SetRangeUser(0.0,200.0);
		MConfigHist(cmmKnoMcut3[ii],"","","",14,0);
		//With cointime cut
		cmmK3[ii] = new TH1F(Form("cmmK3%d",ii),Form("cmmK3%d",ii),100,1.0,1.4);
		EX->Draw(Form("missmass>>+cmmK3%d",ii),call2,"same");//(cWQ2_k && cssshsum && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && Form("abs(cointime-(%f))<%f",center4j3, width4j3) && cmissmass1),"same");
		MConfigHist(cmmK3[ii],"","","",4,0);
		// Latex legend
		l4j3[ii] = new TLatex(1.1,185.0,Form("|cointime - (%.5f)| < %.5f",center4j3,width4j3));
		l4j3[ii]->SetTextSize(0.05);
		l4j3[ii]->Draw();
		cmmK3[ii]->SetDirectory(0);
		cmmKnoMcut3[ii]->SetDirectory(0);
	}
	c4j3->Update();
	pad2png(c4j3,N+"_cointime_background_sliceCointime_FIT.png");
}

//********************************
//********************************
//********************************

	//scatter miss mass (M-kaon vs M-Pi)
if(type==2 && 1==0) {
	cout << "c4i2" << endl;
	// EXP: hcer_npe vs calorimeter
	TCanvas *c4i2 = new TCanvas("c4i2","c4i2",600,600);
	MConfigCanvas(c4i2,0,0);
	//-> EXP.
	EX->Draw("M_misspi:missmass",call2,"goff");
	TGraph *gr4i2 = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	MConfigPoints(gr4i2,4,0.5,8);
	gr4i2->GetXaxis()->SetTitle("Missing mass (hadron: pion)");
	gr4i2->GetYaxis()->SetTitle("Missing mass (hadron: kaon)");
	gr4i2->SetTitle("Missmass calculated from Pi and K mass");
	gr4i2->Draw("ap");
	//// Other things...
	Double_t c4i2xlim[2], c4i2ylim[2];
	c4i2->Update();
	c4i2xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c4i2xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c4i2ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c4i2ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l4i2;
	l4i2.SetTextSize(0.04);
	//l4i2.DrawLatex(c4i2xlim[0]+0.7*(c4i2xlim[1]-c4i2xlim[0]),c4i2ylim[0]+0.7*(c4i2ylim[1]-c4i2ylim[0]),"#color[4]{All data}");
	pad2png(c4i2,N+"_missmass_missmassPimass.png");
	if(savePDF) c4i2->Print(N+".pdf(");
}

//********************************
//********************************
//********************************
	// missmass hist (M-kaon vs M-pi)
if(type==2 && 1==0) {
	TCanvas *c4i3 = new TCanvas("c4i3","c4i3",600,600);
	MConfigCanvas(c4i3,0,0);
	TH2F *h4i3 = new TH2F("h4i3","h4i3",40,1.13,1.19,40,1.09,1.14);
	gStyle->SetNumberContours(250);
	//h4i3->GetZaxis()->SetRangeUser(0.0, 200.0);
	gStyle->SetOptStat(0);
	EX->Draw("missmass:M_misspi>>h4i3",call2,"COLZ"); //call_acept
	h4i3->GetXaxis()->SetTitle("Missing mass (hadron: pion)");
	h4i3->GetYaxis()->SetTitle("Missing mass (hadron: kaon)");
	h4i3->SetTitle("");
	c4i3->Update();
	pad2png(c4i3,N+"_missmass_Mpi_Mk_hist.png");
}

//********************************
//********************************
//********************************
	// missmass hist (M-kaon vs M-pi) whole
if(type==2 && 1==0) {
	TCanvas *c4i3a = new TCanvas("c4i3a","c4i3a",600,600);
	MConfigCanvas(c4i3a,0,0);
	TH2F *h4i3a = new TH2F("h4i3a","h4i3a",40,1.13,1.19,40,1.09,1.14);
	gStyle->SetNumberContours(250);
	//h4i3->GetZaxis()->SetRangeUser(0.0, 200.0);
	gStyle->SetOptStat(0);
	EX->Draw("missmass:M_misspi>>h4i3a",(cWQ2_k && cssshsum && cdelta_k && cxptar_k && cyptar_k && chcer_k && chcal_k  && chaero_k && ccointime2),"COLZ"); //call_acept
	h4i3a->GetXaxis()->SetTitle("Missing mass (hadron: pion)");
	h4i3a->GetYaxis()->SetTitle("Missing mass (hadron: kaon)");
	h4i3a->SetTitle("No missmass cut");
	c4i3a->Update();
	pad2png(c4i3a,N+"_missmass_Mpi_Mk_hist_nomissmass.png");
}


//********************************
//********************************
//********************************
	// missmass vs beta
if(type==2 && 0==1) {
	TCanvas *c5 = new TCanvas("c5","c5",1200,600);
	c5->Divide(2,1);
	MConfigCanvas(c5,0,0);
	TString var = "Q2";
	EX->Draw("missmass:"+var,call_acept,"goff");
	TGraph *gr5a = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	gr5a->SetTitle("missmass vs "+var);
	gr5a->GetXaxis()->SetTitle("missmass");
	gr5a->GetYaxis()->SetTitle(var);
	EX->Draw("cointime:"+var,call1,"goff");
	TGraph *gr5b = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	gr5b->SetTitle("cointime vs "+var);
	gr5b->GetXaxis()->SetTitle("cointime");
	gr5b->GetYaxis()->SetTitle(var);
	//gr5->SetDirectory(0);
	//gr5a->GetXaxis()->SetRangeUser(-0.5,1.5);
	//gr5b->GetXaxis()->SetRangeUser(-0.5,1.5);
	MConfigPoints(gr5a,1,0.5,8);
	MConfigPoints(gr5b,1,0.5,8);
	c5->cd(1);
	gr5a->Draw("ap");
	c5->cd(2);
	gr5b->Draw("ap");
	pad2png(c5,N+"_missmass_cointime_"+var+".png");
	c5->Print(N+"_missmass_cointime_"+var+".eps");
	//gDirectory->GetList()->Add(gr5);
}

//********************************
//********************************
//********************************
if(1==1) {
	// xptar yptar
	cout << "c6" << endl;
	TCanvas *c6 = new TCanvas("c6","c6",600,600);
	MConfigCanvas(c6,0,0);
	EX->Draw("hsxptar:hsyptar",call2,"goff");
	TGraph *gr6all = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	MConfigPoints(gr6all,3,0.5,8);
	T->Draw("hsxptar:hsyptar",SIMCcall1,"goff");
	TGraph *gr6SIMC = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	MConfigPoints(gr6SIMC,2,0.5,8);
	gr6all->GetXaxis()->SetTitle("hsxptar");
	gr6all->GetYaxis()->SetTitle("hsyptar");
	gr6all->SetTitle(Form("run%d - EXP: hms xptar vs. hms yptar",runN));
	//gr6all->GetXaxis()->SetRangeUser(-20.0,30.0);
	//gr6all->GetYaxis()->SetRangeUser(-20.0,20.0);
	gr6all->Draw("ap");
	gr6SIMC->Draw("psame");
	c6->Update();
	Double_t c6xlim[2], c6ylim[3];
	c6xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c6xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c6ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c6ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l6;
	l6.SetTextSize(0.04);
	l6.DrawLatex(c6xlim[0]+0.1*(c6xlim[1]-c6xlim[0]),c6ylim[0]+0.8*(c6ylim[1]-c6ylim[0]),"#color[3]{Data}");
	l6.DrawLatex(c6xlim[0]+0.1*(c6xlim[1]-c6xlim[0]),c6ylim[0]+0.75*(c6ylim[1]-c6ylim[0]),"#color[2]{SIMC}");
	//l6.DrawLatex(c6xlim[0]+0.1*(c6xlim[1]-c6xlim[0]),c6ylim[0]+0.75*(c6ylim[1]-c6ylim[0]),"#color[6]{Background from shift in cointime}");
	pad2png(c6,N+"_hsxptar_hsyptar.png");
	if(savePDF) c6->Print(N+".pdf");
}

//********************************
//********************************
//********************************
Double_t factor = 1.0;
if(1==1) {
	// xptar
	cout << "c6a" << endl;
	TCanvas *c6a = new TCanvas("c6a","c6a",600,600);
	MConfigCanvas(c6a,0,0);
	TH1F *gr6a = new TH1F("xptarEXP6a","xptar",20,-0.1,0.1);
	gr6a->SetLineColor(3);
	gr6a->GetXaxis()->SetTitle("hsxptar");
	gr6a->GetYaxis()->SetTitle("counts");
	gr6a->SetTitle("");
	EX->Draw("hsxptar>>xptarEXP6a",call2 * Form("%g",ExpNorm));
	TH1F *gr6aSIMC = new TH1F("xptarSIMC6a","xptar",20,-0.1,0.1);
	T->Draw("hsxptar>>xptarSIMC6a",SIMCcall1 * Form("Weight*%g",factor*normSIMC),"same");
	MConfigHist(gr6aSIMC,"","","",2,0);
	//gr6a->GetXaxis()->SetRangeUser(-20.0,30.0);
	//gr6a->GetYaxis()->SetRangeUser(-20.0,20.0);
	c6a->Update();
	Double_t c6axlim[2], c6aylim[3];
	c6axlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c6axlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c6aylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c6aylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l6a;
	l6a.SetTextSize(0.04);
	l6a.DrawLatex(c6xlim[0]+0.1*(c6xlim[1]-c6xlim[0]),c6ylim[0]+0.8*(c6ylim[1]-c6ylim[0]),"#color[3]{Data}");
	l6a.DrawLatex(c6xlim[0]+0.1*(c6xlim[1]-c6xlim[0]),c6ylim[0]+0.75*(c6ylim[1]-c6ylim[0]),"#color[2]{SIMC}");
	//l6a.DrawLatex(c6xlim[0]+0.1*(c6xlim[1]-c6xlim[0]),c6ylim[0]+0.75*(c6ylim[1]-c6ylim[0]),"#color[6]{Background from shift in cointime}");
	pad2png(c6a,N+"_hsxptar.png");
	if(savePDF) c6a->Print(N+".pdf");
}

//********************************
//********************************
//********************************
if(1==1) {
	// yptar
	cout << "c6b" << endl;
	TCanvas *c6b = new TCanvas("c6b","c6b",600,600);
	MConfigCanvas(c6b,0,0);
	TH1F *gr6b = new TH1F("yptarEXP6b","yptar",20,-0.05,0.05);
	gr6b->SetLineColor(3);
	gr6b->GetXaxis()->SetTitle("hsyptar");
	gr6b->GetYaxis()->SetTitle("counts");
	gr6b->SetTitle("");
	EX->Draw("hsyptar>>yptarEXP6b",call2 * Form("%g",ExpNorm));
	TH1F *gr6bSIMC = new TH1F("yptarSIMC6b","yptar",20,-0.05,0.05);
	T->Draw("hsyptar>>yptarSIMC6b",SIMCcall1 * Form("Weight*%g",factor*normSIMC),"same");
	MConfigHist(gr6bSIMC,"","","",2,0);
	//gr6b->GetXaxis()->SetRangeUser(-20.0,30.0);
	//gr6b->GetYaxis()->SetRangeUser(-20.0,20.0);
	c6b->Update();
	Double_t c6bxlim[2], c6bylim[3];
	c6bxlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c6bxlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c6bylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c6bylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l6b;
	l6b.SetTextSize(0.04);
	l6b.DrawLatex(c6xlim[0]+0.1*(c6xlim[1]-c6xlim[0]),c6ylim[0]+0.8*(c6ylim[1]-c6ylim[0]),"#color[3]{Data}");
	l6b.DrawLatex(c6xlim[0]+0.1*(c6xlim[1]-c6xlim[0]),c6ylim[0]+0.75*(c6ylim[1]-c6ylim[0]),"#color[2]{SIMC}");
	//l6b.DrawLatex(c6xlim[0]+0.1*(c6xlim[1]-c6xlim[0]),c6ylim[0]+0.75*(c6ylim[1]-c6ylim[0]),"#color[6]{Background from shift in cointime}");
	pad2png(c6b,N+"_hsyptar.png");
	if(savePDF) c6b->Print(N+".pdf");
}
//********************************
//********************************
//********************************
if(1==1) {
	// delta
	cout << "c6c" << endl;
	TCanvas *c6c = new TCanvas("c6c","c6c",600,600);
	MConfigCanvas(c6c,0,0);
	TH1F *gr6c = new TH1F("deltaEXP6c","delta",20,-12.0,0.0);
	gr6c->SetLineColor(3);
	gr6c->GetXaxis()->SetTitle("hsdelta");
	gr6c->GetYaxis()->SetTitle("counts");
	gr6c->SetTitle("");
	EX->Draw("hsdelta>>deltaEXP6c",call2 * Form("%g",ExpNorm));
	TH1F *gr6cSIMC = new TH1F("deltaSIMC6c","delta",20,-12.0,0.0);
	T->Draw("hsdelta>>deltaSIMC6c",SIMCcall1 * Form("Weight*%g",factor*normSIMC),"same");
	MConfigHist(gr6cSIMC,"","","",2,0);
	//gr6c->GetXaxis()->SetRangeUser(-20.0,30.0);
	//gr6c->GetYaxis()->SetRangeUser(-20.0,20.0);
	c6c->Update();
	Double_t c6cxlim[2], c6cylim[3];
	c6cxlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c6cxlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c6cylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c6cylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l6c;
	l6c.SetTextSize(0.04);
	l6c.DrawLatex(c6xlim[0]+0.1*(c6xlim[1]-c6xlim[0]),c6ylim[0]+0.8*(c6ylim[1]-c6ylim[0]),"#color[3]{Data}");
	l6c.DrawLatex(c6xlim[0]+0.1*(c6xlim[1]-c6xlim[0]),c6ylim[0]+0.75*(c6ylim[1]-c6ylim[0]),"#color[2]{SIMC}");
	//l6c.DrawLatex(c6xlim[0]+0.1*(c6xlim[1]-c6xlim[0]),c6ylim[0]+0.75*(c6ylim[1]-c6ylim[0]),"#color[6]{Background from shift in cointime}");
	pad2png(c6c,N+"_hsdelta.png");
	if(savePDF) c6c->Print(N+".pdf");
}
//return 0;

//********************************
//********************************
//********************************
	// AEROGEL
if(type==2) {
	TCanvas *c7 = new TCanvas("c7","c7",600,600);
	MConfigCanvas(c7,1,0);
	aero_all = new TH1F("aero_all","aero_all",100,-5.0,50.0);
	EX->Draw("haero_su>>+aero_all",call_acept);
	MConfigHist(aero_all,"Aerogel","NPE","counts",4,1);
	aero_all->GetXaxis()->SetRangeUser(-5.0,50.0);
	aero_all->SetDirectory(0);
	aeroK = new TH1F("aeroK","aeroK",100,-5.0,50.0);
	EX->Draw("haero_su>>+aeroK",call1,"same");//,"goff");
	MConfigHist(aeroK,"Aerogel","NPE","counts",3,1);
	aeroK->SetDirectory(0);
	// Creating text
	c7->Update();
	Double_t c7xlim[2], c7ylim[3];
	c7xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c7xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c7ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c7ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l7;
	l7.SetTextSize(0.04);
	l7.DrawLatex(c7xlim[0]+0.55*(c7xlim[1]-c7xlim[0]),pow(10.0,c7ylim[0]+0.8*(c7ylim[1]-c7ylim[0])),"#color[4]{All events}");
	l7.DrawLatex(c7xlim[0]+0.5*(c7xlim[1]-c7xlim[0]),pow(10.0,c7ylim[0]+0.75*(c7ylim[1]-c7ylim[0])),"#color[3]{Cuts applied:}");
	l7.DrawLatex(c7xlim[0]+0.5*(c7xlim[1]-c7xlim[0]),pow(10.0,c7ylim[0]+0.70*(c7ylim[1]-c7ylim[0])),"#color[3]{  - missmass}");
	l7.DrawLatex(c7xlim[0]+0.5*(c7xlim[1]-c7xlim[0]),pow(10.0,c7ylim[0]+0.65*(c7ylim[1]-c7ylim[0])),"#color[3]{  - beta}");
	l7.DrawLatex(c7xlim[0]+0.5*(c7xlim[1]-c7xlim[0]),pow(10.0,c7ylim[0]+0.60*(c7ylim[1]-c7ylim[0])),"#color[3]{  - cointime}");
	//l7.Draw();
	pad2png(c7,N+"_aero.png");
	c7->Print(N+"_aero.eps");
	
	// AEROGEL 2
	TCanvas *c8 = new TCanvas("c8","c8",600,600);
	MConfigCanvas(c8,0,0);
	aero8all = new TH1F("aero8all","aero8all",300,-5.0,50.0);
	EX->Draw("haero_su>>+aero8all",call_acept);//,"goff");
	MConfigHist(aero8all,"Aerogel","NPE","counts",4,1);
	aero8all->GetXaxis()->SetRangeUser(-10.0,50.0);
	aero8all->SetDirectory(0);
	aero8err = new TH1F("aero8err","aero8err",300,-10.0,50.0);
	EX->Draw("haero_su>>+aero8err",call1,"same");
	MConfigHist(aero8err,"Aerogel","NPE","counts",2,1);
	aero8err->SetDirectory(0);
	aero8 = new TH1F("aero8","aero8",300,-10.0,50.0);
	EX->Draw("haero_su>>+aero8",call1,"same");
	MConfigHist(aero8,"Aerogel","NPE","counts",3,1);
	aero8->SetDirectory(0);
	// Integrating
	TAxis *axis8 = aero8->GetXaxis();
	Double_t lim8[2] = {-10.0, 1.5}; //limits of integration
	Int_t bmin8 = axis8->FindBin(lim8[0]);
	Int_t bmax8 = axis8->FindBin(lim8[1]);
	Double_t integral8 = aero8->Integral(bmin8,bmax8);
	integral8 -= aero8->GetBinContent(bmin8)*(lim8[0]-axis8->GetBinLowEdge(bmin8))/axis8->GetBinWidth(bmin8);
	integral8 -= aero8->GetBinContent(bmax8)*(axis8->GetBinUpEdge(bmax8)-lim8[1])/axis8->GetBinWidth(bmax8);
	TAxis *axis8err = aero8err->GetXaxis();
	Int_t bmin8err = axis8err->FindBin(lim8[0]);
	Int_t bmax8err = axis8err->FindBin(lim8[1]);
	Double_t integral8err = aero8err->Integral(bmin8err,bmax8err);
	integral8err -= aero8err->GetBinContent(bmin8err)*(lim8[0]-axis8err->GetBinLowEdge(bmin8err))/axis8err->GetBinWidth(bmin8err);
	integral8err -= aero8err->GetBinContent(bmax8err)*(axis8err->GetBinUpEdge(bmax8err)-lim8[1])/axis8err->GetBinWidth(bmax8err);
	integral8err *= (cmissmass1_width)/(cmissmass1err_width[0]+cmissmass1err_width[1]);
	cout << "Signal PED Integral = " << integral8 << endl
	     << "Background PED Integral = " << integral8err << endl;

	///////////
	// For rejection analysis
	///////////
	lim8[0]=-10.0;
	lim8[1]=50.0;
	bmin8 = axis8->FindBin(lim8[0]);
        bmax8 = axis8->FindBin(lim8[1]);
        integral8 = aero8->Integral(bmin8,bmax8);
        integral8 -= aero8->GetBinContent(bmin8)*(lim8[0]-axis8->GetBinLowEdge(bmin8))/axis8->GetBinWidth(bmin8);
        integral8 -= aero8->GetBinContent(bmax8)*(axis8->GetBinUpEdge(bmax8)-lim8[1])/axis8->GetBinWidth(bmax8);
	Double_t int8all = integral8;
	cout << "Integrating aerogel... Total: " << int8all << endl;
	//ofstream output("Elastic_Proton_NoGasCher.csv");
	//output	<< "# Runs " << runN << " - " << runNfinal << endl
	//	<< "# Applied cuts" << endl 
	//	<< "# " << chcer->GetTitle() << endl
	//	<< "# " << chbeta->GetTitle() << endl 
	//	<< "# " << ccointime1->GetTitle() << endl 
	//	<< "# " << cmissmass1->GetTitle() << endl
	//	<< "# Minimum integration limit on \"aero_su\" = " << lim8[0] << endl
	//	<< "# max_int_lim, rejection(\%), miss_id(1:x)" << endl;
	for(lim8[1]=0.0;lim8[1]<=50.0;lim8[1]=lim8[1]+1.0) {
		bmin8 = axis8->FindBin(lim8[0]);
		bmax8 = axis8->FindBin(lim8[1]);
		integral8 = aero8->Integral(bmin8,bmax8);
		integral8 -= aero8->GetBinContent(bmin8)*(lim8[0]-axis8->GetBinLowEdge(bmin8))/axis8->GetBinWidth(bmin8);
	        integral8 -= aero8->GetBinContent(bmax8)*(axis8->GetBinUpEdge(bmax8)-lim8[1])/axis8->GetBinWidth(bmax8);
		cout << Form("[%.2f, %.2f]\t",lim8[0],lim8[1]) << Form("%.2f/%.2f = %.2f\%\tMissidentification: 1:%.2f",integral8,int8all,100.0*integral8/int8all,1.0/(1.0-integral8/int8all)) << endl;
	//	if(integral8!=int8all) output << lim8[1] << "," << integral8/int8all << "," << (1.0/(1.0-integral8/int8all)) << endl;
	//	else output << lim8[1] << "," << integral8/int8all << "," << (1.0/(1.0-integral8/int8all)) << endl;
	}
	// Creating text
	c8->Update();
	Double_t c8xlim[2], c8ylim[3];
	c8xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c8xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c8ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c8ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l8;
	l8.SetTextSize(0.04);
	l8.DrawLatex(c8xlim[0]+0.55*(c8xlim[1]-c8xlim[0]),c8ylim[0]+0.8*(c8ylim[1]-c8ylim[0]),"#color[4]{All events}");
	l8.DrawLatex(c8xlim[0]+0.5*(c8xlim[1]-c8xlim[0]),c8ylim[0]+0.75*(c8ylim[1]-c8ylim[0]),"#color[3]{Cuts applied:}");
	l8.DrawLatex(c8xlim[0]+0.5*(c8xlim[1]-c8xlim[0]),c8ylim[0]+0.70*(c8ylim[1]-c8ylim[0]),"#color[3]{  - missmass}");
	l8.DrawLatex(c8xlim[0]+0.5*(c8xlim[1]-c8xlim[0]),c8ylim[0]+0.65*(c8ylim[1]-c8ylim[0]),"#color[3]{  - beta}");
	l8.DrawLatex(c8xlim[0]+0.5*(c8xlim[1]-c8xlim[0]),c8ylim[0]+0.60*(c8ylim[1]-c8ylim[0]),"#color[3]{  - cointime}");
	l8.DrawLatex(c8xlim[0]+0.33*(c8xlim[1]-c8xlim[0]),c8ylim[0]+0.75*(c8ylim[1]-c8ylim[0]),"#color[3]{Cut missmass && cointime && beta}");
	l8.DrawLatex(c8xlim[0]+0.5*(c8xlim[1]-c8xlim[0]),c8ylim[0]+0.7*(c8ylim[1]-c8ylim[0]),Form("#color[2]{PED cut = %.2f PE}",lim8[1]));
	l8.DrawLatex(c8xlim[0]+0.35*(c8xlim[1]-c8xlim[0]),c8ylim[0]+0.55*(c8ylim[1]-c8ylim[0]),Form("#color[2]{#frac{PED-PED_{back}}{Total-Total_{back}} = #frac{%.0f-%.0f}{%.0f-%.0f} = %.2f%}",integral8,integral8err,aero8->GetEntries(),aero8err->GetEntries(),100.0*(integral8-integral8err)/(aero8->GetEntries()-aero8err->GetEntries())));

	l8.Draw();
	TArrow ar8(lim8[1],c8ylim[0]+0.5*(c8ylim[1]-c8ylim[0]),lim8[1],c8ylim[0]+0.3*(c8ylim[1]-c8ylim[0]),0.07,"|>");
	ar8.SetLineColor(2);
	ar8.SetFillColor(2);
	ar8.SetLineWidth(2);
	ar8.SetAngle(20);
	//ar8.Draw();
	pad2png(c8,N+"_aero2.png");
	c8->Print(N+"_aero2.eps");
}

//********************************
//********************************
//********************************
	// SIMC Weight
if(1==0) {
	TCanvas *c9 = new TCanvas("c9","c9",600,600);
	MConfigCanvas(c9,0,0);
	TH1F *mweight = new TH1F("mweight","mweight",200,0.0,0.1);
	T->Draw("Weight>>+mweight",call_acept);
	mweight->SetDirectory(0);
	MConfigHist(mweight,titles+" - SIMC Weight","Weight","counts",4,0);
	//mweight->GetYaxis()->SetRangeUser(0.0,2000.0);
	TH1F *mweight2 = new TH1F("mweight2","mweight2",200,0.0,0.1);
	T->Draw("Weight>>+mweight2",SIMCcall1,"same");
	mweight2->SetDirectory(0);
	MConfigHist(mweight2,"","","",3,0);
	c9->Update();
	Double_t c9xlim[2], c9ylim[3];
	c9xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c9xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c9ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c9ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l9;
	l9.SetTextSize(0.04);
	//l9.DrawLatex(c9xlim[0]+0.1*(c9xlim[1]-c9xlim[0]),c9ylim[0]+0.8*(c9ylim[1]-c9ylim[0]),"#color[4]{All data}");
	//l9.DrawLatex(c9xlim[0]+0.1*(c9xlim[1]-c9xlim[0]),c9ylim[0]+0.75*(c9ylim[1]-c9ylim[0]),"#color[3]{Cut missmass and cointime}");
	//l9.Draw();
	pad2png(c9,N+"_Weight.png");
	if(savePDF) c9->Print(N+".pdf");
}

//********************************
//********************************
//********************************
	// HMS hszbeam
if(1==0) {
	TCanvas *c10 = new TCanvas("c10","c10",600,600);
	MConfigCanvas(c10,0,0);
	TH1F *mhszbeam = new TH1F("mhszbeam","mhszbeam",200,-4.0,4.0);
	EX->Draw("hszbeam>>+mhszbeam",call_acept);
	mhszbeam->SetDirectory(0);
	MConfigHist(mhszbeam,titles+" - EXP z beam","hszbeam","counts",4,0);
	TH1F *mhszbeam2 = new TH1F("mhszbeam2","mhszbeam2",200,-4.0,4.0);
	EX->Draw("hszbeam>>+mhszbeam2",call1,"same");
	mhszbeam2->SetDirectory(0);
	MConfigHist(mhszbeam2,"","","",3,0);
	c10->Update();
	Double_t c10xlim[2], c10ylim[3];
	c10xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c10xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c10ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c10ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l10;
	l10.SetTextSize(0.04);
	//l10.DrawLatex(c10xlim[0]+0.1*(c10xlim[1]-c10xlim[0]),c10ylim[0]+0.8*(c10ylim[1]-c10ylim[0]),"#color[4]{All data}");
	//l10.DrawLatex(c10xlim[0]+0.1*(c10xlim[1]-c10xlim[0]),c10ylim[0]+0.75*(c10ylim[1]-c10ylim[0]),"#color[3]{Cut missmass and cointime}");
	//l10.Draw();
	pad2png(c10,N+"_hszbeam.png");
	if(savePDF) c10->Print(N+".pdf)");
}

//********************************
//********************************
//********************************
	// hsbeta
if(1==0) {
	TCanvas *c11 = new TCanvas("c11","c11",600,600);
	MConfigCanvas(c11,0,0);
	TH1F *mhdelta = new TH1F("mhdelta","mhdelta",200,-0.5,1.5);
	EX->Draw("hsbeta>>+mhdelta",call_acept);//,chcer && ccointime1);
	mhdelta->SetDirectory(0);
	MConfigHist(mhdelta,"EXP hsbeta","hsbeta","counts",4,1);
	// Applying cut
	TH1F *mhdelta2 = new TH1F("mhdelta2","mhdelta2",200,-0.5,1.5);
	EX->Draw("hsbeta>>+mhdelta2",call1,"same");
	mhdelta2->SetDirectory(0);
	MConfigHist(mhdelta2,"","","",3,1);
	c11->Update();
	Double_t c11xlim[2], c11ylim[3];
	c11xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c11xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c11ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c11ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l11;
	l11.SetTextSize(0.04);
	l11.DrawLatex(c11xlim[0]+0.1*(c11xlim[1]-c11xlim[0]),c11ylim[0]+0.8*(c11ylim[1]-c11ylim[0]),"#color[4]{All data}");
	l11.DrawLatex(c11xlim[0]+0.1*(c11xlim[1]-c11xlim[0]),c11ylim[0]+0.75*(c11ylim[1]-c11ylim[0]),"#color[3]{After all cuts}");
	//l3.Draw();
	pad2png(c11,N+"_EXP_hsbeta.png");
	if(savePDF) c11->Print(N+".pdf");
}

//********************************
//********************************
//********************************
	// haero_su
if(1==0) {
	TCanvas *c12 = new TCanvas("c12","c12",600,600);
	MConfigCanvas(c12,0,0);
	TH1F *mhaero_su = new TH1F("mhaero_su","mhaero_su",200,-1.5,25.0);
	EX->Draw("haero_su>>+mhaero_su",call_acept);//,chcer && ccointime1);
	mhaero_su->SetDirectory(0);
	MConfigHist(mhaero_su,"EXP haero_su","haero_su","counts",4,1);
	// Applying cut
	TH1F *mhaero_su2 = new TH1F("mhaero_su2","mhaero_su2",200,-1.5,25.0);
	EX->Draw("haero_su>>+mhaero_su2",call1,"same");
	mhaero_su2->SetDirectory(0);
	MConfigHist(mhaero_su2,"","","",3,1);
	c12->Update();
	Double_t c12xlim[2], c12ylim[3];
	c12xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c12xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c12ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c12ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l12;
	l12.SetTextSize(0.04);
	l12.DrawLatex(c12xlim[0]+0.1*(c12xlim[1]-c12xlim[0]),c12ylim[0]+0.8*(c12ylim[1]-c12ylim[0]),"#color[4]{All data}");
	l12.DrawLatex(c12xlim[0]+0.1*(c12xlim[1]-c12xlim[0]),c12ylim[0]+0.75*(c12ylim[1]-c12ylim[0]),"#color[3]{After all cuts}");
	//l3.Draw();
	pad2png(c12,N+"_EXP_haero_su.png");
	if(savePDF) c12->Print(N+".pdf");
}

//********************************
//********************************
//********************************
	// Q2
if(1==0) {
	TCanvas *c13 = new TCanvas("c13","c13",600,600);
	MConfigCanvas(c13,0,0);
	TH1F *mQ2 = new TH1F("mQ2","mQ2",300,0.0,6.0);
	EX->Draw("Q2>>+mQ2",call2);//,chcer && ccointime1);
	mQ2->SetDirectory(0);
	MConfigHist(mQ2,"EXP Q2","Q2","counts",4,0);
	// Applying cut
	TH1F *mQ22 = new TH1F("mQ22","mQ22",300,0.0,6.0);
	EX->Draw("Q2>>+mQ22",call2back * Form("%g",1.0/3.0),"same");
	mQ22->SetDirectory(0);
	MConfigHist(mQ22,"","","",6,0);
	c13->Update();
	Double_t c13xlim[2], c13ylim[3];
	c13xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c13xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c13ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c13ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l13;
	l13.SetTextSize(0.04);
	l13.DrawLatex(c13xlim[0]+0.1*(c13xlim[1]-c13xlim[0]),c13ylim[0]+0.8*(c13ylim[1]-c13ylim[0]),"#color[4]{Kaon}");
	l13.DrawLatex(c13xlim[0]+0.1*(c13xlim[1]-c13xlim[0]),c13ylim[0]+0.75*(c13ylim[1]-c13ylim[0]),"#color[6]{Background from shift cointime}");
	pad2png(c13,N+"_EXP_Q2.png");
	if(savePDF) c13->Print(N+".pdf");
}
	if(fast) {
		if(saveRes) fPhi.close();
		gROOT->ProcessLine(".! spd-say \"Marco, work completed!\"");
		return 1;
	}

//********************************
//********************************
//********************************
	// hsdedx1
	TCanvas *c14 = new TCanvas("c14","c14",600,600);
	MConfigCanvas(c14,0,0);
//	TH1F *mhsdedx1 = new TH1F("mhsdedx1","mhsdedx1",200,0.0,3.0);
	EX->Draw("hsdedx1>>+mhsdedx1",call_acept);//,chcer && ccointime1);
	mhsdedx1->SetDirectory(0);
	MConfigHist(mhsdedx1,"EXP hsdedx1","hsdedx1","counts",4,1);
	// Applying cut
//	TH1F *mhsdedx12 = new TH1F("mhsdedx12","mhsdedx12",200,0.0,3.0);
	EX->Draw("hsdedx1>>+mhsdedx12",call1,"same");
	mhsdedx12->SetDirectory(0);
	MConfigHist(mhsdedx12,"","","",3,1);
	c14->Update();
	Double_t c14xlim[2], c14ylim[3];
	c14xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c14xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c14ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c14ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l14;
	l14.SetTextSize(0.04);
	l14.DrawLatex(c14xlim[0]+0.1*(c14xlim[1]-c14xlim[0]),c14ylim[0]+0.8*(c14ylim[1]-c14ylim[0]),"#color[4]{All data}");
	l14.DrawLatex(c14xlim[0]+0.1*(c14xlim[1]-c14xlim[0]),c14ylim[0]+0.75*(c14ylim[1]-c14ylim[0]),"#color[3]{After all cuts}");
	//l3.Draw();
	pad2png(c14,N+"_EXP_hsdedx1.png");
	if(savePDF) c14->Print(N+".pdf");

//********************************
//********************************
//********************************
	// theta phi
	TCanvas *c15 = new TCanvas("c15","c15",600,600);
	MConfigCanvas(c15,0,0);
	EX->Draw("(hsphi*180.0/3.1415):(hstheta*180.0/3.1415)",call_acept,"goff");
	TGraph *gr15all = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	EX->Draw("(hsphi*180.0/3.1415):(hstheta*180.0/3.1415)",call1,"goff");
	TGraph *gr15 = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	MConfigPoints(gr15,2,0.5,8);
	MConfigPoints(gr15all,1,0.5,8);
	gr15all->GetYaxis()->SetTitle("hstheta");
	gr15all->GetXaxis()->SetTitle("hsphi");
	gr15all->SetTitle(Form("run%d - EXP: hms theta vs. hms phi",runN));
	//gr15all->GetXaxis()->SetRangeUser(-20.0,30.0);
	//gr15all->GetYaxis()->SetRangeUser(-20.0,20.0);
	gr15all->Draw("ap");
	gr15->Draw("psame");
	c15->Update();
	Double_t c15xlim[2], c15ylim[3];
	c15xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c15xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c15ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c15ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l15;
	l15.SetTextSize(0.04);
	l15.DrawLatex(c15xlim[0]+0.1*(c15xlim[1]-c15xlim[0]),c15ylim[0]+0.8*(c15ylim[1]-c15ylim[0]),"#color[1]{All data}");
	l15.DrawLatex(c15xlim[0]+0.1*(c15xlim[1]-c15xlim[0]),c15ylim[0]+0.75*(c15ylim[1]-c15ylim[0]),"#color[3]{Cut missmass and cointime}");
	pad2png(c15,N+"_EXP_hstheta_hsphi.png");
	if(savePDF) c15->Print(N+".pdf");

//********************************
//********************************
//********************************
	// W hsshsum
	TCanvas *c16 = new TCanvas("c16","c16",600,600);
	MConfigCanvas(c16,0,0);
	EX->Draw("hsshsum:W",call_acept,"goff");
	TGraph *gr16all = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	EX->Draw("hsshsum:W",call1,"goff");
	TGraph *gr16 = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	MConfigPoints(gr16,3,0.5,8);
	MConfigPoints(gr16all,4,0.5,8);
	gr16all->GetYaxis()->SetTitle("W");
	gr16all->GetXaxis()->SetTitle("hsshsum");
	gr16all->SetTitle(Form("run%d - EXP: W vs. calorimeter",runN));
	//gr16all->GetXaxis()->SetRangeUser(-20.0,30.0);
	//gr16all->GetYaxis()->SetRangeUser(-20.0,20.0);
	gr16all->Draw("ap");
	gr16->Draw("psame");
	c16->Update();
	Double_t c16xlim[2], c16ylim[3];
	c16xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c16xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c16ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c16ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l16;
	l16.SetTextSize(0.04);
	l16.DrawLatex(c16xlim[0]+0.1*(c16xlim[1]-c16xlim[0]),c16ylim[0]+0.8*(c16ylim[1]-c16ylim[0]),"#color[1]{All data}");
	l16.DrawLatex(c16xlim[0]+0.1*(c16xlim[1]-c16xlim[0]),c16ylim[0]+0.75*(c16ylim[1]-c16ylim[0]),"#color[3]{Cut}");
	pad2png(c16,N+"_EXP_W_hsshsum.png");
	if(savePDF) c16->Print(N+".pdf");

//********************************
//********************************
//********************************
	// W dedx1
	TCanvas *c17 = new TCanvas("c17","c17",600,600);
	MConfigCanvas(c17,0,0);
	EX->Draw("hsdedx1:W",call_acept,"goff");
	TGraph *gr17all = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	EX->Draw("hsdedx1:W",call1,"goff");
	TGraph *gr17 = new TGraph(EX->GetSelectedRows(),EX->GetV1(),EX->GetV2());
	MConfigPoints(gr17,2,0.5,8);
	MConfigPoints(gr17all,1,0.5,8);
	gr17all->GetYaxis()->SetTitle("W");
	gr17all->GetXaxis()->SetTitle("hsdedx1");
	gr17all->SetTitle(Form("run%d - EXP: W vs. dedx1",runN));
	//gr17all->GetXaxis()->SetRangeUser(-20.0,30.0);
	//gr17all->GetYaxis()->SetRangeUser(-20.0,20.0);
	gr17all->Draw("ap");
	gr17->Draw("psame");
	c17->Update();
	Double_t c17xlim[2], c17ylim[3];
	c17xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c17xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c17ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c17ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l17;
	l17.SetTextSize(0.04);
	l17.DrawLatex(c17xlim[0]+0.1*(c17xlim[1]-c17xlim[0]),c17ylim[0]+0.8*(c17ylim[1]-c17ylim[0]),"#color[1]{All data}");
	l17.DrawLatex(c17xlim[0]+0.1*(c17xlim[1]-c17xlim[0]),c17ylim[0]+0.75*(c17ylim[1]-c17ylim[0]),"#color[3]{Cut missmass and cointime}");
	pad2png(c17,N+"_EXP_W_hsdedx1.png");
	if(savePDF) c17->Print(N+".pdf");

//********************************
//********************************
//********************************
	// background fit of W
	TCanvas *c18 = new TCanvas("c18","c18",1200,600);
	c18->Divide(2,1);
	TPad *p18a = c18->cd(1);
	//MConfigCanvas(p18a,1,0);
	p18a->SetGrid();
	p18a->GetFrame()->SetBorderSize(10);
	p18a->SetLogy();
	TH1F *mWbackFit = new TH1F("mWbackFit","mWbackFit",200,0.0,3.0);
	EX->Draw("W>>+mWbackFit",call_acept);//,chcer && ccointime1);
	mWbackFit->SetDirectory(0);
	MConfigHist(mWbackFit,"EXP W background fit","W","counts",4,1);
	// Applying cut
	TH1F *mWbackFit2 = new TH1F("mWbackFit2","mWbackFit2",200,0.0,3.0);
	EX->Draw("W>>+mWbackFit2",call1backFit,"same");
	mWbackFit2->SetDirectory(0);
	MConfigHist(mWbackFit2,"","","",3,1);
	// Fitting background
	TF1 *f18 = new TF1("f18","exp([0]+[1]*x)",0.6,1.2);
	f18->SetLineColor(2);
	f18->SetParameters(-3.0,10.0);
	//f18->FixParameter(1,0.938272);
	mWbackFit2->Fit("f18","R");
	c18->Update();
	Double_t c18xlim[2], c18ylim[3];
	c18xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c18xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c18ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c18ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l18;
	l18.SetTextSize(0.04);
	l18.DrawLatex(c18xlim[0]+0.1*(c18xlim[1]-c18xlim[0]),c18ylim[0]+0.8*(c18ylim[1]-c18ylim[0]),"#color[4]{All data}");
	l18.DrawLatex(c18xlim[0]+0.1*(c18xlim[1]-c18xlim[0]),c18ylim[0]+0.75*(c18ylim[1]-c18ylim[0]),"#color[3]{After all cuts}");
	//l3.Draw();
	TPad *p18b = c18->cd(2);
	//MConfigCanvas(p18a,1,0);
	p18b->SetGrid();
	p18b->GetFrame()->SetBorderSize(10);
	//p18b->SetLogy();
	TH1F *mWbackFitb = new TH1F("mWbackFitb","mWbackFitb",200,0.0,3.0);
	EX->Draw("W>>+mWbackFitb",call_acept);//,chcer && ccointime1);
	mWbackFit->SetDirectory(0);
	MConfigHist(mWbackFitb,"EXP W background fit","W","counts",4,1);
	// Applying cut
	TH1F *mWbackFit2 = new TH1F("mWbackFit2b","mWbackFit2b",200,0.0,3.0);
	EX->Draw("W>>+mWbackFit2b",call1backFit,"same");
	mWbackFit2->SetDirectory(0);
	MConfigHist(mWbackFit2,"","","",3,1);
	f18->Draw("lsame");
	pad2png(c18,N+"_EXP_WbackFit.png");
	if(savePDF) c18->Print(N+".pdf");

//********************************
//********************************
//********************************
	// SIMC W
	TCanvas *c19 = new TCanvas("c19","c19",600,600);
	MConfigCanvas(c19,0,0);
	// SIMC
	TH1F *mSW = new TH1F("mSW","mSW",200,0.0,2.0);
	T->Draw("W>>+mSW",call_acept * Form("(Weight*%g)",normSIMC));
	mSW->SetDirectory(0);
	MConfigHist(mSW,titles+" - SIMC W","W","counts",4,0);
	//mSW->GetYaxis()->SetRangeUser(0.0,2000.0);
	// SIMC: Applying cut
	TH1F *mSW2 = new TH1F("mSW2","mSW2",200,0.0,2.0);
	T->Draw("W>>+mSW2",SIMCcall1 * Form("(Weight*%g)",normSIMC),"same");
	mSW2->SetDirectory(0);
	MConfigHist(mSW2,"","","",3,0);
	// EXP: Applying cut
	TH1F *mExpW2 = new TH1F("mExpW2","mExpW2",200,0.0,2.0);
	EX->Draw("W>>+mExpW2",call1 * Form("%g",ExpNorm),"same");
	mExpW2->SetDirectory(0);
	MConfigHist(mExpW2,"","","",2,0);
	// DUM: Applying cut
	TH1F *mDumW = new TH1F("mDumW2","mDumW2",200,0.0,2.0);
	DUM->Draw("W>>+mDumW2",call1 * Form("%g",ExpNormD),"same");
	mDumW2->SetDirectory(0);
	MConfigHist(mDumW2,"","","",5,0);
	c19->Update();
	Double_t c19xlim[2], c19ylim[3];
	c19xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c19xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c19ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c19ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l19;
	l19.SetTextSize(0.04);
	l19.DrawLatex(c19xlim[0]+0.1*(c19xlim[1]-c19xlim[0]),c19ylim[0]+0.8*(c19ylim[1]-c19ylim[0]),"#color[7]{SIMC (weighted)}");
	l19.DrawLatex(c19xlim[0]+0.1*(c19xlim[1]-c19xlim[0]),c19ylim[0]+0.75*(c19ylim[1]-c19ylim[0]),"#color[8]{SIMC (weighted + cut}");
	//l19.Draw();
	pad2png(c19,N+"_SIMC_W.png");
	if(savePDF) c19->Print(N+".pdf)");

//********************************
//********************************
//********************************
	// hsxptar and hsyptar
	TCanvas *c20 = new TCanvas("c20","c20",1200,600);
	//MConfigCanvas(c20,0,0);
	c20->Divide(2,1);
	TVirtualPad *p20 = c20->cd(1);
	//-> EXP.
	TH1F *phsxptar_all = new TH1F("phsxptar_all","Delta",100,-0.10,0.10);
	MConfigHist(phsxptar_all,titles + " - xptar distribution","xptar","Counts",41,0);
	EX->Draw("hsxptar>>+phsxptar_all",call_acept * Form("%g",ExpNorm));//,"","");
	//if(runN==47345) phsxptar_all->GetYaxis()->SetRangeUser(0.0,14000.0);
	phsxptar_all->SetDirectory(0);
	TH1F *phsxptar = new TH1F("phsxptar","DeltaCut",100,-0.10,0.10);
	MConfigHist(phsxptar,"","","",3,0);
	EX->Draw("hsxptar>>+phsxptar",call1 * Form("%g",ExpNorm),"goff");
	phsxptar->SetDirectory(0);
	//-> DUMMY
	//TH1F *phsxptarD_all = new TH1F("phsxptarD_all","DeltaD",100,-0.10,0.10);
	//MConfigHist(phsxptarD_all,"","","",1,0);
	//DUM->Draw("hsxptar>>+phsxptarD_all",Form("%g",ExpNormD));//,"","");
	//phsxptarD_all->SetDirectory(0);
	TH1F *phsxptarD = new TH1F("phsxptarD","DeltaCut",100,-0.10,0.10);
	MConfigHist(phsxptarD,"","","",8,0);
	//mhcer->SetFillStyle(3004);
	DUM->Draw("hsxptar>>+phsxptarD",call1 * Form("%g",ExpNormD),"same");
	phsxptarD->SetDirectory(0);
	//-> EXP - DUMMY
	TH1F *edhsxptar = phsxptar->Clone();
	edhsxptar->SetName("edhsxptar");
	edhsxptar->Add(phsxptarD,-1);
	MConfigHist(edhsxptar,"","","",1,0);
	edhsxptar->Draw("same");
	//-> SIMC
	TH1F *mhsxptar_all = new TH1F("mhsxptar_all","Delta",100,-0.10,0.10);
	MConfigHist(mhsxptar_all,titles+"","","",30,0);
	T->Draw("hsxptar>>+mhsxptar_all",call_acept * Form("Weight*%g",normSIMC),"same");//,"","");
	mhsxptar_all->SetDirectory(0);
	TH1F *mhsxptar = new TH1F("mhsxptar","DeltaCut",100,-0.10,0.10);
	MConfigHist(mhsxptar,"","","",2,0);
	//mhcer->SetFillStyle(3004);
	T->Draw("hsxptar>>+mhsxptar",SIMCcall1 * Form("(Weight*%g)",normSIMC),"same");
	mhsxptar->SetDirectory(0);
	///////////
	// YPTAR //
	///////////
	TVirtualPad *p20b = c20->cd(2);
	//-> EXP.
	TH1F *phsyptar_all = new TH1F("phsyptar_all","Delta",100,-0.05,0.05);
	MConfigHist(phsyptar_all,titles + " - yptar distribution","yptar","Counts",41,0);
	EX->Draw("hsyptar>>+phsyptar_all",call_acept * Form("%g",ExpNorm));//,"","");
	//if(runN==47345) phsyptar_all->GetYaxis()->SetRangeUser(0.0,14000.0);
	phsyptar_all->SetDirectory(0);
	TH1F *phsyptar = new TH1F("phsyptar","DeltaCut",100,-0.05,0.05);
	MConfigHist(phsyptar,"","","",3,0);
	EX->Draw("hsyptar>>+phsyptar",call1 * Form("%g",ExpNorm),"goff");
	phsyptar->SetDirectory(0);
	//-> DUMMY
	//TH1F *phsyptarD_all = new TH1F("phsyptarD_all","DeltaD",100,-0.05,0.05);
	//MConfigHist(phsyptarD_all,"","","",1,0);
	//DUM->Draw("hsyptar>>+phsyptarD_all",call_acept * Form("%g",ExpNormD));//,"","");
	//phsyptarD_all->SetDirectory(0);
	TH1F *phsyptarD = new TH1F("phsyptarD","DeltaCut",100,-0.05,0.05);
	MConfigHist(phsyptarD,"","","",8,0);
	//mhcer->SetFillStyle(3004);
	DUM->Draw("hsyptar>>+phsyptarD",call1 * Form("%g",ExpNormD),"same");
	phsyptarD->SetDirectory(0);
	//-> EXP - DUMMY
	TH1F *edhsyptar = phsyptar->Clone();
	edhsyptar->SetName("edhsyptar");
	edhsyptar->Add(phsyptarD,-1);
	MConfigHist(edhsyptar,"","","",1,0);
	edhsyptar->Draw("same");
	//-> SIMC
	TH1F *mhsyptar_all = new TH1F("mhsyptar_all","Delta",100,-0.05,0.05);
	MConfigHist(mhsyptar_all,titles+"","","",30,0);
	T->Draw("hsyptar>>+mhsyptar_all",call_acept * Form("Weight*%g",normSIMC),"same");//,"","");
	mhsyptar_all->SetDirectory(0);
	TH1F *mhsyptar = new TH1F("mhsyptar","DeltaCut",100,-0.05,0.05);
	MConfigHist(mhsyptar,"","","",2,0);
	//mhcer->SetFillStyle(3004);
	T->Draw("hsyptar>>+mhsyptar",SIMCcall1 * Form("(Weight*%g)",normSIMC),"same");
	mhsyptar->SetDirectory(0);
	//// Other things...
	Double_t c20xlim[2], c20ylim[2];
	c20->Update();
	c20xlim[0] = gPad->GetUxmin();//mhcer_all->GetXaxis()->GetXmin();
	c20xlim[1] = gPad->GetUxmax();//mhcer_all->GetXaxis()->GetXmax();
	c20ylim[0] = gPad->GetUymin(); //pow(10.0,gPad->GetUymin()); //pow(10,x) is the real value, but for latex comments we need the relative value
	c20ylim[1] = gPad->GetUymax(); //pow(10.0,gPad->GetUymax());
	TLatex l20;
	l20.SetTextSize(0.04);
	l20.DrawLatex(c20xlim[0]+0.05*(c20xlim[1]-c20xlim[0]),c20ylim[0]+0.95*(c20ylim[1]-c20ylim[0]),"#color[2]{SIMC (weighted + cut)}");
	l20.DrawLatex(c20xlim[0]+0.05*(c20xlim[1]-c20xlim[0]),c20ylim[0]+0.9*(c20ylim[1]-c20ylim[0]),"#color[1]{Experimental - Dummy (weighted + cut)}");
	l20.DrawLatex(c20xlim[0]+0.05*(c20xlim[1]-c20xlim[0]),c20ylim[0]+0.8*(c20ylim[1]-c20ylim[0]),"#color[30]{SIMC (weighted)}");
	l20.DrawLatex(c20xlim[0]+0.05*(c20xlim[1]-c20xlim[0]),c20ylim[0]+0.75*(c20ylim[1]-c20ylim[0]),"#color[41]{EXP (weighted)}");
	l20.DrawLatex(c20xlim[0]+0.05*(c20xlim[1]-c20xlim[0]),c20ylim[0]+0.7*(c20ylim[1]-c20ylim[0]),"#color[8]{DUMMY (weighted + cut)}");
	l20.Draw();
	cout << "EXP:  Selected/Total entries: " << EX->GetSelectedRows() << "/" << EX->GetEntries() << Form(" (%.1f\%)",100.0*EX->GetSelectedRows()/EX->GetEntries()) << endl
	     << "SIMC: Selected/Total entries:" << T->GetSelectedRows() << "/" <<  T->GetEntries() << Form(" (%.1f\%)",100.0*T->GetSelectedRows()/T->GetEntries()) << endl;
	// Printing
	pad2png(c20,N+"_hsxptar_hsyptar.png");
	if(savePDF) c20->Print(N+".pdf");

	logFile.close();

	if(saveRes) fPhi.close();
	return 0;
}

int main(int argc, char *argv[]) {
        TString inData = argv[1];

        return MlookSIMC();//UncertaintyPion(inData);
}

//********************************
//********************************
//********************************
	// missmass vs cointime 2D-histogram with nboxes boxes centered at xc
int MDraw2DMM(TString name, TTree *EX, TCut mcut, Int_t nboxes, Double_t *xc, Double_t *xw, Double_t *yc, Double_t *yw, Double_t *range, TH2F *a, Int_t nbin, Double_t maxZ) {
	TCanvas *c = new TCanvas(name,name,1200,600);
	MConfigCanvas(c,0,0);
	a = new TH2F("a","a",nbin,range[0],range[2],nbin,range[1],range[3]);
	gStyle->SetNumberContours(250);
	if(maxZ>0.0) a->GetZaxis()->SetRangeUser(0.0, maxZ);
	gStyle->SetOptStat(0);
	EX->Draw("missmass:cointime>>a",mcut,"COLZ");
	a->GetYaxis()->SetTitle("Missing mass / GeV");
	a->GetXaxis()->SetTitle("Coincidence time / ns");
	a->SetTitle("");
	c->Update();
	// Drawing boxes
	TLatex l;
	l.SetTextSize(0.03);
	TBox b1[nboxes];
	for(Int_t ii=0; ii<nboxes; ii++) {
		// Drawing box
		b1[ii].SetLineColor(1);
		b1[ii].SetLineWidth(2);
		b1[ii].SetFillStyle(0);
		b1[ii].DrawBox(xc[ii]-xw[ii],yc[ii]-yw[ii],xc[ii]+xw[ii],yc[ii]+yw[ii]);
		// Drawing numbers
		if(nboxes>1) l.DrawLatex(TMath::Max(gPad->GetUxmin(),xc[ii]-xw[ii])+0.1*(TMath::Min(gPad->GetUxmax(),xc[ii]+xw[ii])-TMath::Max(gPad->GetUxmin(),xc[ii]-xw[ii])),TMath::Max(gPad->GetUymin(),yc[ii]-yw[ii])+0.65*(TMath::Min(gPad->GetUymax(),yc[ii]+yw[ii])-TMath::Max(gPad->GetUymin(),yc[ii]-yw[ii])),Form("%d",ii));
	}
	//pad2png(c,name+".png");
	//c4c->Print(name+".eps");
	a->SetDirectory(0);
	return 0;
}

// To read Scaler
Double_t ReadScaler(TString filenameInit, Int_t runN) {
	ifstream f_read;
	f_read.open("./files/Scalars/"+filenameInit+"_"+Form("%d",runN)+".dat");
	Double_t val;
	f_read >> val;
	f_read.close();
	return val;
}
