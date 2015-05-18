#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>

// old
/*
Double_t fminim(Double_t *x, Double_t *par, Bool_t evF=1);
double minim(const double *x);
//double minim(const double *x, const double *par=0);
Double_t ModelSig(Double_t thetacm, Double_t phicm, Double_t t, Double_t q2_gev, Double_t s_gev, Double_t w2_gev, Double_t eps, Double_t mtar_gev, Double_t *par, Bool_t fit=1, Bool_t print_separated=0);

Double_t minEvts;
*/

// start here
Double_t ModelSigUsedSIMC(Double_t thetacm, Double_t phicm, Double_t t, Double_t q2_gev, Double_t s_gev, Double_t w2_gev, Double_t eps, Double_t mtar_gev);
Double_t FitModel(Double_t thetacm, Double_t phicm, Double_t t, Double_t q2_gev, Double_t s_gev, Double_t w2_gev, Double_t eps, Double_t mtar_gev, Double_t *par);
Double_t funcFitModel(Double_t *x, Double_t *par);
Double_t *MDot (Int_t n, Bool_t *a_ok, Double_t *a_error);

Int_t npts = 0;

const int nvariables = 7; // number of inputs in the cross section, i.e., q2, w2, eps, ...
const int nfitpars = 9; // number of parameters to be fitted by the function

Bool_t *OK;
Int_t *runN, *folderN;
Double_t *centerPhi, *meanTheta, *meanT, *meanQ2, *meanS, *meanW2, *meanEPS, *mtar, *NevtsEXP, *YieldEXP, *ratioEXPSIMC, *EXP, *errorEXP;

//______________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
   Int_t nevts = npts;
//calculate chisquare
   Double_t chisq = 0;
   Double_t delta;
   Double_t val;
   for (Int_t i=0;i<nevts; i++) {
//	cout << "PHI: " << centerPhi[i] << endl;
     if(OK[i]) {
	     val = FitModel(meanTheta[i], centerPhi[i], meanT[i], meanQ2[i], meanS[i], meanW2[i], meanEPS[i], mtar[i], par);
	     delta  = (EXP[i]-val)/errorEXP[i];
	     //cout << "i=" << i << "\tEXP = " << EXP[i] << "\tval = " << val << endl;
	     chisq += delta*delta;
     }
   }
   cout << "chisq = " << chisq << endl;
   f = chisq;
}
 
int Mfit(TString infolder="Rosenbluth2015_006", Int_t minNumEvt=50) {

	//////////////////////
	/// UPDATE HERE!!! ///
	//////////////////////
	Double_t FPI2_lowQ2_tcuts[]  = {0.180,0.270,0.400};
	Double_t FPI2_highQ2_tcuts[] = {0.330,0.415,0.530};
	Int_t nphi = 10; // VERY important!
	//////////////////////
	//    until here   ///
	//////////////////////

	// Creating phi vector
	Double_t *phi = new Double_t[nphi]; // center value of phi
	for(Int_t i=0; i<nphi; i++) phi[i] = (360.0/nphi/2.0)+(i*360.0/nphi);

	// Getting all data into arrays
	OK           = new Bool_t[1000];
	runN         = new Int_t[1000];
	folderN      = new Int_t[1000];
	centerPhi    = new Float_t[1000];
	meanTheta    = new Float_t[1000];
	meanT        = new Float_t[1000];
	meanQ2       = new Float_t[1000];
	meanS        = new Float_t[1000];
	meanW2       = new Float_t[1000];
	meanEPS      = new Float_t[1000];
	mtar         = new Float_t[1000];
	NevtsEXP     = new Float_t[1000];
	YieldEXP     = new Float_t[1000];
	ratioEXPSIMC = new Float_t[1000];
	EXP          = new Float_t[1000];
	errorEXP     = new Float_t[1000];
	TTree *T = new TTree("T","T");

	// Getting all csv files from infolder
	TSystemDirectory dir(infolder,infolder);
	TList *files = dir.GetListOfFiles();
	if(files==0) {
		cout << "Could not find folder: " << infolder << endl;
		return 1;
	}
	TSystemFile *file, *file1, *file2;
	TString fname, fname1, fname2;
	TIter next(files);
	Int_t nfiles=0;
	Int_t runNread[1000];
	Int_t countFolder=0;
	Double_t **meanAll = new Double_t*[1000];
	TString *Folder = new TString[1000];
	TString tmpStr;
	Int_t fdigit, ldigit; // first and last numeric digit in fname
	// Getting total number of runs, and run numbers
	while ((file=(TSystemFile*)next())) {
		fname = file->GetName();
		if(file->IsDirectory() && !fname.BeginsWith(".")) { // opening higheps and loweps folders
			TSystemDirectory dir1(infolder+"/"+fname,infolder+"/"+fname);
			TList *files1 = dir1.GetListOfFiles();
			TIter next1(files1);
			while ((file1=(TSystemFile*)next1())) {
				fname1 = file1->GetName();
				if(file1->IsDirectory() && fname1.BeginsWith("t_")) { // opening t folders
					TSystemDirectory dir2(infolder+"/"+fname+"/"+fname1,infolder+"/"+fname+"/"+fname1);
					TList *files2 = dir2.GetListOfFiles();
					TIter next2(files2);
					TTree *ntu = new TTree("ntu","ntu");//,"runN:centerPhi:meanTheta:meanT:meanQ2:meanS:meanW2:meanEPS:mtar:NevtsEXP:YieldEXP:YieldSIMC:ratioEXPSIMC");
					while ((file2=(TSystemFile*)next2())) {
						fname2 = file2->GetName();
						if (!file2->IsDirectory() && fname2.EndsWith("scaleOutput.csv")) {
							cout << fname2 << endl;
							// Extract run number
							for(fdigit=0; fdigit<fname2.Length() && (fname2(fdigit)<'0' || fname2(fdigit)>'9'); fdigit++);
							for(ldigit=fdigit; ldigit<fname2.Length() && (fname2(ldigit)>='0' && fname2(ldigit)<='9'); ldigit++);
							tmpStr = fname2(fdigit,ldigit-fdigit);
							runNread[nfiles] = tmpStr.Atoi();
							// Reading file
							T->ReadFile(infolder+"/"+fname+"/"+fname1+"/"+fname2,"runN/I:centerPhi/F:meanTheta:meanT:meanQ2:meanS:meanW2:meanEPS:mtar:NevtsEXP:YieldEXP:YieldSIMC:ratioEXPSIMC",',');
							ntu->ReadFile(infolder+"/"+fname+"/"+fname1+"/"+fname2,"runN/I:centerPhi/F:meanTheta:meanT:meanQ2:meanS:meanW2:meanEPS:mtar:NevtsEXP:YieldEXP:YieldSIMC:ratioEXPSIMC",','); // for mean value for this kinematic point (temporary)
							cout << runNread[nfiles] << endl;
							nfiles++;
						}
						else if (!file2->IsDirectory() && fname2.EndsWith(".csv") && fname2.BeginsWith("summary_")) {
							// read summary here, if needed
						}
					}
					// Getting mean values
					Folder[countFolder] = fname+" - "+fname1;
					cout << "Done reading: " << Folder[countFolder] << endl;
					for(Int_t ii=(countFolder*nphi); ii<(countFolder*nphi+nphi); ii++) {
						ntu->Draw("NevtsEXP:meanTheta:meanT:meanQ2",Form("NevtsEXP>=%d && abs(centerPhi-%f)<0.1",minNumEvt,phi[ii%nphi]),"goff");
//						ntu->Scan("NevtsEXP:meanTheta:meanT:meanQ2",Form("NevtsEXP>=%d && abs(centerPhi-%f)<0.1",minNumEvt,phi[ii%nphi]));
						folderN[ii] = countFolder;
						centerPhi[ii] = phi[ii%nphi]; // phicm
						if(ntu->GetSelectedRows()>0) {
							meanTheta[ii] = TMath::Mean(ntu->GetSelectedRows(),ntu->GetV2(),ntu->GetV1()); // thetacm
							meanT[ii] = TMath::Mean(ntu->GetSelectedRows(),ntu->GetV3(),ntu->GetV1()); // t
							meanQ2[ii] = TMath::Mean(ntu->GetSelectedRows(),ntu->GetV4(),ntu->GetV1()); // Q2
							ntu->Draw("NevtsEXP:meanS:meanW2:meanEPS",Form("NevtsEXP>=%d && abs(centerPhi-%f)<0.1",minNumEvt,phi[ii%nphi]),"goff");
//							ntu->Scan("NevtsEXP:meanS:meanW2:meanEPS",Form("NevtsEXP>=%d && abs(centerPhi-%f)<0.1",minNumEvt,phi[ii%nphi]));
							meanS[ii] = TMath::Mean(ntu->GetSelectedRows(),ntu->GetV2(),ntu->GetV1()); // S
							meanW2[ii] = TMath::Mean(ntu->GetSelectedRows(),ntu->GetV3(),ntu->GetV1()); // W2
							meanEPS[ii] = TMath::Mean(ntu->GetSelectedRows(),ntu->GetV4(),ntu->GetV1()); // eps
							ntu->Draw("NevtsEXP:mtar:YieldEXP:ratioEXPSIMC",Form("NevtsEXP>=%d && abs(centerPhi-%f)<0.1",minNumEvt,phi[ii%nphi]),"goff");
							mtar[ii] = TMath::Mean(ntu->GetSelectedRows(),ntu->GetV2(),ntu->GetV1()); // mtar
							YieldEXP[ii] = TMath::Mean(ntu->GetSelectedRows(),ntu->GetV3(),ntu->GetV1()); // YieldEXP
							ratioEXPSIMC[ii] = TMath::Mean(ntu->GetSelectedRows(),ntu->GetV4(),ntu->GetV1()); // ratioEXPSIMC
							NevtsEXP[ii] = ntu->GetSelectedRows()*TMath::Mean(ntu->GetSelectedRows(),ntu->GetV1()); // sum of NevtsEXP
							EXP[ii] = ratioEXPSIMC[ii] * ModelSigUsedSIMC(meanTheta[ii], centerPhi[ii], meanT[ii], meanQ2[ii], meanS[ii], meanW2[ii], meanEPS[ii], mtar[ii]);
							errorEXP[ii] = EXP[ii] * (1.0/sqrt(NevtsEXP[ii])); // should also include std of data used in the mean
							OK[ii]=1;
						}
						else {
							meanTheta[ii] = 0;
							meanT[ii] = 0;
							meanQ2[ii] = 0;
							meanS[ii] = 0;
							meanW2[ii] = 0;
							meanEPS[ii] = 0;
							mtar[ii] = 0;
							YieldEXP[ii] = 0;
							ratioEXPSIMC[ii] = 0;
							NevtsEXP[ii] = 0;
							EXP[ii] = -1;
							errorEXP[ii] = 0;
							OK[ii]=0;
						}
						// Printing to check
//						cout << "-------" << endl << "  -> phi = " << centerPhi[ii] << endl << "  -> folderN = " << folderN[ii] << endl << "  -> meanTheta = " << meanTheta[countFolder*nphi+ii] << endl << "  -> meanT = " << meanT[ii] << endl << "  -> meanQ2 = " << meanQ2[ii] << endl << "  -> meanS = " << meanS[ii] << endl << "  -> meanW2 = " << meanQ2[ii] << endl << "  -> meanEPS = " << meanEPS[ii] << endl << "  -> EXP = " << EXP[ii] << endl << "  -> errorEXP = " << errorEXP[ii] << endl << "  -> OK = " << OK[ii] << endl;;
					}
					countFolder++;
					delete ntu;
				}
			}
		}
	}
	npts = countFolder*nphi;

	/////////////////
	////// FIT //////
	/////////////////

	// https://root.cern.ch/drupal/content/how-fit-histograms-or-data-points

	// Documentation:
	// http://cern-tex.web.cern.ch/cern-tex/minuit/node18.html

	// initialize TMinuit with a maximum of nfitpars params
	TMinuit *gMinuit = new TMinuit(nfitpars);  
	gMinuit->SetFCN(fcn);
 
	Double_t arglist[10];
	Int_t ierflg = 0;
 
	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
 
	// Set starting values and step sizes for parameters
	Double_t vstart[nfitpars] = {1.2,0.0001,0.00001,1,1,1,1,1,1};
	//Double_t vstart[nfitpars] = {1.2,0.0,0.0,0.5,2.1,0.5,2,2,2};
	Double_t step[nfitpars] = {0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001,0.00001};
	gMinuit->mnparm(0, "par0", vstart[0], step[0], 0,5,ierflg);
	gMinuit->mnparm(1, "par1", vstart[1], step[1], 0,20,ierflg);
	gMinuit->mnparm(2, "par2", vstart[2], step[2], 0,20,ierflg);
	gMinuit->mnparm(3, "par3", vstart[3], step[3], 0,10,ierflg);
	gMinuit->mnparm(4, "par4", vstart[4], step[4], 0,10,ierflg);
	gMinuit->mnparm(5, "par5", vstart[5], step[5], 0,10,ierflg);
	gMinuit->mnparm(6, "par6", vstart[6], step[6], 0,10,ierflg);
	gMinuit->mnparm(7, "par7", vstart[7], step[7], 0,10,ierflg);
	gMinuit->mnparm(8, "par8", vstart[8], step[8], 0,10,ierflg);
 
	// Now ready for minimization step
	arglist[0] = 500000;
	arglist[1] = 1;

	// Only sigA first
	gMinuit->FixParameter(1);
	gMinuit->FixParameter(2);
	gMinuit->FixParameter(5);
	gMinuit->FixParameter(6);
	gMinuit->FixParameter(7);
	gMinuit->FixParameter(8);
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	gMinuit->Release(1);
	gMinuit->Release(2);
	gMinuit->Release(5);
	gMinuit->Release(6);
	gMinuit->Release(7);
	gMinuit->Release(8);

	// Fix sigA, and fit the others
	gMinuit->FixParameter(0);
	gMinuit->FixParameter(3);
	gMinuit->FixParameter(4);
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	gMinuit->Release(0);
	gMinuit->Release(3);
	gMinuit->Release(4);
	
	// Fit all together now
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	arglist[0]=0;
	arglist[1]=51;
	arglist[2]=0;
	arglist[3]=1;
	//gMinuit->mnexcm("SCA", arglist, 4, ierflg);
	arglist[0] = 500000;

	// Evaluate function
	Double_t grad[nfitpars];
	Double_t fval;
	Double_t parx[] = {2.74616e-01,0,0,0,0,0,0,0,0};
	gMinuit->Eval(1, grad, fval, parx, 2);
	cout << "fval = " << fval << endl;

//	for(Int_t i=3; i<=8; i++) gMinuit->FixParameter(i);
//	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
//	for(Int_t i=3; i<=8; i++) gMinuit->Release(i);
//	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
 
	// Print results
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	gMinuit->mnprin(3,amin);

	// Get parameters
	Double_t *fittedVal = new Double_t[nfitpars];
	Double_t *fittedValErr = new Double_t[nfitpars];
	for(Int_t i=0; i<nfitpars; i++)	gMinuit->GetParameter(i,fittedVal[i],fittedValErr[i]);
//	cout << "P3 = " << fittedVal[2] << endl;
	
	// DONE WITH FIT
	cout << "END OF FIT" << endl;
	for(Int_t i=0; i<36; i++) cout << "*";
	cout << endl << endl;




//	// Adding branch with original cross section, used in simc before this fit
//	Float_t val;
//	TBranch *sigModel = T->Branch("sigModel",&val,"sigModel/F");
//	for(Int_t i=0; i<T->GetEntries(); i++) {
//		T->GetEntry(i);
//		//cout << "ok? ";
//		val = ModelSigUsedSIMC(a_thetacm, a_phicm, a_t, a_q2, a_s, a_w2, a_eps, a_mtar);
//		sigModel->Fill();
////		cout << val << endl;
//	}

	//T->Scan("runN:centerPhi:NevtsEXP:sigModel");

	// Preparing plots
/*
	///////////
	// Plot 0 - individual runs
	///////////
	TGraphErrors **gr0 = new TGraphErrors[nfiles];
	TCanvas *c0 = new TCanvas("c0","c0",970,600);
	c0->SetTicks(1,1);
	TLegend *leg0 = new TLegend(0.7,0.15,0.88,0.43);
	for(i=0; i<nfiles; i++) {
		T->Draw("centerPhi:(ratioEXPSIMC*sigModel):(ratioEXPSIMC*sigModel/TMath::Sqrt(NevtsEXP))",Form("runN==%d && NevtsEXP>=%d",runN[i],minNumEvt),"goff");
		gr0[i] = new TGraphErrors(T->GetSelectedRows(),T->GetV1(),T->GetV2(),0,T->GetV3());
		gr0[i]->SetMarkerStyle(20);
		gr0[i]->SetMarkerColor(i+1);
		if(i==0) {
			gr0[i]->GetXaxis()->SetTitle("Phi (deg)");
			gr0[i]->GetXaxis()->SetLimits(0.0,360.0);
			gr0[i]->GetYaxis()->SetTitle("#sigma");
			T->Draw("ratioEXPSIMC*sigModel","","goff");
			gr0[i]->GetYaxis()->SetRangeUser(0.0,1.1*TMath::MaxElement(T->GetSelectedRows(),T->GetV1()));
			gr0[i]->SetTitle("");
			gr0[i]->Draw("ap");
		}
		else gr0[i]->Draw("psame");
		leg0->AddEntry(gr0[i],Form("run %d",runN[i]),"p");
	}
	leg0->Draw();
	c0->Print(infolder+"/runs.png");
*/
	///////////
	// Plot 1 - all together
	///////////
	Int_t nfit = countFolder;

	// FIRST FIT - ALL SEPARATED
	TCanvas **c1 = new TCanvas*[nfit];
	TGraphErrors **gr1 = new TGraphErrors*[nfit];
	TF1 **f1 = new TF1*[nfit];
	Double_t **parf1 = new Double_t*[nfit];
	TString title;
	for(Int_t i=0; i<nfit; i++) {
		cout << "-> Plotting " << Folder[i] << endl;
		title = Folder[i];
		title.ReplaceAll(' ',"");
		cout << "TITLE: " << title << endl;
		c1[i] = new TCanvas("c1","c1",970,600);
		gr1[i] = new TGraphErrors(nphi,phi,&EXP[i*nphi],0,&errorEXP[i*nphi]);
		gr1[i]->SetMarkerStyle(7);
		gr1[i]->SetMarkerSize(12);
		gr1[i]->SetMarkerColor(1);
		gr1[i]->SetTitle(title+";Phi (deg);Cross section (nb/GeV^{2})");
		gr1[i]->GetXaxis()->SetLimits(0.0,360.0);
		gr1[i]->GetYaxis()->SetRangeUser(0.0,150);//1.1*TMath::MaxElement(nphi,sigAll));
		gr1[i]->Draw("ap");
		f1[i] = new TF1("f1",&funcFitModel,0.0,360.0,nvariables+nfitpars);
		parf1[i] = new Double_t[nvariables+nfitpars];
		parf1[i][0] = TMath::Mean(nphi,&meanTheta[i*nphi],MDot(nphi,&OK[i*nphi],&NevtsEXP[i*nphi])); // meanTheta
		parf1[i][1] = TMath::Mean(nphi,&meanT[i*nphi],MDot(nphi,&OK[i*nphi],&NevtsEXP[i*nphi])); // t
		parf1[i][2] = TMath::Mean(nphi,&meanQ2[i*nphi],MDot(nphi,&OK[i*nphi],&NevtsEXP[i*nphi])); // q2
		parf1[i][3] = TMath::Mean(nphi,&meanS[i*nphi],MDot(nphi,&OK[i*nphi],&NevtsEXP[i*nphi])); // s
		parf1[i][4] = TMath::Mean(nphi,&meanW2[i*nphi],MDot(nphi,&OK[i*nphi],&NevtsEXP[i*nphi])); // w2
		parf1[i][5] = TMath::Mean(nphi,&meanEPS[i*nphi],MDot(nphi,&OK[i*nphi],&NevtsEXP[i*nphi])); // eps
		parf1[i][6] = TMath::Mean(nphi,&mtar[i*nphi],MDot(nphi,&OK[i*nphi],&NevtsEXP[i*nphi])); // mtar
		for(Int_t ii=nvariables; ii<(nvariables+nfitpars); ii++) {
			parf1[i][ii] = fittedVal[ii-nvariables];
			cout << parf1[i][ii] << "//";
		}
		cout << endl;
		f1[i]->SetParameters(parf1[i]);
		cout << "PARS: ";
		for(Int_t ii=0; ii<(nvariables+nfitpars); ii++) cout << parf1[i][ii] << ", ";
		cout << endl;
		//Double_t ModelSig(Double_t thetacm, Double_t phicm, Double_t t, Double_t q2_gev, Double_t s_gev, Double_t w2_gev, Double_t eps, Double_t mtar_gev, Double_t *par, Bool_t fit) {
		f1[i]->SetLineColor(2);
		//gr1->Fit("f1","RM+");
		f1[i]->Draw("lsame");
		c1[i]->Print(infolder+"/"+title+".png");
		// Printing separated cross sections
		//ModelSig(parf1[0], 0, parf1[1], parf1[2], parf1[3], parf1[4], parf1[5], parf1[6], &parf1[nfitpars], 0, 1);

		// Printing
		//cout << "sigA  = " << f1->GetParameter(0)*(2.0*TMath::Pi()*1.0e6) << endl;
		//cout << "sigLT = " << f1->GetParameter(1)/(TMath::Sqrt(2.0*eps*(1.0+eps))/(2.0*TMath::Pi()*1.0e6)) << endl;
		//cout << "sigTT = " << f1->GetParameter(2)/(eps/(2.0*TMath::Pi()*1.0e6)) << endl;
	}

	// SECOND FIT - ALL IN SAME CANVAS
	TCanvas *c2 = new TCanvas("c2","c2",3000,1000);
	c2->Divide(4,2,0,0);
	TGraphErrors **gr2 = new TGraphErrors*[nfit];
	TF1 **f2 = new TF1*[nfit];
	Double_t **parf2 = new Double_t*[nfit];
	TString title;
	TVirtualPad *p2;// = new TVirtualPad[nfit];
	for(Int_t i=0; i<nfit; i++) {
		p2 = c2->cd(i+1);
		p2->SetTicks(1,1);
		title = Folder[i];
		title.ReplaceAll(' ',"");
		gr2[i] = new TGraphErrors(nphi,phi,&EXP[i*nphi],0,&errorEXP[i*nphi]);
		gr2[i]->SetMarkerStyle(21);
		gr2[i]->SetMarkerColor(1);
		gr2[i]->SetTitle(title+";Phi (deg);Cross section (nb/GeV^{2})");
		gr2[i]->GetXaxis()->SetLimits(0.0,360.0);
		gr2[i]->GetYaxis()->SetRangeUser(0.0,150);//1.1*TMath::MaxElement(nphi,sigAll));
		gr2[i]->Draw("ap");
		f2[i] = new TF1("f1",&funcFitModel,0.0,360.0,nvariables+nfitpars);
		parf2[i] = new Double_t[nvariables+nfitpars];
		for(Int_t ii=0; ii<(nvariables+nfitpars); ii++) parf2[i][ii] = parf1[i][ii];
		f2[i]->SetParameters(parf2[i]);
		f2[i]->SetLineColor(2);
		f2[i]->Draw("lsame");
	}
	c2->Print(infolder+"/all.png");

	return 0;
}

Double_t ModelSigUsedSIMC(Double_t thetacm, Double_t phicm, Double_t t, Double_t q2_gev, Double_t s_gev, Double_t w2_gev, Double_t eps, Double_t mtar_gev) {
	// Rosenbluth2015_000
        //Float_t sigA = 1200.0/(1.0+0.53*q2_gev);
        //Float_t sigLT = 1.0*exp(-2.1*(t-2.0))*(1.0/(q2_gev+0.785)**2)*(1.0/(sqrt(w2_gev)*(w2_gev-0.93827**2)))*(2.0+eps)/3.0;
        //Float_t sigTT = 1.0*exp(-2.1*(t-2.0))*(1.0/(q2_gev+0.785)**2)*(1.0/(sqrt(w2_gev)*(w2_gev-0.93827**2)))*(2.0+eps)/3.0;

	// Rosenbluth2015_004
        //Float_t sigA = 1.495660*1000.0*1200.0/(1.0+0.53*q2_gev);
        //Float_t sigLT = 3.927480*100000.0*1.0*exp(-2.1*(t-2.0))*(1.0/(q2_gev+0.785)**2)*(1.0/(sqrt(w2_gev)*(w2_gev-0.93827**2)))*(2.0+eps)/3.0;
        //Float_t sigTT = 4.5764808*100000.0*1.0*exp(-2.1*(t-2.0))*(1.0/(q2_gev+0.785)**2)*(1.0/(sqrt(w2_gev)*(w2_gev-0.93827**2)))*(2.0+eps)/3.0;

	// Rosenbluth2015_005
//	if(!fit && !print_separated) {
//		cout << "YES!" << endl;
//	        Double_t sigA = 1.495660*1200.0/(1.0+0.53*q2_gev);
//        	Double_t sigLT = 0.3927480*exp(-2.1*(t-2.0))*(1.0/(q2_gev+0.785)**2)*(1.0/(sqrt(w2_gev)*(w2_gev-0.93827**2)))*(2.0+eps)/3.0;
//	        Double_t sigTT = 0.45764808*1.0*exp(-2.1*(t-2.0))*(1.0/(q2_gev+0.785)**2)*(1.0/(sqrt(w2_gev)*(w2_gev-0.93827**2)))*(2.0+eps)/3.0;
//	}
//	else { // FIT to Rosenbluth2015_005
//	        Double_t sigA = par[0]*1200.0/(1.0+0.53*q2_gev);
//	        Double_t sigLT = par[1]*exp(-2.1*(t-2.0))*(1.0/(q2_gev+0.785)**2)*(1.0/(sqrt(w2_gev)*(w2_gev-0.93827**2)))*(2.0+eps)/3.0;
//	        Double_t sigTT = par[2]*exp(-2.1*(t-2.0))*(1.0/(q2_gev+0.785)**2)*(1.0/(sqrt(w2_gev)*(w2_gev-0.93827**2)))*(2.0+eps)/3.0;
//		cout << par[0] << " -- " << par[1] << " -- " << par[2] << endl;
//	}

	// Rosenbluth2015_006
//	cout << "YES!" << endl;
        Double_t sigA = 1.2/(1.0+0.53*q2_gev);
       	Double_t sigLT = 0.3927480*exp(-2.1*(t-2.0))*(1.0/(q2_gev+0.785)**2)*(1.0/(sqrt(w2_gev)*(w2_gev-0.93827**2)))*(2.0+eps)/3.0;
        Double_t sigTT = 0.45764808*1.0*exp(-2.1*(t-2.0))*(1.0/(q2_gev+0.785)**2)*(1.0/(sqrt(w2_gev)*(w2_gev-0.93827**2)))*(2.0+eps)/3.0;

//	if(print_separated) {
//		cout << "sigA  = " << sigA  << endl;
//		cout << "sigLT = " << sigLT << endl;
//		cout << "sigTT = " << sigTT << endl;
//	}

	// Deg to Rad
	phicm = phicm*TMath::DegToRad();

	//Together
	Double_t sig219=(sigA+eps*cos(2.*phicm)*sigTT+sqrt(2.0*eps*(1.+eps))*cos(phicm)*sigLT);
	Double_t sig=sig219/2./TMath::Pi()/1.e06; 

	sig *= 1.e9; // ub/MeV^2 --> nb/GeV^2
	
	return sig;
}

Double_t FitModel(Double_t thetacm, Double_t phicm, Double_t t, Double_t q2_gev, Double_t s_gev, Double_t w2_gev, Double_t eps, Double_t mtar_gev, Double_t *par) {

	// Rosenbluth2015_006
        Double_t sigA = par[0]*exp(-par[4]*t)/(1.0+par[3]*q2_gev);
       	Double_t sigLT = par[1]*exp(-par[5]*t)/(1.0+q2_gev*par[6]);
        Double_t sigTT = par[2]*exp(-par[7]*t)/(1.0+q2_gev*par[8]);

//	if(print_separated) {
//		cout << "sigA  = " << sigA  << endl;
//		cout << "sigLT = " << sigLT << endl;
//		cout << "sigTT = " << sigTT << endl;
//	}

	// Deg to Rad
	phicm = phicm*TMath::DegToRad();

	//Together
	Double_t sig219=(sigA+eps*cos(2.*phicm)*sigTT+sqrt(2.0*eps*(1.+eps))*cos(phicm)*sigLT);
	Double_t sig=sig219/2./TMath::Pi()/1.e06; 
	
	sig *= 1.e9; // ub/MeV^2 --> nb/GeV^2

	return sig;
}

Double_t funcFitModel(Double_t *x, Double_t *par) {
//Double_t FitModel(Double_t thetacm, Double_t phicm, Double_t t, Double_t q2_gev, Double_t s_gev, Double_t w2_gev, Double_t eps, Double_t mtar_gev, Double_t *par) {
	return FitModel(par[0], x[0], par[1], par[2], par[3], par[4], par[5], par[6], &par[7]);
}

Double_t *MDot (Int_t n, Bool_t *a_ok, Double_t *a_error) {
	Double_t *res = new Double_t[n];
	for(Int_t i=0; i<n; i++) {
		if(a_ok[i]) res[i] = a_error[i];
		else res[i] = 0.0;
	}
	return res;
}
