int Q2eps(Int_t run1a=47056, Int_t run1b=47154, Int_t run2a=47358, Int_t run2b=47428, Int_t run3a=47587, Int_t run3b=47613, Int_t run4a=47687, Int_t run4b=47736) {

	// Important defininitions
	TString EXinFolder = "../files/";
	TString EXtype = "coin";
	TString treename = "h9500";
	Int_t NotAnalyze[] = {47368, 47372, 47373, 47374, 47375, 47376, 47377, 47378, 47379, 47380, 47381, 47382, 47383, 47384, 47385, 47386, 47387, 47388, 47389, 47390, 47391, 47392, 47397, 47403, 47406, 47407, 47408, 47409, 47410, 47411, 47412, 47413, 47414, 47415, 47416, 47417, 47418, 47420, 47425, 47088, 47089, 47090, 47091, 47092, 47093, 47094, 47095, 47096, 47097, 47098, 47099, 47100, 47101, 47102, 47103, 47104, 47105, 47106, 47107, 47108, 47109, 47110, 47111, 47112, 47113, 47114, 47115, 47116, 47117, 47118, 47368, 47372, 47373, 47374, 47375, 47376, 47377, 47378, 47379, 47380, 47381, 47382, 47383, 47384, 47385, 47386, 47387, 47388, 47389, 47390, 47391, 47392, 47397, 47403, 47406, 47407, 47408, 47409, 47410, 47411, 47412, 47413, 47414, 47415, 47416, 47417, 47418, 47420, 47425, 47592, 47595, 47604, 47694, 47697, 47729, 47736, 47598, 47599, 47614, 47615, 47616, 47617, 47618, 47619, 47620, 47621, 47622, 47623, 47624, 47625, 47626, 47627, 47628, 47629, 47630, 47631, 47632, 47633, 47634, 47635, 47636, 47637, 47638, 47639, 47640, 47641, 47642, 47643, 47644, 47645, 47646, 47647, 47648, 47649, 47650, 47651, 47652, 47653, 47654, 47655, 47656, 47657, 47658, 47659, 47660, 47661, 47662, 47663, 47664, 47665, 47666, 47667, 47668, 47669, 47670, 47671, 47672, 47673, 47674, 47675, 47676, 47677, 47678, 47679, 47680, 47681, 47682, 47683, 47684, 47685, 47693, 47694, 47702, 47703, 47704, 47705, 47706, 47707, 47708, 47709, 47710, 47711, 47712, 47713, 47714, 47715, 47716, 47717, 47718, 47719, 47720, 47721, 47722, 47723};

	// Variables
	Bool_t found;
	TString tempname;

	// Loading run1
	cout << "**************" << endl << "*** RUN 1 ***" << endl << "*************" << endl;
	TChain *t1 = new TChain(treename);
	for(Int_t i=run1a; i<=run1b; i++) {
		found=0;
		for(Int_t j=0;j<(sizeof(NotAnalyze)/sizeof(NotAnalyze[0])) && !found;j++)
			if(i==NotAnalyze[j]) found = 1;
		if(!found) {
			tempname = EXinFolder + EXtype + Form("%d",i) + ".root";
			if(!gSystem->AccessPathName(tempname)) {
				// to count raw number of events of individual files
				f = TFile::Open(tempname);
				ttemp = (TTree*)f->Get(treename);
				// Chaining file
				cout << "\t- File: " << tempname << " - " << t1->Add(tempname) << " - " << ttemp->GetEntries() << " entries" << endl;
			}
			Int_t split=1;
			tempname = EXinFolder + EXtype + Form("%d.%d",i,split) + ".root";
			while(!gSystem->AccessPathName(tempname)) {
				// to count raw number of events of individual files
				f = TFile::Open(tempname);
				ttemp = (TTree*)f->Get(treename);
				// Chaining file
				cout << "\t- File: " << tempname << " - " << t1->Add(tempname) << " - " << ttemp->GetEntries() << " entries" << endl;
				// Preparing next split
				tempname = EXinFolder + EXtype + Form("%d.%d",i,++split) + ".root";
			}
		}
	}

	// Loading run2
	cout << "**************" << endl << "*** RUN 2 ***" << endl << "*************" << endl;
	TChain *t2 = new TChain(treename);
	for(Int_t i=run2a; i<=run2b; i++) {
		found=0;
		for(Int_t j=0;j<(sizeof(NotAnalyze)/sizeof(NotAnalyze[0])) && !found;j++)
			if(i==NotAnalyze[j]) found = 1;
		if(!found) {
			tempname = EXinFolder + EXtype + Form("%d",i) + ".root";
			if(!gSystem->AccessPathName(tempname)) {
				// to count raw number of events of individual files
				f = TFile::Open(tempname);
				ttemp = (TTree*)f->Get(treename);
				// Chaining file
				cout << "\t- File: " << tempname << " - " << t2->Add(tempname) << " - " << ttemp->GetEntries() << " entries" << endl;
			}
			Int_t split=1;
			tempname = EXinFolder + EXtype + Form("%d.%d",i,split) + ".root";
			while(!gSystem->AccessPathName(tempname)) {
				// to count raw number of events of individual files
				f = TFile::Open(tempname);
				ttemp = (TTree*)f->Get(treename);
				// Chaining file
				cout << "\t- File: " << tempname << " - " << t2->Add(tempname) << " - " << ttemp->GetEntries() << " entries" << endl;
				// Preparing next split
				tempname = EXinFolder + EXtype + Form("%d.%d",i,++split) + ".root";
			}
		}
	}

	// Loading run3
	cout << "**************" << endl << "*** RUN 3 ***" << endl << "*************" << endl;
	TChain *t3 = new TChain(treename);
	for(Int_t i=run3a; i<=run3b; i++) {
		found=0;
		for(Int_t j=0;j<(sizeof(NotAnalyze)/sizeof(NotAnalyze[0])) && !found;j++)
			if(i==NotAnalyze[j]) found = 1;
		if(!found) {
			tempname = EXinFolder + EXtype + Form("%d",i) + ".root";
			if(!gSystem->AccessPathName(tempname)) {
				// to count raw number of events of individual files
				f = TFile::Open(tempname);
				ttemp = (TTree*)f->Get(treename);
				// Chaining file
				cout << "\t- File: " << tempname << " - " << t3->Add(tempname) << " - " << ttemp->GetEntries() << " entries" << endl;
			}
			Int_t split=1;
			tempname = EXinFolder + EXtype + Form("%d.%d",i,split) + ".root";
			while(!gSystem->AccessPathName(tempname)) {
				// to count raw number of events of individual files
				f = TFile::Open(tempname);
				ttemp = (TTree*)f->Get(treename);
				// Chaining file
				cout << "\t- File: " << tempname << " - " << t3->Add(tempname) << " - " << ttemp->GetEntries() << " entries" << endl;
				// Preparing next split
				tempname = EXinFolder + EXtype + Form("%d.%d",i,++split) + ".root";
			}
		}
	}

	// Loading run4
	cout << "**************" << endl << "*** RUN 4 ***" << endl << "*************" << endl;
	TChain *t4 = new TChain(treename);
	for(Int_t i=run4a; i<=run4b; i++) {
		found=0;
		for(Int_t j=0;j<(sizeof(NotAnalyze)/sizeof(NotAnalyze[0])) && !found;j++)
			if(i==NotAnalyze[j]) found = 1;
		if(!found) {
			tempname = EXinFolder + EXtype + Form("%d",i) + ".root";
			if(!gSystem->AccessPathName(tempname)) {
				// to count raw number of events of individual files
				f = TFile::Open(tempname);
				ttemp = (TTree*)f->Get(treename);
				// Chaining file
				cout << "\t- File: " << tempname << " - " << t4->Add(tempname) << " - " << ttemp->GetEntries() << " entries" << endl;
			}
			Int_t split=1;
			tempname = EXinFolder + EXtype + Form("%d.%d",i,split) + ".root";
			while(!gSystem->AccessPathName(tempname)) {
				// to count raw number of events of individual files
				f = TFile::Open(tempname);
				ttemp = (TTree*)f->Get(treename);
				// Chaining file
				cout << "\t- File: " << tempname << " - " << t4->Add(tempname) << " - " << ttemp->GetEntries() << " entries" << endl;
				// Preparing next split
				tempname = EXinFolder + EXtype + Form("%d.%d",i,++split) + ".root";
			}
		}
	}

//	TCut cut1 = "abs(hsxptar)<0.09 && abs(hsyptar)<0.055 && abs(missmass-1.121683)<0.026 && ssshsum > 0.8 && haero_su > 1.5 && abs(cointime-0.15)<0.75";
//	TCut cut1 = "abs(hsxptar)<0.09 && abs(hsyptar)<0.055 && abs(missmass-1.121683)<0.026 && ssshsum > 0.8 && cointime<0.9 && cointime>-1.6";
        TCut cut1 = "(abs(ssdelta)<20.0) && (ssshsum > 0.7) && (scer_npe>0.5) && (abs(hsdelta)<8.5) && (abs(hsxptar)<0.09) && (abs(hsyptar)<0.055) && (hcer_npe<2.0) && (haero_su > 1.5) && (abs(missmass-1.121683)<0.026000) && (abs(cointime-0.150000)<0.750000)";

	// Line 1
//	Double_t WQ2cutW[] = {2.322, 2.267, 2.33, 2.422};//{2.322, 2.422, 2.33, 2.267}; //y
//        Double_t WQ2cutQ2[] = {2.02, 2.362, 2.218, 1.75};//{2.02, 1.75, 2.218, 2.362}; //x
	Double_t WQ2cutW[] = {2.322, 2.422, 2.33, 2.267}; //y
        //Double_t WQ2cutQ2[] = {1.25, 1.06, 1.35, 1.5}; //x
        Double_t WQ2cutQ2[] = {2.02, 1.75, 2.218, 2.362}; //x
        Double_t p0_WQ2cut[4], p1_WQ2cut[4];
        for(Int_t i=0; i<4; i++) {
                p1_WQ2cut[i] = (WQ2cutW[(i+1)%4]-WQ2cutW[i])/(WQ2cutQ2[(i+1)%4]-WQ2cutQ2[i]);
                p0_WQ2cut[i] = WQ2cutW[i] - p1_WQ2cut[i]*WQ2cutQ2[i];
        }

	// Line 2
//	Double_t WQ2cutWa[] = {2.322, 2.267, 2.33, 2.422};//{2.322, 2.422, 2.33, 2.267}; //y
//        Double_t WQ2cutQ2a[] = {2.02, 2.362, 2.218, 1.75};//{2.02, 1.75, 2.218, 2.362}; //x
	Double_t WQ2cutWa[] = {2.305, 2.37, 2.33, 2.265}; //y
        //Double_t WQ2cutQ2a[] = {2.02, 1.75, 2.218, 2.362}; //x
        Double_t WQ2cutQ2a[] = {1.275, 1.15, 1.37, 1.59}; //x
        Double_t p0_WQ2cuta[4], p1_WQ2cuta[4];
        for(Int_t i=0; i<4; i++) {
                p1_WQ2cuta[i] = (WQ2cutWa[(i+1)%4]-WQ2cutWa[i])/(WQ2cutQ2a[(i+1)%4]-WQ2cutQ2a[i]);
                p0_WQ2cuta[i] = WQ2cutWa[i] - p1_WQ2cuta[i]*WQ2cutQ2a[i];
        }

	gStyle->SetOptStat(0);

	// Plot Q2 vs eps
	TCanvas *c0 = new TCanvas("c0","c0",600,600);
	c0->SetTicks(1,1);
	// run2
	Double_t lims0[] = {0.5,3,2.1,2.7};
	TH2F *h0a = new TH2F("h0a",";Q^{2} (GeV^{2});W (GeV})",100,lims0[0],lims0[1],100,lims0[2],lims0[3]);//,100,2,2.8,100,1.4,3);
	h0a->SetMarkerColor(2);
	t2->Draw("W:Q2>>h0a",cut1);//,"same");
	// run1
	TH2F *h0b = new TH2F("h0b","",100,lims0[0],lims0[1],100,lims0[2],lims0[3]);
	h0b->SetMarkerColor(4);
	t1->Draw("W:Q2>>h0b",cut1,"same");
	// run4
	TH2F *h0c = new TH2F("h0c","",100,lims0[0],lims0[1],100,lims0[2],lims0[3]);
	h0c->SetMarkerColor(3);
	t4->Draw("W:Q2>>h0c",cut1,"same");
	// run3
	TH2F *h0d = new TH2F("h0d","",100,lims0[0],lims0[1],100,lims0[2],lims0[3]);
	h0d->SetMarkerColor(1);
	t3->Draw("W:Q2>>h0d",cut1,"same");
	// Drawing lines
        TLine *l2aab[4];
        for(Int_t i=0;i<4;i++) {
                l2aab[i] = new TLine(WQ2cutQ2[i],WQ2cutW[i],WQ2cutQ2[(i+1)%4],WQ2cutW[(i+1)%4]);
                l2aab[i]->SetLineColor(1);
                l2aab[i]->SetLineWidth(2);
                l2aab[i]->Draw();
        }
	// Drawing lines 2
        TLine *l2aab2a[4];
        for(Int_t i=0;i<4;i++) {
                l2aab2a[i] = new TLine(WQ2cutQ2a[i],WQ2cutWa[i],WQ2cutQ2a[(i+1)%4],WQ2cutWa[(i+1)%4]);
                l2aab2a[i]->SetLineColor(6);
                l2aab2a[i]->SetLineWidth(2);
                l2aab2a[i]->Draw();
        }
	TLatex l;
	l.SetTextSize(0.03);
	l.DrawLatex(lims0[0]+0.6*(lims0[1]-lims0[0]),lims0[2]+0.90*(lims0[3]-lims0[2]),"#color[4]{Q^{2}=2.1 GeV^{2} #epsilon = 0.27}");
	l.DrawLatex(lims0[0]+0.6*(lims0[1]-lims0[0]),lims0[2]+0.85*(lims0[3]-lims0[2]),"#color[2]{Q^{2}=2.1 GeV^{2} #epsilon = 0.54}");
	l.DrawLatex(lims0[0]+0.6*(lims0[1]-lims0[0]),lims0[2]+0.80*(lims0[3]-lims0[2]),"#color[3]{Q^{2}=1.4 GeV^{2} #epsilon = 0.33}");
	l.DrawLatex(lims0[0]+0.6*(lims0[1]-lims0[0]),lims0[2]+0.75*(lims0[3]-lims0[2]),"#color[1]{Q^{2}=1.4 GeV^{2} #epsilon = 0.58}");
	//l.DrawLatex(2.1,1.6,"#color[2]{#epsilon = 0.54}");
	c0->Print("c0.png");

	return 0;
}
