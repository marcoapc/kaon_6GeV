int teps(Int_t run1a=47056, Int_t run1b=47154, Int_t run2a=47358, Int_t run2b=47428, Int_t run3a=47587, Int_t run3b=47613, Int_t run4a=47687, Int_t run4b=47736) {

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
//	TCut cut1 = "abs(hsxptar)<0.09 && abs(hsyptar)<0.055 && abs(missmass-1.121683)<0.026 && ssshsum > 0.8 && abs(cointime-0.15)<0.75";
	TCut cut1 = "((((((((((((Q2>(-2.700000*W+8.289400)) && (Q2<(-5.086957*W+14.070609)) && (Q2<(-2.285714*W+7.543714)) && (Q2>(-6.218182*W+16.458618)))&&(abs(ssdelta)<20.0))&&(ssshsum > 0.7))&&(scer_npe>0.5))&&(abs(hsdelta)<8.5))&&(abs(hsxptar)<0.09))&&(abs(hsyptar)<0.055))&&(hcer_npe<2.0))&&(haero_su > 1.5))&&(abs(missmass-1.121683)<0.026000))&&(abs(cointime-0.150000)<0.750000))"; // FPI-2 - high Q2
	TCut cut2 = "((((((((((((Q2>(-1.923077*W+5.707692)) && (Q2<(-5.500000*W+14.185000)) && (Q2<(-3.384615*W+9.256154)) && (Q2>(-7.875000*W+19.426875)))&&(abs(ssdelta)<20.0))&&(ssshsum > 0.7))&&(scer_npe>0.5))&&(abs(hsdelta)<8.5))&&(abs(hsxptar)<0.09))&&(abs(hsyptar)<0.055))&&(hcer_npe<2.0))&&(haero_su > 1.5))&&(abs(missmass-1.121683)<0.026000))&&(abs(cointime-0.150000)<0.750000))"; // FPI-2 - low Q2

	gStyle->SetOptStat(0);

	// Plot t dist FPI-2 highQ2
	TCanvas *c0 = new TCanvas("c0","c0",600,600);
	c0->SetTicks(1,1);
	TH1F *h0a = new TH1F("h0a",";t (GeV^{2});Normalized counts",100,0,1);
	h0a->SetLineColor(2);
	t2->Draw("t>>h0a",cut1);
	h0a->Scale(1.0/h0a->GetBinContent(h0a->GetMaximumBin()));
	TH1F *h0b = new TH1F("h0b","",100,0,1);
	h0b->SetLineColor(4);
	t1->Draw("t>>h0b",cut1,"same");
	h0b->Scale(1.0/h0b->GetBinContent(h0b->GetMaximumBin()));
	TLatex l;
	l.SetTextSize(0.05);
	c0->Update();
	l.DrawLatex(gPad->GetUxmin()+0.8*(gPad->GetUxmax()-gPad->GetUxmin()),gPad->GetUymin()+0.9*(gPad->GetUymax()-gPad->GetUymin()),"#color[4]{#epsilon = 0.27}");
	l.DrawLatex(gPad->GetUxmin()+0.8*(gPad->GetUxmax()-gPad->GetUxmin()),gPad->GetUymin()+0.85*(gPad->GetUymax()-gPad->GetUymin()),"#color[2]{#epsilon = 0.54}");
	c0->Print("FPI2_highQ2_t_eps.png");

	// Plot t dist FPI-2 lowQ2
	TCanvas *c1 = new TCanvas("c1","c1",600,600);
	c1->SetTicks(1,1);
	TH1F *h1a = new TH1F("h1a",";t (GeV^{2});Normalized counts",100,0,1);
	h1a->SetLineColor(2);
	t4->Draw("t>>h1a",cut2);
	h1a->Scale(1.0/h1a->GetBinContent(h1a->GetMaximumBin()));
	TH1F *h1b = new TH1F("h1b","",100,0,1);
	h1b->SetLineColor(4);
	t3->Draw("t>>h1b",cut2,"same");
	h1b->Scale(1.0/h1b->GetBinContent(h1b->GetMaximumBin()));
	TLatex l;
	l.SetTextSize(0.05);
	c1->Update();
	l.DrawLatex(gPad->GetUxmin()+0.8*(gPad->GetUxmax()-gPad->GetUxmin()),gPad->GetUymin()+0.9*(gPad->GetUymax()-gPad->GetUymin()),"#color[4]{#epsilon = 0.27}");
	l.DrawLatex(gPad->GetUxmin()+0.8*(gPad->GetUxmax()-gPad->GetUxmin()),gPad->GetUymin()+0.85*(gPad->GetUymax()-gPad->GetUymin()),"#color[2]{#epsilon = 0.54}");
	c1->Print("FPI2_lowQ2_t_eps.png");

	return 0;
}
