void plotIndiFitResults()
{
	double temp;

	int netLambda_intvl_exp, netProton_intvl_exp, netCharge_intvl_exp, netKaon_intvl_exp;

	double *netLambda_energy_exp, *netProton_energy_exp, *netCharge_energy_exp, *netKaon_energy_exp;

	double *netLambda_sigma2byM_exp_0to5, *netLambda_sigma2byM_exp_0to5_staterr, *netLambda_sigma2byM_exp_0to5_syserr;

	double *netProton_sigma2byM_exp_0to5, *netProton_sigma2byM_exp_0to5_staterr, *netProton_sigma2byM_exp_0to5_syserr;

	double *netCharge_sigma2byM_exp_0to5, *netCharge_sigma2byM_exp_0to5_staterr, *netCharge_sigma2byM_exp_0to5_syserr;

	double *netKaon_sigma2byM_exp_0to5, *netKaon_sigma2byM_exp_0to5_staterr, *netKaon_sigma2byM_exp_0to5_syserr;

	ifstream infile; string fname;

	fname = "data/cumulantRatios_netLambda.dat";
	infile.open(fname);
	if (infile)
	{
		infile >> netLambda_intvl_exp;

		netLambda_energy_exp = new double[netLambda_intvl_exp];

		netLambda_sigma2byM_exp_0to5 = new double[netLambda_intvl_exp];
		netLambda_sigma2byM_exp_0to5_syserr = new double[netLambda_intvl_exp]; netLambda_sigma2byM_exp_0to5_staterr = new double[netLambda_intvl_exp];

		for (int i = 0; i < netLambda_intvl_exp; ++i)
		{
			infile >> netLambda_energy_exp[i] >> netLambda_sigma2byM_exp_0to5[i] >> netLambda_sigma2byM_exp_0to5_staterr[i] >> netLambda_sigma2byM_exp_0to5_syserr[i]
			 >> temp >> temp >> temp
			 >> temp >> temp >> temp
			 >> temp >> temp >> temp;

			 netLambda_sigma2byM_exp_0to5_staterr[i] = sqrt(pow(netLambda_sigma2byM_exp_0to5_staterr[i],2)+pow(netLambda_sigma2byM_exp_0to5_syserr[i],2));
		}
	}
	else{
		cout << "Unable to load " << fname << "...\n";
		netLambda_intvl_exp = 0;
	}
	infile.close();

	fname = "data/cumulantRatios_netProton.dat";
	infile.open(fname);
	if (infile)
	{
		infile >> netProton_intvl_exp;

		netProton_intvl_exp--;	// remove 7.7 GeV
		infile >> temp >> temp >> temp >> temp
			 >> temp >> temp >> temp >> temp >> temp >> temp
			 >> temp >> temp >> temp
			 >> temp >> temp >> temp >> temp >> temp >> temp;

		netProton_intvl_exp--;	// remove 11 GeV
		infile >> temp >> temp >> temp >> temp
			 >> temp >> temp >> temp >> temp >> temp >> temp
			 >> temp >> temp >> temp
			 >> temp >> temp >> temp >> temp >> temp >> temp;

		netProton_energy_exp = new double[netProton_intvl_exp];

		netProton_sigma2byM_exp_0to5 = new double[netProton_intvl_exp];
		netProton_sigma2byM_exp_0to5_syserr = new double[netProton_intvl_exp]; netProton_sigma2byM_exp_0to5_staterr = new double[netProton_intvl_exp];

		for (int i = 0; i < netProton_intvl_exp; ++i)
		{
			infile >> netProton_energy_exp[i] >> netProton_sigma2byM_exp_0to5[i] >> netProton_sigma2byM_exp_0to5_staterr[i] >> netProton_sigma2byM_exp_0to5_syserr[i]
			 >> temp >> temp >> temp >> temp >> temp >> temp
			 >> temp >> temp >> temp
			 >> temp >> temp >> temp >> temp >> temp >> temp;

			 netProton_sigma2byM_exp_0to5_staterr[i] = sqrt(pow(netProton_sigma2byM_exp_0to5_staterr[i],2)+pow(netProton_sigma2byM_exp_0to5_syserr[i],2));
		}
	}
	else{
		cout << "Unable to load " << fname << "...\n";
		netProton_intvl_exp = 0;
	}
	infile.close();

	fname = "data/cumulantRatios_netKaon.dat";
	infile.open(fname);
	if (infile)
	{
		infile >> netKaon_intvl_exp;

		netKaon_energy_exp = new double[netKaon_intvl_exp];

		netKaon_sigma2byM_exp_0to5 = new double[netKaon_intvl_exp];
		netKaon_sigma2byM_exp_0to5_syserr = new double[netKaon_intvl_exp]; netKaon_sigma2byM_exp_0to5_staterr = new double[netKaon_intvl_exp];

		for (int i = 0; i < netKaon_intvl_exp; ++i)
		{
			infile >> netKaon_energy_exp[i] >> netKaon_sigma2byM_exp_0to5[i] >> netKaon_sigma2byM_exp_0to5_staterr[i] >> netKaon_sigma2byM_exp_0to5_syserr[i]
			 >> temp >> temp >> temp >> temp >> temp >> temp
			 >> temp >> temp >> temp
			 >> temp >> temp >> temp >> temp >> temp >> temp;

			 netKaon_sigma2byM_exp_0to5_staterr[i] = sqrt(pow(netKaon_sigma2byM_exp_0to5_staterr[i],2)+pow(netKaon_sigma2byM_exp_0to5_syserr[i],2));
		}
	}
	else{
		cout << "Unable to load " << fname << "...\n";
		netKaon_intvl_exp = 0;
	}
	infile.close();

	fname = "data/momentRatios_netCharge.dat";
	infile.open(fname);
	if (infile)
	{
		infile >> netCharge_intvl_exp;
		netCharge_intvl_exp--;	// remove 7.7 GeV
		infile >> temp >> temp >> temp >> temp
			 >> temp >> temp >> temp >> temp >> temp >> temp
			 >> temp >> temp >> temp
			 >> temp >> temp >> temp >> temp >> temp >> temp;

		netCharge_energy_exp = new double[netCharge_intvl_exp];

		netCharge_sigma2byM_exp_0to5 = new double[netCharge_intvl_exp];
		netCharge_sigma2byM_exp_0to5_syserr = new double[netCharge_intvl_exp]; netCharge_sigma2byM_exp_0to5_staterr = new double[netCharge_intvl_exp];

		for (int i = 0; i < netCharge_intvl_exp; ++i)
		{
			infile >> netCharge_energy_exp[i] >> netCharge_sigma2byM_exp_0to5[i] >> netCharge_sigma2byM_exp_0to5_staterr[i] >> netCharge_sigma2byM_exp_0to5_syserr[i]
			 >> temp >> temp >> temp >> temp >> temp >> temp
			 >> temp >> temp >> temp
			 >> temp >> temp >> temp >> temp >> temp >> temp;

			netCharge_sigma2byM_exp_0to5_staterr[i] = sqrt(pow(netCharge_sigma2byM_exp_0to5_staterr[i],2)+pow(netCharge_sigma2byM_exp_0to5_syserr[i],2));
		}
	}
	else{
		cout << "Unable to load " << fname << "...\n";
		netCharge_intvl_exp = 0;
	}
	infile.close();

	double T, mub, muq, mus;

	int intvl_fo_P, intvl_fo_K, intvl_fo_L;

	double *energy_fo_P, *T_fo_P, *mub_fo_P, *muq_fo_P, *mus_fo_P; 
	double *netProton_sigma2byM_fo_P, *netKaon_sigma2byM_fo_P, *netLambda_sigma2byM_fo_P;

	double *energy_fo_K, *T_fo_K, *mub_fo_K, *muq_fo_K, *mus_fo_K; 
	double *netProton_sigma2byM_fo_K, *netKaon_sigma2byM_fo_K, *netLambda_sigma2byM_fo_K;

	double *energy_fo_L, *T_fo_L, *mub_fo_L, *muq_fo_L, *mus_fo_L; 
	double *netProton_sigma2byM_fo_L, *netKaon_sigma2byM_fo_L, *netLambda_sigma2byM_fo_L;

	DataFileReader infile_self;

	
	STAR_netProton_sigma2byM netProton_sigma2byM_hrg;
	STAR_netKaon_sigma2byM netKaon_sigma2byM_hrg;
	STAR_netLambda_sigma2byM netLambda_sigma2byM_hrg;

	infile_self.read("data/freezeout_netProton_fluc.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_fo_P = infile_self.nRows();
		energy_fo_P = infile_self.getColumn(0);
		T_fo_P = infile_self.getColumn(1);
		mub_fo_P = infile_self.getColumn(3);
		muq_fo_P = infile_self.getColumn(5);
		mus_fo_P = infile_self.getColumn(7);
		infile_self.close();
		
		netProton_sigma2byM_fo_P = new double[intvl_fo_P];
		netKaon_sigma2byM_fo_P = new double[intvl_fo_P];
		netLambda_sigma2byM_fo_P = new double[intvl_fo_P];

		for (int i = 0; i < intvl_fo_P; ++i)
		{
			T = T_fo_P[i];
			mub = mub_fo_P[i];
			muq = muq_fo_P[i];
			mus = mus_fo_P[i];

			netProton_sigma2byM_fo_P[i] = netProton_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netKaon_sigma2byM_fo_P[i] = netKaon_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netLambda_sigma2byM_fo_P[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
		}
	}

	infile_self.read("data/freezeout_netKaon_fluc.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_fo_K = infile_self.nRows();
		energy_fo_K = infile_self.getColumn(0);
		T_fo_K = infile_self.getColumn(1);
		mub_fo_K = infile_self.getColumn(3);
		muq_fo_K = infile_self.getColumn(5);
		mus_fo_K = infile_self.getColumn(7);
		infile_self.close();
		
		netProton_sigma2byM_fo_K = new double[intvl_fo_K];
		netKaon_sigma2byM_fo_K = new double[intvl_fo_K];
		netLambda_sigma2byM_fo_K = new double[intvl_fo_K];

		for (int i = 0; i < intvl_fo_K; ++i)
		{
			T = T_fo_K[i];
			mub = mub_fo_K[i];
			muq = muq_fo_K[i];
			mus = mus_fo_K[i];

			netProton_sigma2byM_fo_K[i] = netProton_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netKaon_sigma2byM_fo_K[i] = netKaon_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netLambda_sigma2byM_fo_K[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
		}
	}

	infile_self.read("data/freezeout_netLambda_fluc.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_fo_L = infile_self.nRows();
		energy_fo_L = infile_self.getColumn(0);
		T_fo_L = infile_self.getColumn(1);
		mub_fo_L = infile_self.getColumn(3);
		muq_fo_L = infile_self.getColumn(5);
		mus_fo_L = infile_self.getColumn(7);
		infile_self.close();
		
		netProton_sigma2byM_fo_L = new double[intvl_fo_L];
		netKaon_sigma2byM_fo_L = new double[intvl_fo_L];
		netLambda_sigma2byM_fo_L = new double[intvl_fo_L];

		for (int i = 0; i < intvl_fo_K; ++i)
		{
			T = T_fo_L[i];
			mub = mub_fo_L[i];
			muq = muq_fo_L[i];
			mus = mus_fo_L[i];

			netProton_sigma2byM_fo_L[i] = netProton_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netKaon_sigma2byM_fo_L[i] = netKaon_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netLambda_sigma2byM_fo_L[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
		}
	}

	double topMargin = 0.01;
	double bottomMargin = 0.1;
	double rightMargin = 0.01;
	double leftMargin = 0.1;

	double textSize = 0.05;
	double markerSize = 1.6;
	double lineWidth = 1;
	double opaque = 1;
	double transparent = 0.3;

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetEndErrorSize(3);

	int Color_P = TColor::GetColor("#e31e1e");
	int Color_K = TColor::GetColor("#11a600");
	int Color_L = TColor::GetColor("#0835c7");

	int Marker_P_exp = kFullCircle;
	int Marker_K_exp = kFullCross;
	int Marker_L_exp = kFullSquare;

	int Marker_P_hrg = kOpenCircle;
	int Marker_K_hrg = kOpenCross;
	int Marker_L_hrg = kOpenSquare;

	double markerSizeRed = 0.63;

	TCanvas *cnvs = new TCanvas("cnvs_IndiFitResults","IndiFitResults",400,350);

	TPad *pad = new TPad("pad_IndiFitResults", "pad",0,0,1,1);
	pad->SetLogx();
	pad->SetLogy();
	pad->SetBottomMargin(bottomMargin);
	pad->SetTopMargin(topMargin);
	pad->SetLeftMargin(leftMargin);
	pad->SetRightMargin(rightMargin);
	pad->Draw();

	pad->cd();

	auto *gr111 = new TGraphErrors(netProton_intvl_exp,netProton_energy_exp,netProton_sigma2byM_exp_0to5,(double*)0,netProton_sigma2byM_exp_0to5_staterr);
	gr111->SetMarkerSize(markerSize-markerSizeRed);
	gr111->SetMarkerColorAlpha(Color_P,opaque);
	gr111->SetLineColorAlpha(Color_P,opaque);
	gr111->SetFillStyle(1001);
	gr111->SetFillColorAlpha(Color_P,transparent);
	gr111->SetLineWidth(lineWidth);
	gr111->SetMarkerStyle(Marker_P_exp);
	auto *gr112 = new TGraphErrors(netProton_intvl_exp,netProton_energy_exp,netProton_sigma2byM_exp_0to5,(double*)0,netProton_sigma2byM_exp_0to5_syserr);
	gr112->SetLineWidth(lineWidth);
	gr112->SetLineColorAlpha(Color_P,opaque);

	auto *gr113 = new TGraph(intvl_fo_K,energy_fo_K,netProton_sigma2byM_fo_K);
	gr113->SetMarkerSize(markerSize);
	gr113->SetMarkerColorAlpha(Color_P,opaque);
	gr113->SetMarkerStyle(Marker_K_hrg);
	auto *gr114 = new TGraph(intvl_fo_L,energy_fo_L,netProton_sigma2byM_fo_L);
	gr114->SetMarkerSize(markerSize);
	gr114->SetMarkerColorAlpha(Color_P,opaque);
	gr114->SetMarkerStyle(Marker_L_hrg);
	auto *gr115 = new TGraph(intvl_fo_P,energy_fo_P,netProton_sigma2byM_fo_P);
	gr115->SetMarkerSize(markerSize);
	gr115->SetMarkerColorAlpha(Color_P,opaque);
	gr115->SetMarkerStyle(Marker_P_hrg);

	auto *gr121 = new TGraphErrors(netKaon_intvl_exp,netKaon_energy_exp,netKaon_sigma2byM_exp_0to5,(double*)0,netKaon_sigma2byM_exp_0to5_staterr);
	gr121->SetMarkerSize(markerSize-markerSizeRed);
	gr121->SetMarkerColorAlpha(Color_K,opaque);
	gr121->SetLineColorAlpha(Color_K,opaque);
	gr121->SetFillColorAlpha(Color_K,transparent);
	gr121->SetLineWidth(lineWidth);
	gr121->SetMarkerStyle(Marker_K_exp);
	auto *gr122 = new TGraphErrors(netKaon_intvl_exp,netKaon_energy_exp,netKaon_sigma2byM_exp_0to5,(double*)0,netKaon_sigma2byM_exp_0to5_syserr);
	gr122->SetLineWidth(lineWidth);
	gr122->SetLineColorAlpha(Color_K,opaque);

	auto *gr123 = new TGraph(intvl_fo_L,energy_fo_L,netKaon_sigma2byM_fo_L);
	gr123->SetMarkerSize(markerSize);
	gr123->SetMarkerColorAlpha(Color_K,opaque);
	gr123->SetMarkerStyle(Marker_L_hrg);
	auto *gr124 = new TGraph(intvl_fo_P,energy_fo_P,netKaon_sigma2byM_fo_P);
	gr124->SetMarkerSize(markerSize);
	gr124->SetMarkerColorAlpha(Color_K,opaque);
	gr124->SetMarkerStyle(Marker_P_hrg);
	auto *gr125 = new TGraph(intvl_fo_K,energy_fo_K,netKaon_sigma2byM_fo_K);
	gr124->SetMarkerSize(markerSize);
	gr124->SetMarkerColorAlpha(Color_K,opaque);
	gr124->SetMarkerStyle(Marker_K_hrg);

	auto *gr131 = new TGraphErrors(netLambda_intvl_exp,netLambda_energy_exp,netLambda_sigma2byM_exp_0to5,(double*)0,netLambda_sigma2byM_exp_0to5_staterr);
	gr131->SetMarkerSize(markerSize-markerSizeRed);
	gr131->SetMarkerColorAlpha(Color_L,opaque);
	gr131->SetLineColorAlpha(Color_L,opaque);
	gr131->SetFillColorAlpha(Color_L,transparent);
	gr131->SetLineWidth(lineWidth);
	gr131->SetMarkerStyle(Marker_L_exp);
	auto *gr132 = new TGraphErrors(netLambda_intvl_exp,netLambda_energy_exp,netLambda_sigma2byM_exp_0to5,(double*)0,netLambda_sigma2byM_exp_0to5_syserr);
	gr132->SetLineWidth(lineWidth);
	gr132->SetLineColorAlpha(Color_L,opaque);

	auto *gr133 = new TGraph(intvl_fo_P,energy_fo_P,netLambda_sigma2byM_fo_P);
	gr133->SetMarkerSize(markerSize);
	gr133->SetMarkerColorAlpha(Color_L,opaque);
	gr133->SetMarkerStyle(Marker_P_hrg);
	auto *gr134 = new TGraph(intvl_fo_K,energy_fo_K,netLambda_sigma2byM_fo_K);
	gr134->SetMarkerSize(markerSize);
	gr134->SetMarkerColorAlpha(Color_L,opaque);
	gr134->SetMarkerStyle(Marker_K_hrg);
	auto *gr135 = new TGraph(intvl_fo_L,energy_fo_L,netLambda_sigma2byM_fo_L);
	gr135->SetMarkerSize(markerSize);
	gr135->SetMarkerColorAlpha(Color_L,opaque);
	gr135->SetMarkerStyle(Marker_L_hrg);

	TMultiGraph *mg1  = new TMultiGraph();

	mg1->Add(gr111,"pl3");
	// mg1->Add(gr112,"||");
	mg1->Add(gr113,"p");
	mg1->Add(gr114,"p");
	mg1->Add(gr115,"p");
	
	mg1->Add(gr121,"pl3");
	// mg1->Add(gr122,"||");
	mg1->Add(gr123,"p");
	mg1->Add(gr124,"p");
	mg1->Add(gr125,"p");

	mg1->Add(gr131,"pl3");
	// mg1->Add(gr132,"||");
	mg1->Add(gr133,"p");
	mg1->Add(gr134,"p");
	mg1->Add(gr135,"p");

	mg1->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg1->GetYaxis()->SetTitle("#sigma^{2}/M");
	mg1->GetXaxis()->SetLimits(15,250);
	mg1->GetYaxis()->SetRangeUser(0.7,80);
	// mg1->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg1->GetXaxis()->CenterTitle(true);
	mg1->GetYaxis()->CenterTitle(true);
	mg1->GetXaxis()->SetNoExponent();
	mg1->GetXaxis()->SetLabelSize(textSize);
	mg1->GetXaxis()->SetTitleSize(textSize);
	mg1->GetXaxis()->SetTitleOffset(0.75);
	mg1->GetYaxis()->SetLabelSize(textSize);
	mg1->GetYaxis()->SetTitleSize(textSize);
	mg1->GetYaxis()->SetTickLength(0.02);
	mg1->GetYaxis()->SetTitleOffset(0.95);
	mg1->Draw("a");

	auto *grP = new TGraph();
	grP->SetMarkerSize(markerSize);
	grP->SetMarkerColorAlpha(kBlack,opaque);
	grP->SetMarkerStyle(Marker_P_hrg);

	auto *grK = new TGraph();
	grK->SetMarkerSize(markerSize);
	grK->SetMarkerColorAlpha(kBlack,opaque);
	grK->SetMarkerStyle(Marker_K_hrg);

	auto *grL = new TGraph();
	grL->SetMarkerSize(markerSize);
	grL->SetMarkerColorAlpha(kBlack,opaque);
	grL->SetMarkerStyle(Marker_L_hrg);

	TLegend *lgnd1 = new TLegend(0.15,0.63,0.55,0.85);
	lgnd1->SetNColumns(2);
	lgnd1->SetHeader("HRG Fits(data used)","C");
	lgnd1->AddEntry(grP,"net-p","p");
	lgnd1->AddEntry(grK,"net-k","p");
	lgnd1->AddEntry(grL,"net-#Lambda","p");
	lgnd1->SetTextSize(textSize);
	lgnd1->SetBorderSize(0);
	lgnd1->SetFillStyle(0);
	lgnd1->Draw();

	TLegend *lgnd2 = new TLegend(0.6,0.12,0.95,0.35);
	lgnd2->SetNColumns(2);
	lgnd2->SetHeader("STAR Au-Au 0-5%","C");
	lgnd2->AddEntry(gr111,"net-p","p");
	lgnd2->AddEntry(gr121,"net-k","p");
	lgnd2->AddEntry(gr131,"net-#Lambda","p");
	lgnd2->SetTextSize(textSize);
	lgnd2->SetBorderSize(0);
	lgnd2->SetFillStyle(0);
	lgnd2->Draw();

	cnvs->SaveAs("plots/cumulantRatios_IndiFitResults.pdf");
}