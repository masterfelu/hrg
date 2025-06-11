#include <iostream>
#include <cmath>
#include <fstream>

#include "Particles.h"
#include "ThermalFunctions.h"

#include "TROOT.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TStyle.h"

using namespace std;

class DataFileReader
{
	private:

	int n_rows {0}, n_columns {0};
	double **data {nullptr};
	bool dataLoaded {false};
	std::string fileName {};

	public:

	DataFileReader(){}

	~DataFileReader()
	{
		for (int i = 0; i < n_rows; ++i)
		{
			delete[] data[i];
		}
		delete[] data;
		data = nullptr;
		n_rows = 0; n_columns = 0;
	}

	int nRows()			{ return n_rows; }
	int nColumns()		{ return n_columns; }
	bool isDataLoaded()	{ return dataLoaded; }

	void read(std::string fname, bool showData = false)
	{
		fileName = fname;
		std::ifstream infile;
		infile.open(fileName);
		if(infile)
		{
			infile >> n_rows >> n_columns;
			std::cout << "Reading data from " << fileName << " with " << n_rows << " rows and " << n_columns << " columns..." << '\n';
			data = new double*[n_rows];
			for (int i = 0; i < n_rows; ++i)
			{
				data[i] = new double[n_columns];
			}
			for (int i = 0; i < n_rows; ++i)
			{
				for (int j = 0; j < n_columns; ++j)
				{
					infile >> data[i][j];
					if (showData)
					{
						std::cout << data[i][j] << '\t';
					}
				}
				if (showData)
				{
					std::cout << '\n';
				}
			}
			dataLoaded = true;
			std::cout << "Data loaded successfully..." << std::endl;
		}
		else
		{
			std::cout << "Unable to load data file " << fname << "..." << std::endl;
		}
		infile.close();
	}

	void close()
	{
		dataLoaded = false;

		for (int i = 0; i < n_rows; ++i)
		{
			delete[] data[i];
		}
		delete[] data;
		data = nullptr;
		n_rows = 0; n_columns = 0;
	}

	double get(int r, int c)
	{
		if (dataLoaded)
		{
			if(r >= 0 && r < n_rows && c >= 0 && c < n_columns)
			{
				return data[r][c];
			}
		}
		return 0;
	}

	double* getColumn(int c)
	{
		if (dataLoaded)
		{
			if(c >= 0 && c < n_columns)
			{
				double *col_data = new double[n_rows];
				for (int i = 0; i < n_rows; ++i)
				{
					col_data[i] = data[i][c];
				}
				return col_data;
			}
		}
		return 0;
	}

	double* getRow(int r)
	{
		if (dataLoaded)
		{
			if(r >= 0 && r < n_rows)
			{
				double *row_data = new double[n_columns];
				for (int i = 0; i < n_columns; ++i)
				{
					row_data[i] = data[r][i];
				}
				return row_data;
			}
		}
		return 0;
	}
};

class STAR_netLambda_sigma2byM:public ThermalFunction
{

	Susceptibility1S *chi1S;
	Susceptibility2S *chi2S;

	public:

	STAR_netLambda_sigma2byM()
	{
		ParticleSystem system;
		system.loadDefaultData();

		chi1S = new Susceptibility1S {system};		
		chi1S->setDetectorRapidityCoordinates();
		chi1S->setIntegralLimits(0.9,2,-0.5,0.5);
		chi1S->setDecay(true);

		chi2S = new Susceptibility2S {system};
		chi2S->setDetectorRapidityCoordinates();
		chi2S->setIntegralLimits(0.9,2,-0.5,0.5);
		chi2S->setDecay(true);
	}
	~STAR_netLambda_sigma2byM()
	{
		delete chi1S;
		delete chi2S;
		chi1S = nullptr;
		chi2S = nullptr;
	}

	double getValue(double T, double mub, double muq, double mus)
	{
		double num = chi2S->getValueEach(33,T,mub,muq,mus)+chi2S->getValueEach(34,T,mub,muq,mus);
		double denom = -(chi1S->getValueEach(33,T,mub,muq,mus)+chi1S->getValueEach(34,T,mub,muq,mus));
		return num/denom;
	}
};

class STAR_netLambda_Ssigma:public ThermalFunction
{

	Susceptibility3S *chi3S;
	Susceptibility2S *chi2S;

	public:

	STAR_netLambda_Ssigma()
	{
		ParticleSystem system;
		system.loadDefaultData();

		chi3S = new Susceptibility3S {system};		
		chi3S->setDetectorRapidityCoordinates();
		chi3S->setIntegralLimits(0.9,2,-0.5,0.5);
		chi3S->setDecay(true);

		chi2S = new Susceptibility2S {system};
		chi2S->setDetectorRapidityCoordinates();
		chi2S->setIntegralLimits(0.9,2,-0.5,0.5);
		chi2S->setDecay(true);
	}
	~STAR_netLambda_Ssigma()
	{
		delete chi3S;
		delete chi2S;
		chi3S = nullptr;
		chi2S = nullptr;
	}

	double getValue(double T, double mub, double muq, double mus)
	{
		double denom = chi2S->getValueEach(33,T,mub,muq,mus)+chi2S->getValueEach(34,T,mub,muq,mus);
		double num = -(chi3S->getValueEach(33,T,mub,muq,mus)+chi3S->getValueEach(34,T,mub,muq,mus));
		return num/denom;
	}
};

class STAR_netProton_sigma2byM:public ThermalFunction
{

	Susceptibility1B *chi1B;
	Susceptibility2B *chi2B;

	public:

	STAR_netProton_sigma2byM()
	{
		ParticleSystem system;
		system.loadDefaultData();

		chi1B = new Susceptibility1B {system};		
		chi1B->setDetectorRapidityCoordinates();
		chi1B->setIntegralLimits(0.4,0.8,-0.5,0.5);
		chi1B->setDecay(true);

		chi2B = new Susceptibility2B {system};
		chi2B->setDecay(true);
		chi2B->setIntegralLimits(0.4,0.8,-0.5,0.5);
		chi2B->setDetectorRapidityCoordinates();
	}
	~STAR_netProton_sigma2byM()
	{
		delete chi1B;
		delete chi2B;
		chi1B = nullptr;
		chi2B = nullptr;
	}

	double getValue(double T, double mub, double muq, double mus)
	{
		double num = chi2B->getValueEach(23,T,mub,muq,mus)+chi2B->getValueEach(24,T,mub,muq,mus);
		double denom = chi1B->getValueEach(23,T,mub,muq,mus)+chi1B->getValueEach(24,T,mub,muq,mus);
		return num/denom;
	}
};

class STAR_netCharge_sigma2byM:public ThermalFunction
{

	Susceptibility1Q *chi1Q;
	Susceptibility2Q *chi2Q;

	public:

	STAR_netCharge_sigma2byM()
	{
		ParticleSystem system;
		system.loadDefaultData();

		chi1Q = new Susceptibility1Q {system};		
		chi1Q->setDetectorPseudorapidityCoordinates();
		chi1Q->setIntegralLimits(0.2,2,-0.5,0.5);
		chi1Q->setDecay(true);

		chi2Q = new Susceptibility2Q {system};
		chi2Q->setDetectorPseudorapidityCoordinates();
		chi2Q->setIntegralLimits(0.2,2,-0.5,0.5);
		chi2Q->setDecay(true);
	}
	~STAR_netCharge_sigma2byM()
	{
		delete chi1Q;
		delete chi2Q;
		chi1Q = nullptr;
		chi2Q = nullptr;
	}

	double getValue(double T, double mub, double muq, double mus)
	{
		double num = chi2Q->getValue(T,mub,muq,mus);
		double denom = chi1Q->getValue(T,mub,muq,mus);
		return num/denom;
	}
};

class STAR_netKaon_sigma2byM:public ThermalFunction
{

	Susceptibility1S *chi1S;
	Susceptibility2S *chi2S;

	public:

	STAR_netKaon_sigma2byM()
	{
		ParticleSystem system;
		system.loadDefaultData();

		chi1S = new Susceptibility1S {system};		
		chi1S->setDetectorRapidityCoordinates();
		chi1S->setIntegralLimits(0.2,1.6,-0.5,0.5);
		chi1S->setDecay(true);

		chi2S = new Susceptibility2S {system};
		chi2S->setDecay(true);
		chi2S->setIntegralLimits(0.2,1.6,-0.5,0.5);
		chi2S->setDetectorRapidityCoordinates();
	}
	~STAR_netKaon_sigma2byM()
	{
		delete chi1S;
		delete chi2S;
		chi1S = nullptr;
		chi2S = nullptr;
	}

	double getValue(double T, double mub, double muq, double mus)
	{
		double num = chi2S->getValueEach(4,T,mub,muq,mus)+chi2S->getValueEach(5,T,mub,muq,mus);
		double denom = chi1S->getValueEach(4,T,mub,muq,mus)+chi1S->getValueEach(5,T,mub,muq,mus);
		// std::cout << "chi2P/chi1P(" << T << ", " << mub << ", " << muq << ", " << mus << ") : " << num/denom << '\n';
		return num/denom;
	}
};

namespace MyStyle
{
	int Color1 = TColor::GetColor("#8c8c8c");
	int Color2 = TColor::GetColor("#12c415");
	int Color3 = TColor::GetColor("#ff9900");
	int Color4 = TColor::GetColor("#ff0000");
	int Color5 = TColor::GetColor("#0000ff");
	int Color6 = TColor::GetColor("#0c9fcc");
	int Color7 = TColor::GetColor("#cd36ff");

	int Marker1 = kFullCircle;
	int Marker2 = kFullSquare;
	int Marker3 = kFullTriangleUp;	
	int Marker4 = kFullTriangleDown;

	int Line1 = 1;
	int Line2 = 2;
	int Line3 = 9;
};

namespace MyStyle1
{
	int Color1 = TColor::GetColor("#ff0000");
	int Color2 = TColor::GetColor("#0000ff");
	int Color3 = TColor::GetColor("#0abf0a");
	int Color4 = TColor::GetColor("#ba0bb4");

	int Marker1 = kFullCircle;

	int Marker_P = kFullCircle;
	int Marker_L = kFullCircle;	
	int Marker_K = kFullCircle;
	int Marker_Q = kFullCircle;

	int Marker_PQ = kOpenSquare;
	int Marker_PL = kOpenCircle;	
	int Marker_KQ = kOpenCircle;
	int Marker_KL = kOpenSquare;

	int Line1 = 2;
};

void plotNetLambdaResults()
{
	double temp;

	int intv_exp;

	double *energy_exp;

	double *netLambda_sigma2byM_exp_0to5, *netLambda_sigma2byM_exp_0to5_staterr, *netLambda_sigma2byM_exp_0to5_syserr;
	double *netLambda_Ssigma_exp_0to5, *netLambda_Ssigma_exp_0to5_syserr, *netLambda_Ssigma_exp_0to5_staterr;

	double *netLambda_sigma2byM_exp_70to80, *netLambda_sigma2byM_exp_70to80_staterr, *netLambda_sigma2byM_exp_70to80_syserr;
	double *netLambda_Ssigma_exp_70to80, *netLambda_Ssigma_exp_70to80_syserr, *netLambda_Ssigma_exp_70to80_staterr;

	ifstream infile;
	infile.open("data/cumulantRatios_netLambda.dat");
	if (infile)
	{
		infile >> intv_exp;

		energy_exp = new double[intv_exp];

		netLambda_sigma2byM_exp_0to5 = new double[intv_exp]; netLambda_sigma2byM_exp_0to5_syserr = new double[intv_exp]; netLambda_sigma2byM_exp_0to5_staterr = new double[intv_exp];
		netLambda_Ssigma_exp_0to5 = new double[intv_exp]; netLambda_Ssigma_exp_0to5_syserr = new double[intv_exp]; netLambda_Ssigma_exp_0to5_staterr = new double[intv_exp];

		netLambda_sigma2byM_exp_70to80 = new double[intv_exp]; netLambda_sigma2byM_exp_70to80_syserr = new double[intv_exp]; netLambda_sigma2byM_exp_70to80_staterr = new double[intv_exp];
		netLambda_Ssigma_exp_70to80 = new double[intv_exp]; netLambda_Ssigma_exp_70to80_syserr = new double[intv_exp]; netLambda_Ssigma_exp_70to80_staterr = new double[intv_exp];

		for (int i = 0; i < intv_exp; ++i)
		{
			infile >> energy_exp[i] >> netLambda_sigma2byM_exp_0to5[i] >> netLambda_sigma2byM_exp_0to5_staterr[i] >> netLambda_sigma2byM_exp_0to5_syserr[i]
			 >> netLambda_Ssigma_exp_0to5[i] >> netLambda_Ssigma_exp_0to5_staterr[i] >> netLambda_Ssigma_exp_0to5_syserr[i]
			 >> netLambda_sigma2byM_exp_70to80[i] >> netLambda_sigma2byM_exp_70to80_staterr[i] >> netLambda_sigma2byM_exp_70to80_syserr[i]
			 >> netLambda_Ssigma_exp_70to80[i] >> netLambda_Ssigma_exp_70to80_staterr[i] >> netLambda_Ssigma_exp_70to80_syserr[i];
		}
	}
	else{
		cout << "Unable to load cumulant/moment ratios...\n";
		intv_exp = 0;
	}
	infile.close();

	double E_l = 5;
	double E_u = 250;
	int intvl = 100;
	double dE = pow(E_u/E_l,1.0/(intvl-1));

	double x;
	double T, mub, muq, mus;

	double energy_cleymans[intvl];

	double netLambda_sigma2byM_cleymans[intvl], netLambda_Ssigma_cleymans[intvl];

	STAR_netLambda_sigma2byM netLambda_sigma2byM_hrg;
	STAR_netLambda_Ssigma netLambda_Ssigma_hrg;

	for (int i = 0; i < intvl; ++i)
	{
		x = E_l*pow(dE,i);
		energy_cleymans[i] = x;

		mub = 1.308/(1+0.273*x);
		muq = -0.0202/(1+0.126*x);
		mus = 0.214/(1+0.184*x);

		T = 0.166-0.139*mub*mub-0.053*pow(mub,4);

		netLambda_sigma2byM_cleymans[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
		netLambda_Ssigma_cleymans[i] = netLambda_Ssigma_hrg.getValue(T,mub,muq,mus);
	}

	int intvl_fo_PQ, intvl_fo_KQ;

	double *energy_fo_PQ, *T_fo_PQ, *mub_fo_PQ, *muq_fo_PQ, *mus_fo_PQ; 
	double *netLambda_sigma2byM_fo_PQ, *netLambda_Ssigma_fo_PQ;

	double *energy_fo_KQ, *T_fo_KQ, *mub_fo_KQ, *muq_fo_KQ, *mus_fo_KQ; 
	double *netLambda_sigma2byM_fo_KQ, *netLambda_Ssigma_fo_KQ;

	DataFileReader infile_self;

	infile_self.read("data/freezeout_netProtonCharge_fluc.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_fo_PQ = infile_self.nRows();
		energy_fo_PQ = infile_self.getColumn(0);
		T_fo_PQ = infile_self.getColumn(1);
		mub_fo_PQ = infile_self.getColumn(3);
		muq_fo_PQ = infile_self.getColumn(5);
		mus_fo_PQ = infile_self.getColumn(7);
		infile_self.close();

		netLambda_sigma2byM_fo_PQ = new double[intvl_fo_PQ];
		netLambda_Ssigma_fo_PQ = new double[intvl_fo_PQ];
		for (int i = 0; i < intvl_fo_PQ; ++i)
		{
			T = T_fo_PQ[i];
			mub = mub_fo_PQ[i];
			muq = muq_fo_PQ[i];
			mus = mus_fo_PQ[i];

			netLambda_sigma2byM_fo_PQ[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netLambda_Ssigma_fo_PQ[i] = netLambda_Ssigma_hrg.getValue(T,mub,muq,mus);
		}
	}

	infile_self.read("data/freezeout_netKaonCharge_fluc.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_fo_KQ = infile_self.nRows();
		energy_fo_KQ = infile_self.getColumn(0);
		T_fo_KQ = infile_self.getColumn(1);
		mub_fo_KQ = infile_self.getColumn(3);
		muq_fo_KQ = infile_self.getColumn(5);
		mus_fo_KQ = infile_self.getColumn(7);
		infile_self.close();

		netLambda_sigma2byM_fo_KQ = new double[intvl_fo_KQ];
		netLambda_Ssigma_fo_KQ = new double[intvl_fo_KQ];
		for (int i = 0; i < intvl_fo_KQ; ++i)
		{
			T = T_fo_KQ[i];
			mub = mub_fo_KQ[i];
			muq = muq_fo_KQ[i];
			mus = mus_fo_KQ[i];

			netLambda_sigma2byM_fo_KQ[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netLambda_Ssigma_fo_KQ[i] = netLambda_Ssigma_hrg.getValue(T,mub,muq,mus);
		}
	}	

	int intvl_bellwied;
	double *energy_bellwied;
	double *netLambda_sigma2byM_bellwied_kaon, *netLambda_sigma2byM_bellwied_qp;
	double *netLambda_Ssigma_bellwied_kaon, *netLambda_Ssigma_bellwied_qp;

	infile_self.read("data/cumulantRatios_bellwied.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_bellwied = infile_self.nRows();
		energy_bellwied = infile_self.getColumn(0);
		netLambda_sigma2byM_bellwied_kaon = infile_self.getColumn(1);
		netLambda_sigma2byM_bellwied_qp = infile_self.getColumn(2);
		netLambda_Ssigma_bellwied_kaon = infile_self.getColumn(3);
		netLambda_Ssigma_bellwied_qp = infile_self.getColumn(4);
		infile_self.close();
	}

	double topMargin = 0.01;
	double bottomMargin = 0.07;
	double rightMargin = 0.02;
	double leftMargin = 0.1;

	double textSize = 0.03;
	double markerSize = 1.5;
	double lineWidth = 2;
	double transparency = 1;

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetEndErrorSize(5);

	using namespace MyStyle;

	TCanvas *cnvs1 = new TCanvas("cnvs1_Lambda","netLambda_fluc",400,500);

	double ySizeEach = (1 - topMargin - bottomMargin)/2;

	TPad *pad11 = new TPad("pad11", "pad11",0,ySizeEach+bottomMargin,1,1);
	pad11->SetLogx();
	pad11->SetTickx();
	pad11->SetBottomMargin(0);
	pad11->SetTopMargin(topMargin/pad11->GetHNDC());
	pad11->SetRightMargin(rightMargin);
	pad11->Draw();

	TPad *pad12 = new TPad("pad12", "pad12",0,0,1,ySizeEach+bottomMargin);
	pad12->SetLogx();
	pad12->SetTickx();
	pad12->SetTopMargin(0);
	pad12->SetBottomMargin(bottomMargin/pad12->GetHNDC());
	pad12->SetRightMargin(rightMargin);
	pad12->Draw();

	pad11->cd();

	auto *gr111 = new TGraph(intvl,energy_cleymans,netLambda_sigma2byM_cleymans);
	gr111->SetLineStyle(Line1);
	gr111->SetLineColorAlpha(Color1,transparency);
	gr111->SetLineWidth(lineWidth);

	auto *gr112 = new TGraphErrors(intv_exp,energy_exp,netLambda_sigma2byM_exp_0to5,(double*)0,netLambda_sigma2byM_exp_0to5_staterr);
	gr112->SetMarkerSize(markerSize);
	gr112->SetMarkerColorAlpha(Color2,transparency);
	gr112->SetLineColorAlpha(Color2,transparency);
	gr112->SetLineWidth(lineWidth);
	gr112->SetMarkerStyle(Marker1);
	auto *gr113 = new TGraphErrors(intv_exp,energy_exp,netLambda_sigma2byM_exp_0to5,(double*)0,netLambda_sigma2byM_exp_0to5_syserr);
	gr113->SetLineWidth(lineWidth);
	gr113->SetLineColorAlpha(Color2,transparency);

	auto *gr114 = new TGraphErrors(intv_exp,energy_exp,netLambda_sigma2byM_exp_70to80,(double*)0,netLambda_sigma2byM_exp_70to80_staterr);
	gr114->SetMarkerSize(markerSize);
	gr114->SetMarkerColorAlpha(Color3,transparency);
	gr114->SetLineColorAlpha(Color3,transparency);
	gr114->SetLineWidth(lineWidth);
	gr114->SetMarkerStyle(Marker2);
	auto *gr115 = new TGraphErrors(intv_exp,energy_exp,netLambda_sigma2byM_exp_70to80,(double*)0,netLambda_sigma2byM_exp_70to80_syserr);
	gr115->SetLineWidth(lineWidth);
	gr115->SetLineColorAlpha(Color3,transparency);

	auto *gr116 = new TGraph(intvl_fo_PQ,energy_fo_PQ,netLambda_sigma2byM_fo_PQ);
	gr116->SetMarkerSize(markerSize);
	gr116->SetMarkerColorAlpha(Color4,transparency);
	gr116->SetMarkerStyle(Marker3);

	auto *gr117 = new TGraph(intvl_fo_KQ,energy_fo_KQ,netLambda_sigma2byM_fo_KQ);
	gr117->SetMarkerSize(markerSize);
	gr117->SetMarkerColorAlpha(Color5,transparency);
	gr117->SetMarkerStyle(Marker4);

	auto *gr118 = new TGraph(intvl_bellwied,energy_bellwied,netLambda_sigma2byM_bellwied_kaon);
	gr118->SetLineStyle(Line2);
	gr118->SetLineWidth(lineWidth);
	gr118->SetLineColorAlpha(Color6,transparency);

	auto *gr119 = new TGraph(intvl_bellwied,energy_bellwied,netLambda_sigma2byM_bellwied_qp);
	gr119->SetLineStyle(Line3);
	gr119->SetLineWidth(lineWidth);
	gr119->SetLineColorAlpha(Color7,transparency);

	TMultiGraph *mg11  = new TMultiGraph();

	mg11->Add(gr111,"l");
	mg11->Add(gr112,"pz");	// 0-5% centrality
	mg11->Add(gr113,"||");	// 0-5% centrality err
	mg11->Add(gr114,"pz");	// 70-80% centraliity
	mg11->Add(gr115,"||");	// 70-80% centraliity err
	mg11->Add(gr116,"p");
	mg11->Add(gr117,"p");
	mg11->Add(gr118,"l");
	// mg11->Add(gr119,"l");	// bellwied using alba (qp kit)

	// mg11->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg11->GetYaxis()->SetTitle("#sigma^{2}/M");
	// mg11->GetXaxis()->SetRangeUser(2,250);
	mg11->GetYaxis()->SetRangeUser(0.15,13.25);
	mg11->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg11->GetXaxis()->CenterTitle(true);
	mg11->GetYaxis()->CenterTitle(true);
	mg11->GetXaxis()->SetNoExponent();
	mg11->GetXaxis()->SetLabelSize(textSize/pad11->GetHNDC());
	mg11->GetXaxis()->SetTitleSize(textSize/pad11->GetHNDC());
	mg11->GetXaxis()->SetTitleOffset(0.75);
	mg11->GetYaxis()->SetLabelSize(textSize/pad11->GetHNDC());
	mg11->GetYaxis()->SetTitleSize(textSize/pad11->GetHNDC());
	mg11->GetYaxis()->SetTickLength(0.02);
	mg11->GetYaxis()->SetTitleOffset(0.75);
	mg11->Draw("a");

	TLatex *txt11 = new TLatex(0.4,0.8,"#splitline{net-#Lambda results}{and calculations}");
	txt11->SetNDC();
	txt11->SetTextSize(textSize/pad11->GetHNDC());
	txt11->Draw();

	TLegend *lgnd11 = new TLegend(0.12,0.25,0.4,0.6);
	lgnd11->AddEntry(gr116,"HRG (net-qp fit)","p");
	lgnd11->AddEntry(gr117,"HRG (net-kq fit)","p");
	lgnd11->AddEntry(gr118,"Bellwied et al. (kaon FO)","l");
	// lgnd11->AddEntry(gr119,"Bellwied et al. (light FO)","l");
	lgnd11->SetTextSize(textSize/pad12->GetHNDC());
	lgnd11->SetBorderSize(0);
	lgnd11->SetFillStyle(0);
	lgnd11->Draw();

	pad12->cd();

	auto *gr121 = new TGraphErrors(intv_exp,energy_exp,netLambda_Ssigma_exp_0to5,(double*)0,netLambda_Ssigma_exp_0to5_staterr);
	gr121->SetMarkerSize(markerSize);
	gr121->SetMarkerColorAlpha(Color2,transparency);
	gr121->SetLineColorAlpha(Color2,transparency);
	gr121->SetLineWidth(lineWidth);
	gr121->SetMarkerStyle(Marker1);
	auto *gr122 = new TGraphErrors(intv_exp,energy_exp,netLambda_Ssigma_exp_0to5,(double*)0,netLambda_Ssigma_exp_0to5_syserr);
	gr122->SetLineWidth(lineWidth);
	gr122->SetLineColorAlpha(Color2,transparency);

	auto *gr123 = new TGraphErrors(intv_exp,energy_exp,netLambda_Ssigma_exp_70to80,(double*)0,netLambda_Ssigma_exp_70to80_staterr);
	gr123->SetMarkerSize(markerSize);
	gr123->SetMarkerColorAlpha(Color3,transparency);
	gr123->SetLineColorAlpha(Color3,transparency);
	gr123->SetLineWidth(lineWidth);
	gr123->SetMarkerStyle(Marker2);
	auto *gr124 = new TGraphErrors(intv_exp,energy_exp,netLambda_Ssigma_exp_70to80,(double*)0,netLambda_Ssigma_exp_70to80_syserr);
	gr124->SetLineWidth(lineWidth);
	gr124->SetLineColorAlpha(Color3,transparency);

	auto *gr125 = new TGraph(intvl,energy_cleymans,netLambda_Ssigma_cleymans);
	gr125->SetLineStyle(1);
	gr125->SetLineColorAlpha(Color1,transparency);
	gr125->SetLineWidth(lineWidth);
	auto *gr126 = new TGraph(intvl_fo_PQ,energy_fo_PQ,netLambda_Ssigma_fo_PQ);
	gr126->SetMarkerSize(markerSize);
	gr126->SetMarkerColorAlpha(Color4,transparency);
	gr126->SetMarkerStyle(Marker3);

	auto *gr127 = new TGraph(intvl_fo_KQ,energy_fo_KQ,netLambda_Ssigma_fo_KQ);
	gr127->SetMarkerSize(markerSize);
	gr127->SetMarkerColorAlpha(Color5,transparency);
	gr127->SetMarkerStyle(Marker4);

	auto *gr128 = new TGraph(intvl_bellwied,energy_bellwied,netLambda_Ssigma_bellwied_kaon);
	gr128->SetLineStyle(Line2);
	gr128->SetLineWidth(lineWidth);
	gr128->SetLineColorAlpha(Color6,transparency);

	auto *gr129 = new TGraph(intvl_bellwied,energy_bellwied,netLambda_Ssigma_bellwied_qp);
	gr129->SetLineStyle(Line3);
	gr129->SetLineWidth(lineWidth);
	gr129->SetLineColorAlpha(Color7,transparency);

	TMultiGraph *mg12  = new TMultiGraph();

	mg12->Add(gr121,"pz");	// 0-5% centrality
	mg12->Add(gr122,"||");	// 0-5% centrality err
	mg12->Add(gr123,"pz");	// 70-80% centraliity
	mg12->Add(gr124,"||");	// 70-80% centraliity err
	mg12->Add(gr125,"l");
	mg12->Add(gr126,"p");
	mg12->Add(gr127,"p");
	mg12->Add(gr128,"l");
	// mg12->Add(gr129,"l");	// bellwied using alba (qp fit)

	mg12->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg12->GetYaxis()->SetTitle("S#sigma");
	// mg12->GetXaxis()->SetLimits(5,250);
	mg12->GetYaxis()->SetRangeUser(-0.15,1.25);
	mg12->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg12->GetXaxis()->CenterTitle(true);
	mg12->GetYaxis()->CenterTitle(true);
	mg12->GetXaxis()->SetNoExponent();
	mg12->GetXaxis()->SetLabelSize(textSize/pad12->GetHNDC());
	mg12->GetXaxis()->SetTitleSize(textSize/pad12->GetHNDC());
	mg12->GetXaxis()->SetTitleOffset(0.75);
	mg12->GetYaxis()->SetLabelSize(textSize/pad12->GetHNDC());
	mg12->GetYaxis()->SetTitleSize(textSize/pad12->GetHNDC());
	mg12->GetYaxis()->SetTickLength(0.02);
	mg12->GetYaxis()->SetTitleOffset(0.75);
	mg12->Draw("a");

	TLatex *txt12 = new TLatex(0.6,0.8,"#splitline{0.9 < p_{T} < 2 (GeV/c)}{#left|y#right| #leq 0.5}");
	txt12->SetNDC();
	txt12->SetTextSize(textSize/pad12->GetHNDC());
	txt12->Draw();

	TLegend *lgnd12 = new TLegend(0.12,0.17,0.4,0.43);
	lgnd12->AddEntry(gr112,"STAR Au-Au 0 - 5%","p");
	lgnd12->AddEntry(gr114,"STAR Au-Au 70 - 80%","p");
	lgnd12->AddEntry(gr111,"HRG (Cleymans et al.)","l");
	lgnd12->SetTextSize(textSize/pad12->GetHNDC());
	lgnd12->SetBorderSize(0);
	lgnd12->SetFillStyle(0);
	lgnd12->Draw();

	cnvs1->SaveAs("plots/cumulantRatios_netLambda.pdf");
}

void plotFitResults_Calib()
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

	int intvl_fo_PQ, intvl_fo_KQ, intvl_fo_KL, intvl_fo_PL;

	double *energy_fo_PQ, *T_fo_PQ, *mub_fo_PQ, *muq_fo_PQ, *mus_fo_PQ; 
	double *netLambda_sigma2byM_fo_PQ, *netProton_sigma2byM_fo_PQ, *netCharge_sigma2byM_fo_PQ, *netKaon_sigma2byM_fo_PQ;

	double *energy_fo_KQ, *T_fo_KQ, *mub_fo_KQ, *muq_fo_KQ, *mus_fo_KQ; 
	double *netLambda_sigma2byM_fo_KQ, *netProton_sigma2byM_fo_KQ, *netCharge_sigma2byM_fo_KQ, *netKaon_sigma2byM_fo_KQ;

	double *energy_fo_KL, *T_fo_KL, *mub_fo_KL, *muq_fo_KL, *mus_fo_KL; 
	double *netLambda_sigma2byM_fo_KL, *netProton_sigma2byM_fo_KL, *netCharge_sigma2byM_fo_KL, *netKaon_sigma2byM_fo_KL;

	double *energy_fo_PL, *T_fo_PL, *mub_fo_PL, *muq_fo_PL, *mus_fo_PL; 
	double *netLambda_sigma2byM_fo_PL, *netProton_sigma2byM_fo_PL, *netCharge_sigma2byM_fo_PL, *netKaon_sigma2byM_fo_PL;

	DataFileReader infile_self;

	STAR_netLambda_sigma2byM netLambda_sigma2byM_hrg;
	STAR_netProton_sigma2byM netProton_sigma2byM_hrg;
	STAR_netCharge_sigma2byM netCharge_sigma2byM_hrg;
	STAR_netKaon_sigma2byM netKaon_sigma2byM_hrg;

	infile_self.read("data/freezeout_netProtonCharge_fluc.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_fo_PQ = infile_self.nRows();
		energy_fo_PQ = infile_self.getColumn(0);
		T_fo_PQ = infile_self.getColumn(1);
		mub_fo_PQ = infile_self.getColumn(3);
		muq_fo_PQ = infile_self.getColumn(5);
		mus_fo_PQ = infile_self.getColumn(7);
		infile_self.close();

		netLambda_sigma2byM_fo_PQ = new double[intvl_fo_PQ];
		netProton_sigma2byM_fo_PQ = new double[intvl_fo_PQ];
		netCharge_sigma2byM_fo_PQ = new double[intvl_fo_PQ];
		netKaon_sigma2byM_fo_PQ = new double[intvl_fo_PQ];
		for (int i = 0; i < intvl_fo_PQ; ++i)
		{
			T = T_fo_PQ[i];
			mub = mub_fo_PQ[i];
			muq = muq_fo_PQ[i];
			mus = mus_fo_PQ[i];

			netLambda_sigma2byM_fo_PQ[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netProton_sigma2byM_fo_PQ[i] = netProton_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netCharge_sigma2byM_fo_PQ[i] = netCharge_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netKaon_sigma2byM_fo_PQ[i] = netKaon_sigma2byM_hrg.getValue(T,mub,muq,mus);
		}
	}

	infile_self.read("data/freezeout_netKaonCharge_fluc.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_fo_KQ = infile_self.nRows();
		energy_fo_KQ = infile_self.getColumn(0);
		T_fo_KQ = infile_self.getColumn(1);
		mub_fo_KQ = infile_self.getColumn(3);
		muq_fo_KQ = infile_self.getColumn(5);
		mus_fo_KQ = infile_self.getColumn(7);
		infile_self.close();

		netLambda_sigma2byM_fo_KQ = new double[intvl_fo_KQ];
		netProton_sigma2byM_fo_KQ = new double[intvl_fo_KQ];
		netCharge_sigma2byM_fo_KQ = new double[intvl_fo_KQ];
		netKaon_sigma2byM_fo_KQ = new double[intvl_fo_KQ];
		for (int i = 0; i < intvl_fo_KQ; ++i)
		{
			T = T_fo_KQ[i];
			mub = mub_fo_KQ[i];
			muq = muq_fo_KQ[i];
			mus = mus_fo_KQ[i];

			netLambda_sigma2byM_fo_KQ[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netProton_sigma2byM_fo_KQ[i] = netProton_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netCharge_sigma2byM_fo_KQ[i] = netCharge_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netKaon_sigma2byM_fo_KQ[i] = netKaon_sigma2byM_hrg.getValue(T,mub,muq,mus);
		}
	}

	infile_self.read("data/freezeout_netKaonLambda_fluc.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_fo_KL = infile_self.nRows();
		energy_fo_KL = infile_self.getColumn(0);
		T_fo_KL = infile_self.getColumn(1);
		mub_fo_KL = infile_self.getColumn(3);
		muq_fo_KL = infile_self.getColumn(5);
		mus_fo_KL = infile_self.getColumn(7);
		infile_self.close();

		netLambda_sigma2byM_fo_KL = new double[intvl_fo_KL];
		netProton_sigma2byM_fo_KL = new double[intvl_fo_KL];
		netCharge_sigma2byM_fo_KL = new double[intvl_fo_KL];
		netKaon_sigma2byM_fo_KL = new double[intvl_fo_KL];
		for (int i = 0; i < intvl_fo_KL; ++i)
		{
			T = T_fo_KL[i];
			mub = mub_fo_KL[i];
			muq = muq_fo_KL[i];
			mus = mus_fo_KL[i];

			netLambda_sigma2byM_fo_KL[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netProton_sigma2byM_fo_KL[i] = netProton_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netCharge_sigma2byM_fo_KL[i] = netCharge_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netKaon_sigma2byM_fo_KL[i] = netKaon_sigma2byM_hrg.getValue(T,mub,muq,mus);
		}
	}

	infile_self.read("data/freezeout_netProtonLambda_fluc.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_fo_PL = infile_self.nRows();
		energy_fo_PL = infile_self.getColumn(0);
		T_fo_PL = infile_self.getColumn(1);
		mub_fo_PL = infile_self.getColumn(3);
		muq_fo_PL = infile_self.getColumn(5);
		mus_fo_PL = infile_self.getColumn(7);
		infile_self.close();

		netLambda_sigma2byM_fo_PL = new double[intvl_fo_PL];
		netProton_sigma2byM_fo_PL = new double[intvl_fo_PL];
		netCharge_sigma2byM_fo_PL = new double[intvl_fo_PL];
		netKaon_sigma2byM_fo_PL = new double[intvl_fo_PL];
		for (int i = 0; i < intvl_fo_PL; ++i)
		{
			T = T_fo_PL[i];
			mub = mub_fo_PL[i];
			muq = muq_fo_PL[i];
			mus = mus_fo_PL[i];

			netLambda_sigma2byM_fo_PL[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netProton_sigma2byM_fo_PL[i] = netProton_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netCharge_sigma2byM_fo_PL[i] = netCharge_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netKaon_sigma2byM_fo_PL[i] = netKaon_sigma2byM_hrg.getValue(T,mub,muq,mus);
		}
	}

	double topMargin = 0.01;
	double bottomMargin = 0.1;
	double rightMargin = 0.01;
	double leftMargin = 0.1;

	double textSize = 0.05;
	double markerSize = 1.5;
	double lineWidth = 2;
	double transparency = 1;

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetEndErrorSize(3);

	using namespace MyStyle1;

	TCanvas *cnvs = new TCanvas("cnvs_FitResults_Calib","FitResults_Calib",400,350);

	TPad *pad = new TPad("pad_FitResults_Calib", "pad",0,0,1,1);
	pad->SetLogx();
	pad->SetLogy();
	pad->SetBottomMargin(bottomMargin);
	pad->SetTopMargin(topMargin);
	pad->SetLeftMargin(leftMargin);
	pad->SetRightMargin(rightMargin);
	pad->Draw();

	pad->cd();

	auto *gr111 = new TGraphErrors(netLambda_intvl_exp,netLambda_energy_exp,netLambda_sigma2byM_exp_0to5,(double*)0,netLambda_sigma2byM_exp_0to5_staterr);
	gr111->SetMarkerSize(markerSize-0.5);
	gr111->SetMarkerColorAlpha(Color1,transparency);
	gr111->SetLineColorAlpha(Color1,transparency);
	gr111->SetFillColorAlpha(Color1,transparency);
	gr111->SetLineWidth(lineWidth);
	gr111->SetMarkerStyle(Marker_L);
	auto *gr112 = new TGraphErrors(netLambda_intvl_exp,netLambda_energy_exp,netLambda_sigma2byM_exp_0to5,(double*)0,netLambda_sigma2byM_exp_0to5_syserr);
	gr112->SetLineWidth(lineWidth);
	gr112->SetLineColorAlpha(Color1,transparency);

	auto *gr113 = new TGraph(intvl_fo_PL,energy_fo_PL,netLambda_sigma2byM_fo_PL);
	gr113->SetMarkerSize(markerSize);
	gr113->SetMarkerColorAlpha(Color1,transparency);
	gr113->SetMarkerStyle(Marker_PL);

	auto *gr114 = new TGraph(intvl_fo_KL,energy_fo_KL,netLambda_sigma2byM_fo_KL);
	gr114->SetMarkerSize(markerSize);
	gr114->SetMarkerColorAlpha(Color1,transparency);
	gr114->SetMarkerStyle(Marker_KL);

	auto *gr121 = new TGraphErrors(netProton_intvl_exp,netProton_energy_exp,netProton_sigma2byM_exp_0to5,(double*)0,netProton_sigma2byM_exp_0to5_staterr);
	gr121->SetMarkerSize(markerSize-0.5);
	gr121->SetMarkerColorAlpha(Color2,transparency);
	gr121->SetLineColorAlpha(Color2,transparency);
	gr121->SetFillColorAlpha(Color2,transparency);
	gr121->SetLineWidth(lineWidth);
	gr121->SetMarkerStyle(Marker_P);
	auto *gr122 = new TGraphErrors(netProton_intvl_exp,netProton_energy_exp,netProton_sigma2byM_exp_0to5,(double*)0,netProton_sigma2byM_exp_0to5_syserr);
	gr122->SetLineWidth(lineWidth);
	gr122->SetLineColorAlpha(Color2,transparency);

	auto *gr123 = new TGraph(intvl_fo_PL,energy_fo_PL,netProton_sigma2byM_fo_PL);
	gr123->SetMarkerSize(markerSize);
	gr123->SetMarkerColorAlpha(Color2,transparency);
	gr123->SetMarkerStyle(Marker_PL);

	auto *gr124 = new TGraph(intvl_fo_PQ,energy_fo_PQ,netProton_sigma2byM_fo_PQ);
	gr124->SetMarkerSize(markerSize);
	gr124->SetMarkerColorAlpha(Color2,transparency);
	gr124->SetMarkerStyle(Marker_PQ);

	auto *gr131 = new TGraphErrors(netCharge_intvl_exp,netCharge_energy_exp,netCharge_sigma2byM_exp_0to5,(double*)0,netCharge_sigma2byM_exp_0to5_staterr);
	gr131->SetMarkerSize(markerSize-0.5);
	gr131->SetMarkerColorAlpha(Color3,transparency);
	gr131->SetLineColorAlpha(Color3,transparency);
	gr131->SetFillColorAlpha(Color3,transparency);
	gr131->SetLineWidth(lineWidth);
	gr131->SetMarkerStyle(Marker_Q);
	auto *gr132 = new TGraphErrors(netCharge_intvl_exp,netCharge_energy_exp,netCharge_sigma2byM_exp_0to5,(double*)0,netCharge_sigma2byM_exp_0to5_syserr);
	gr132->SetLineWidth(lineWidth);
	gr132->SetLineColorAlpha(Color3,transparency);

	auto *gr133 = new TGraph(intvl_fo_PQ,energy_fo_PQ,netCharge_sigma2byM_fo_PQ);
	gr133->SetMarkerSize(markerSize);
	gr133->SetMarkerColorAlpha(Color3,transparency);
	gr133->SetMarkerStyle(Marker_PQ);

	auto *gr134 = new TGraph(intvl_fo_KQ,energy_fo_KQ,netCharge_sigma2byM_fo_KQ);
	gr134->SetMarkerSize(markerSize);
	gr134->SetMarkerColorAlpha(Color3,transparency);
	gr134->SetMarkerStyle(Marker_KQ);

	auto *gr141 = new TGraphErrors(netKaon_intvl_exp,netKaon_energy_exp,netKaon_sigma2byM_exp_0to5,(double*)0,netKaon_sigma2byM_exp_0to5_staterr);
	gr141->SetMarkerSize(markerSize-0.5);
	gr141->SetMarkerColorAlpha(Color4,transparency);
	gr141->SetLineColorAlpha(Color4,transparency);
	gr141->SetFillColorAlpha(Color4,transparency);
	gr141->SetLineWidth(lineWidth);
	gr141->SetMarkerStyle(Marker_K);
	auto *gr142 = new TGraphErrors(netKaon_intvl_exp,netKaon_energy_exp,netKaon_sigma2byM_exp_0to5,(double*)0,netKaon_sigma2byM_exp_0to5_syserr);
	gr142->SetLineWidth(lineWidth);
	gr142->SetLineColorAlpha(Color4,transparency);

	auto *gr143 = new TGraph(intvl_fo_KL,energy_fo_KL,netKaon_sigma2byM_fo_KL);
	gr143->SetMarkerSize(markerSize);
	gr143->SetMarkerColorAlpha(Color4,transparency);
	gr143->SetMarkerStyle(Marker_KL);

	auto *gr144 = new TGraph(intvl_fo_KQ,energy_fo_KQ,netKaon_sigma2byM_fo_KQ);
	gr144->SetMarkerSize(markerSize);
	gr144->SetMarkerColorAlpha(Color4,transparency);
	gr144->SetMarkerStyle(Marker_KQ);

	TMultiGraph *mg1  = new TMultiGraph();

	mg1->Add(gr111,"p");
	// mg1->Add(gr112,"||");
	mg1->Add(gr113,"p");
	mg1->Add(gr114,"p");
	
	mg1->Add(gr121,"p");
	// mg1->Add(gr122,"||");
	mg1->Add(gr123,"p");
	mg1->Add(gr124,"p");

	mg1->Add(gr131,"p");
	// mg1->Add(gr132,"||");
	mg1->Add(gr133,"p");
	mg1->Add(gr134,"p");

	mg1->Add(gr141,"p");
	// mg1->Add(gr142,"||");
	mg1->Add(gr143,"p");
	mg1->Add(gr144,"p");

	mg1->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg1->GetYaxis()->SetTitle("#sigma^{2}/M");
	mg1->GetXaxis()->SetLimits(8,250);
	mg1->GetYaxis()->SetRangeUser(0.7,200);
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

	TLegend *lgnd1 = new TLegend(0.68,0.12,0.97,0.33);
	lgnd1->SetNColumns(2);
	lgnd1->AddEntry(gr111,"net-#Lambda","p");
	lgnd1->AddEntry(gr121,"net-p","p");
	lgnd1->AddEntry(gr131,"net-q","p");
	lgnd1->AddEntry(gr141,"net-k","p");
	lgnd1->SetTextSize(textSize);
	lgnd1->SetBorderSize(0);
	lgnd1->SetFillStyle(0);
	lgnd1->Draw();

	auto *grM1 = new TGraph(); grM1->SetMarkerStyle(Marker1); grM1->SetMarkerSize(markerSize-0.5); grM1->SetMarkerColorAlpha(kBlack,transparency);
	auto *grM_PL = new TGraph(); grM_PL->SetMarkerStyle(Marker_PL); grM_PL->SetMarkerSize(markerSize); grM_PL->SetMarkerColorAlpha(kBlack,transparency);
	auto *grM_PQ = new TGraph(); grM_PQ->SetMarkerStyle(Marker_PQ);	grM_PQ->SetMarkerSize(markerSize); grM_PQ->SetMarkerColorAlpha(kBlack,transparency);
	auto *grM_KL = new TGraph(); grM_KL->SetMarkerStyle(Marker_KL); grM_KL->SetMarkerSize(markerSize); grM_KL->SetMarkerColorAlpha(kBlack,transparency);
	auto *grM_KQ = new TGraph(); grM_KQ->SetMarkerStyle(Marker_KQ); grM_KQ->SetMarkerSize(markerSize); grM_KQ->SetMarkerColorAlpha(kBlack,transparency);

	TLegend *lgnd2 = new TLegend(0.13,0.63,0.43,0.83);
	lgnd2->SetNColumns(2);
	lgnd2->SetColumnSeparation(-0.5);
	lgnd2->SetMargin(0.35);
	lgnd2->AddEntry(grM1,"STAR Au-Au 0-5%","p");
	lgnd2->AddEntry((TObject*)0,0,0);
	lgnd2->AddEntry(grM_PL,"p#Lambda fit","p");
	lgnd2->AddEntry(grM_PQ,"pq fit","p");
	lgnd2->AddEntry(grM_KL,"k#Lambda fit","p");
	lgnd2->AddEntry(grM_KQ,"kq fit","p");
	lgnd2->SetTextSize(textSize);
	lgnd2->SetBorderSize(0);
	lgnd2->SetFillStyle(0);
	lgnd2->Draw();

	TLatex *txt1 = new TLatex(0.35,0.9,"Fit results (calib.)");
	txt1->SetNDC();
	txt1->SetTextSize(textSize);
	txt1->Draw();

	cnvs->SaveAs("plots/cumulantRatios_fitResults_Calib.pdf");
}

void plotFitResults_Estimate()
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

	int intvl_fo_PQ, intvl_fo_KQ, intvl_fo_KL, intvl_fo_PL;

	double *energy_fo_PQ, *T_fo_PQ, *mub_fo_PQ, *muq_fo_PQ, *mus_fo_PQ; 
	double *netLambda_sigma2byM_fo_PQ, *netProton_sigma2byM_fo_PQ, *netCharge_sigma2byM_fo_PQ, *netKaon_sigma2byM_fo_PQ;

	double *energy_fo_KQ, *T_fo_KQ, *mub_fo_KQ, *muq_fo_KQ, *mus_fo_KQ; 
	double *netLambda_sigma2byM_fo_KQ, *netProton_sigma2byM_fo_KQ, *netCharge_sigma2byM_fo_KQ, *netKaon_sigma2byM_fo_KQ;

	double *energy_fo_KL, *T_fo_KL, *mub_fo_KL, *muq_fo_KL, *mus_fo_KL; 
	double *netLambda_sigma2byM_fo_KL, *netProton_sigma2byM_fo_KL, *netCharge_sigma2byM_fo_KL, *netKaon_sigma2byM_fo_KL;

	double *energy_fo_PL, *T_fo_PL, *mub_fo_PL, *muq_fo_PL, *mus_fo_PL; 
	double *netLambda_sigma2byM_fo_PL, *netProton_sigma2byM_fo_PL, *netCharge_sigma2byM_fo_PL, *netKaon_sigma2byM_fo_PL;

	DataFileReader infile_self;

	STAR_netLambda_sigma2byM netLambda_sigma2byM_hrg;
	STAR_netProton_sigma2byM netProton_sigma2byM_hrg;
	STAR_netCharge_sigma2byM netCharge_sigma2byM_hrg;
	STAR_netKaon_sigma2byM netKaon_sigma2byM_hrg;

	infile_self.read("data/freezeout_netProtonCharge_fluc.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_fo_PQ = infile_self.nRows();
		energy_fo_PQ = infile_self.getColumn(0);
		T_fo_PQ = infile_self.getColumn(1);
		mub_fo_PQ = infile_self.getColumn(3);
		muq_fo_PQ = infile_self.getColumn(5);
		mus_fo_PQ = infile_self.getColumn(7);
		infile_self.close();

		netLambda_sigma2byM_fo_PQ = new double[intvl_fo_PQ];
		netProton_sigma2byM_fo_PQ = new double[intvl_fo_PQ];
		netCharge_sigma2byM_fo_PQ = new double[intvl_fo_PQ];
		netKaon_sigma2byM_fo_PQ = new double[intvl_fo_PQ];
		for (int i = 0; i < intvl_fo_PQ; ++i)
		{
			T = T_fo_PQ[i];
			mub = mub_fo_PQ[i];
			muq = muq_fo_PQ[i];
			mus = mus_fo_PQ[i];

			netLambda_sigma2byM_fo_PQ[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netProton_sigma2byM_fo_PQ[i] = netProton_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netCharge_sigma2byM_fo_PQ[i] = netCharge_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netKaon_sigma2byM_fo_PQ[i] = netKaon_sigma2byM_hrg.getValue(T,mub,muq,mus);
		}
	}

	infile_self.read("data/freezeout_netKaonCharge_fluc.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_fo_KQ = infile_self.nRows();
		energy_fo_KQ = infile_self.getColumn(0);
		T_fo_KQ = infile_self.getColumn(1);
		mub_fo_KQ = infile_self.getColumn(3);
		muq_fo_KQ = infile_self.getColumn(5);
		mus_fo_KQ = infile_self.getColumn(7);
		infile_self.close();

		netLambda_sigma2byM_fo_KQ = new double[intvl_fo_KQ];
		netProton_sigma2byM_fo_KQ = new double[intvl_fo_KQ];
		netCharge_sigma2byM_fo_KQ = new double[intvl_fo_KQ];
		netKaon_sigma2byM_fo_KQ = new double[intvl_fo_KQ];
		for (int i = 0; i < intvl_fo_KQ; ++i)
		{
			T = T_fo_KQ[i];
			mub = mub_fo_KQ[i];
			muq = muq_fo_KQ[i];
			mus = mus_fo_KQ[i];

			netLambda_sigma2byM_fo_KQ[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netProton_sigma2byM_fo_KQ[i] = netProton_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netCharge_sigma2byM_fo_KQ[i] = netCharge_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netKaon_sigma2byM_fo_KQ[i] = netKaon_sigma2byM_hrg.getValue(T,mub,muq,mus);
		}
	}

	infile_self.read("data/freezeout_netKaonLambda_fluc.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_fo_KL = infile_self.nRows();
		energy_fo_KL = infile_self.getColumn(0);
		T_fo_KL = infile_self.getColumn(1);
		mub_fo_KL = infile_self.getColumn(3);
		muq_fo_KL = infile_self.getColumn(5);
		mus_fo_KL = infile_self.getColumn(7);
		infile_self.close();

		netLambda_sigma2byM_fo_KL = new double[intvl_fo_KL];
		netProton_sigma2byM_fo_KL = new double[intvl_fo_KL];
		netCharge_sigma2byM_fo_KL = new double[intvl_fo_KL];
		netKaon_sigma2byM_fo_KL = new double[intvl_fo_KL];
		for (int i = 0; i < intvl_fo_KL; ++i)
		{
			T = T_fo_KL[i];
			mub = mub_fo_KL[i];
			muq = muq_fo_KL[i];
			mus = mus_fo_KL[i];

			netLambda_sigma2byM_fo_KL[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netProton_sigma2byM_fo_KL[i] = netProton_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netCharge_sigma2byM_fo_KL[i] = netCharge_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netKaon_sigma2byM_fo_KL[i] = netKaon_sigma2byM_hrg.getValue(T,mub,muq,mus);
		}
	}

	infile_self.read("data/freezeout_netProtonLambda_fluc.dat");
	if(infile_self.isDataLoaded())
	{
		intvl_fo_PL = infile_self.nRows();
		energy_fo_PL = infile_self.getColumn(0);
		T_fo_PL = infile_self.getColumn(1);
		mub_fo_PL = infile_self.getColumn(3);
		muq_fo_PL = infile_self.getColumn(5);
		mus_fo_PL = infile_self.getColumn(7);
		infile_self.close();

		netLambda_sigma2byM_fo_PL = new double[intvl_fo_PL];
		netProton_sigma2byM_fo_PL = new double[intvl_fo_PL];
		netCharge_sigma2byM_fo_PL = new double[intvl_fo_PL];
		netKaon_sigma2byM_fo_PL = new double[intvl_fo_PL];
		for (int i = 0; i < intvl_fo_PL; ++i)
		{
			T = T_fo_PL[i];
			mub = mub_fo_PL[i];
			muq = muq_fo_PL[i];
			mus = mus_fo_PL[i];

			netLambda_sigma2byM_fo_PL[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netProton_sigma2byM_fo_PL[i] = netProton_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netCharge_sigma2byM_fo_PL[i] = netCharge_sigma2byM_hrg.getValue(T,mub,muq,mus);
			netKaon_sigma2byM_fo_PL[i] = netKaon_sigma2byM_hrg.getValue(T,mub,muq,mus);
		}
	}

	double E_l = 5;
	double E_u = 250;
	int intvl = 100;
	double dE = pow(E_u/E_l,1.0/(intvl-1));
	double x;

	double energy_cleymans[intvl];

	double netLambda_sigma2byM_cleymans[intvl];
	double netProton_sigma2byM_cleymans[intvl];
	double netCharge_sigma2byM_cleymans[intvl];
	double netKaon_sigma2byM_cleymans[intvl];

	for (int i = 0; i < intvl; ++i)
	{
		x = E_l*pow(dE,i);
		energy_cleymans[i] = x;

		mub = 1.308/(1+0.273*x);
		muq = -0.0202/(1+0.126*x);
		mus = 0.214/(1+0.184*x);

		T = 0.166-0.139*mub*mub-0.053*pow(mub,4);

		netLambda_sigma2byM_cleymans[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
		netProton_sigma2byM_cleymans[i] = netProton_sigma2byM_hrg.getValue(T,mub,muq,mus);
		netCharge_sigma2byM_cleymans[i] = netCharge_sigma2byM_hrg.getValue(T,mub,muq,mus);
		netKaon_sigma2byM_cleymans[i] = netKaon_sigma2byM_hrg.getValue(T,mub,muq,mus);
	}

	double topMargin = 0.01;
	double bottomMargin = 0.1;
	double rightMargin = 0.01;
	double leftMargin = 0.1;

	double textSize = 0.05;
	double markerSize = 1.5;
	double lineWidth = 2;
	double transparency = 1;

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetEndErrorSize(3);

	using namespace MyStyle1;

	TCanvas *cnvs = new TCanvas("cnvs_FitResultsEstimate","FitResults_Estimate",400,350);

	TPad *pad = new TPad("pad_FitResultsEstimate", "pad",0,0,1,1);
	pad->SetLogx();
	pad->SetLogy();
	pad->SetBottomMargin(bottomMargin);
	pad->SetTopMargin(topMargin);
	pad->SetLeftMargin(leftMargin);
	pad->SetRightMargin(rightMargin);
	pad->Draw();

	pad->cd();

	auto *gr110 = new TGraph(intvl,energy_cleymans,netLambda_sigma2byM_cleymans);
	gr110->SetLineStyle(Line1);
	gr110->SetLineColorAlpha(Color1,transparency);
	gr110->SetLineWidth(lineWidth);

	auto *gr111 = new TGraphErrors(netLambda_intvl_exp,netLambda_energy_exp,netLambda_sigma2byM_exp_0to5,(double*)0,netLambda_sigma2byM_exp_0to5_staterr);
	gr111->SetMarkerSize(markerSize-0.5);
	gr111->SetMarkerColorAlpha(Color1,transparency);
	gr111->SetLineColorAlpha(Color1,transparency);
	gr111->SetFillColorAlpha(Color1,transparency);
	gr111->SetLineWidth(lineWidth);
	gr111->SetMarkerStyle(Marker_L);
	auto *gr112 = new TGraphErrors(netLambda_intvl_exp,netLambda_energy_exp,netLambda_sigma2byM_exp_0to5,(double*)0,netLambda_sigma2byM_exp_0to5_syserr);
	gr112->SetLineWidth(lineWidth);
	gr112->SetLineColorAlpha(Color1,transparency);

	auto *gr113 = new TGraph(intvl_fo_PQ,energy_fo_PQ,netLambda_sigma2byM_fo_PQ);
	gr113->SetMarkerSize(markerSize);
	gr113->SetMarkerColorAlpha(Color1,transparency);
	gr113->SetMarkerStyle(Marker_PQ);

	auto *gr114 = new TGraph(intvl_fo_KQ,energy_fo_KQ,netLambda_sigma2byM_fo_KQ);
	gr114->SetMarkerSize(markerSize);
	gr114->SetMarkerColorAlpha(Color1,transparency);
	gr114->SetMarkerStyle(Marker_KQ);


	auto *gr120 = new TGraph(intvl,energy_cleymans,netProton_sigma2byM_cleymans);
	gr120->SetLineStyle(Line1);
	gr120->SetLineColorAlpha(Color2,transparency);
	gr120->SetLineWidth(lineWidth);

	auto *gr121 = new TGraphErrors(netProton_intvl_exp,netProton_energy_exp,netProton_sigma2byM_exp_0to5,(double*)0,netProton_sigma2byM_exp_0to5_staterr);
	gr121->SetMarkerSize(markerSize-0.5);
	gr121->SetMarkerColorAlpha(Color2,transparency);
	gr121->SetLineColorAlpha(Color2,transparency);
	gr121->SetFillColorAlpha(Color2,transparency);
	gr121->SetLineWidth(lineWidth);
	gr121->SetMarkerStyle(Marker_P);
	auto *gr122 = new TGraphErrors(netProton_intvl_exp,netProton_energy_exp,netProton_sigma2byM_exp_0to5,(double*)0,netProton_sigma2byM_exp_0to5_syserr);
	gr122->SetLineWidth(lineWidth);
	gr122->SetLineColorAlpha(Color2,transparency);

	auto *gr123 = new TGraph(intvl_fo_KL,energy_fo_KL,netProton_sigma2byM_fo_KL);
	gr123->SetMarkerSize(markerSize);
	gr123->SetMarkerColorAlpha(Color2,transparency);
	gr123->SetMarkerStyle(Marker_KL);

	auto *gr124 = new TGraph(intvl_fo_KQ,energy_fo_KQ,netProton_sigma2byM_fo_KQ);
	gr124->SetMarkerSize(markerSize);
	gr124->SetMarkerColorAlpha(Color2,transparency);
	gr124->SetMarkerStyle(Marker_KQ);


	auto *gr130 = new TGraph(intvl,energy_cleymans,netCharge_sigma2byM_cleymans);
	gr130->SetLineStyle(Line1);
	gr130->SetLineColorAlpha(Color3,transparency);
	gr130->SetLineWidth(lineWidth);

	auto *gr131 = new TGraphErrors(netCharge_intvl_exp,netCharge_energy_exp,netCharge_sigma2byM_exp_0to5,(double*)0,netCharge_sigma2byM_exp_0to5_staterr);
	gr131->SetMarkerSize(markerSize-0.5);
	gr131->SetMarkerColorAlpha(Color3,transparency);
	gr131->SetLineColorAlpha(Color3,transparency);
	gr131->SetFillColorAlpha(Color3,transparency);
	gr131->SetLineWidth(lineWidth);
	gr131->SetMarkerStyle(Marker_Q);
	auto *gr132 = new TGraphErrors(netCharge_intvl_exp,netCharge_energy_exp,netCharge_sigma2byM_exp_0to5,(double*)0,netCharge_sigma2byM_exp_0to5_syserr);
	gr132->SetLineWidth(lineWidth);
	gr132->SetLineColorAlpha(Color3,transparency);

	auto *gr133 = new TGraph(intvl_fo_PL,energy_fo_PL,netCharge_sigma2byM_fo_PL);
	gr133->SetMarkerSize(markerSize);
	gr133->SetMarkerColorAlpha(Color3,transparency);
	gr133->SetMarkerStyle(Marker_PL);

	auto *gr134 = new TGraph(intvl_fo_KL,energy_fo_KL,netCharge_sigma2byM_fo_KL);
	gr134->SetMarkerSize(markerSize);
	gr134->SetMarkerColorAlpha(Color3,transparency);
	gr134->SetMarkerStyle(Marker_KL);


	auto *gr140 = new TGraph(intvl,energy_cleymans,netKaon_sigma2byM_cleymans);
	gr140->SetLineStyle(Line1);
	gr140->SetLineColorAlpha(Color4,transparency);
	gr140->SetLineWidth(lineWidth);

	auto *gr141 = new TGraphErrors(netKaon_intvl_exp,netKaon_energy_exp,netKaon_sigma2byM_exp_0to5,(double*)0,netKaon_sigma2byM_exp_0to5_staterr);
	gr141->SetMarkerSize(markerSize-0.5);
	gr141->SetMarkerColorAlpha(Color4,transparency);
	gr141->SetLineColorAlpha(Color4,transparency);
	gr141->SetFillColorAlpha(Color4,transparency);
	gr141->SetLineWidth(lineWidth);
	gr141->SetMarkerStyle(Marker_K);
	auto *gr142 = new TGraphErrors(netKaon_intvl_exp,netKaon_energy_exp,netKaon_sigma2byM_exp_0to5,(double*)0,netKaon_sigma2byM_exp_0to5_syserr);
	gr142->SetLineWidth(lineWidth);
	gr142->SetLineColorAlpha(Color4,transparency);

	auto *gr143 = new TGraph(intvl_fo_PL,energy_fo_PL,netKaon_sigma2byM_fo_PL);
	gr143->SetMarkerSize(markerSize);
	gr143->SetMarkerColorAlpha(Color4,transparency);
	gr143->SetMarkerStyle(Marker_PL);

	auto *gr144 = new TGraph(intvl_fo_PQ,energy_fo_PQ,netKaon_sigma2byM_fo_PQ);
	gr144->SetMarkerSize(markerSize);
	gr144->SetMarkerColorAlpha(Color4,transparency);
	gr144->SetMarkerStyle(Marker_PQ);

	TMultiGraph *mg1  = new TMultiGraph();


	mg1->Add(gr110,"l");
	mg1->Add(gr111,"p");
	// mg1->Add(gr112,"||");
	mg1->Add(gr113,"p");
	mg1->Add(gr114,"p");
	
	mg1->Add(gr120,"l");
	mg1->Add(gr121,"p");
	// mg1->Add(gr122,"||");
	mg1->Add(gr123,"p");
	mg1->Add(gr124,"p");

	mg1->Add(gr130,"l");
	mg1->Add(gr131,"p");
	// mg1->Add(gr132,"||");
	mg1->Add(gr133,"p");
	mg1->Add(gr134,"p");

	mg1->Add(gr140,"l");
	mg1->Add(gr141,"p");
	// mg1->Add(gr142,"||");
	mg1->Add(gr143,"p");
	mg1->Add(gr144,"p");

	mg1->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg1->GetYaxis()->SetTitle("#sigma^{2}/M");
	mg1->GetXaxis()->SetLimits(8,250);
	mg1->GetYaxis()->SetRangeUser(0.7,200);
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

	TLegend *lgnd1 = new TLegend(0.68,0.12,0.97,0.32);
	lgnd1->SetNColumns(2);
	lgnd1->AddEntry(gr111,"net-#Lambda","p");
	lgnd1->AddEntry(gr121,"net-p","p");
	lgnd1->AddEntry(gr131,"net-q","p");
	lgnd1->AddEntry(gr141,"net-k","p");
	lgnd1->SetTextSize(textSize);
	lgnd1->SetBorderSize(0);
	lgnd1->SetFillStyle(0);
	lgnd1->Draw();

	auto *grL1 = new TGraph(); grL1->SetLineStyle(Line1); grL1->SetLineWidth(lineWidth); grL1->SetLineColorAlpha(kBlack,transparency);
	auto *grM1 = new TGraph(); grM1->SetMarkerStyle(Marker1); grM1->SetMarkerSize(markerSize-0.5); grM1->SetMarkerColorAlpha(kBlack,transparency);
	auto *grM_PL = new TGraph(); grM_PL->SetMarkerStyle(Marker_PL); grM_PL->SetMarkerSize(markerSize); grM_PL->SetMarkerColorAlpha(kBlack,transparency);
	auto *grM_PQ = new TGraph(); grM_PQ->SetMarkerStyle(Marker_PQ);	grM_PQ->SetMarkerSize(markerSize); grM_PQ->SetMarkerColorAlpha(kBlack,transparency);
	auto *grM_KL = new TGraph(); grM_KL->SetMarkerStyle(Marker_KL); grM_KL->SetMarkerSize(markerSize); grM_KL->SetMarkerColorAlpha(kBlack,transparency);
	auto *grM_KQ = new TGraph(); grM_KQ->SetMarkerStyle(Marker_KQ); grM_KQ->SetMarkerSize(markerSize); grM_KQ->SetMarkerColorAlpha(kBlack,transparency);

	TLegend *lgnd2 = new TLegend(0.13,0.67,0.43,0.89);
	lgnd2->SetNColumns(2);
	lgnd2->SetColumnSeparation(-0.5);
	lgnd2->SetMargin(0.35);
	lgnd2->AddEntry(grM1,"STAR Au-Au 0-5%","p");
	lgnd2->AddEntry((TObject*)0,0,0);
	lgnd2->AddEntry(grM_PL,"p#Lambda fit","p");
	lgnd2->AddEntry(grM_PQ,"pq fit","p");
	lgnd2->AddEntry(grM_KL,"k#Lambda fit","p");
	lgnd2->AddEntry(grM_KQ,"kq fit","p");
	lgnd2->AddEntry(grL1,"yield fit (Cleymans et. al)","l");
	lgnd2->SetTextSize(textSize);
	lgnd2->SetBorderSize(0);
	lgnd2->SetFillStyle(0);
	lgnd2->Draw();

	TLatex *txt1 = new TLatex(0.35,0.93,"Fit results (estimate)");
	txt1->SetNDC();
	txt1->SetTextSize(textSize);
	txt1->Draw();

	cnvs->SaveAs("plots/cumulantRatios_fitResults_Estimate.pdf");
}

void plotIndiFitResults_Calib()
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

	TCanvas *cnvs = new TCanvas("cnvs_IndiFitResults_Calib","IndiFitResults_Calib",400,350);

	TPad *pad = new TPad("pad_IndiFitResults_Calib", "pad",0,0,1,1);
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

	auto *gr113 = new TGraph(intvl_fo_P,energy_fo_P,netProton_sigma2byM_fo_P);
	gr113->SetMarkerSize(markerSize);
	gr113->SetMarkerColorAlpha(Color_P,opaque);
	gr113->SetMarkerStyle(Marker_P_hrg);

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

	auto *gr123 = new TGraph(intvl_fo_K,energy_fo_K,netKaon_sigma2byM_fo_K);
	gr123->SetMarkerSize(markerSize);
	gr123->SetMarkerColorAlpha(Color_K,opaque);
	gr123->SetMarkerStyle(Marker_K_hrg);

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

	auto *gr133 = new TGraph(intvl_fo_L,energy_fo_L,netLambda_sigma2byM_fo_L);
	gr133->SetMarkerSize(markerSize);
	gr133->SetMarkerColorAlpha(Color_L,opaque);
	gr133->SetMarkerStyle(Marker_L_hrg);

	TMultiGraph *mg1  = new TMultiGraph();

	mg1->Add(gr111,"pl3");
	// mg1->Add(gr112,"||");
	mg1->Add(gr113,"p");
	
	mg1->Add(gr121,"pl3");
	// mg1->Add(gr122,"||");
	mg1->Add(gr123,"p");

	mg1->Add(gr131,"pl3");
	// mg1->Add(gr132,"||");
	mg1->Add(gr133,"p");

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

	TText *txt1 = new TText(0.35,0.9,"Fit results (calib.)");
	txt1->SetNDC();
	txt1->SetTextSize(textSize);
	txt1->Draw();

	// TLatex *txt2 = new TLatex(0.16,0.73,"#splitline{solid #rightarrow STAR Au-Au 0-5%}{hollow #rightarrow HRG fit}");
	// txt2->SetNDC();
	// txt2->SetTextSize(textSize);
	// txt2->Draw();

	cnvs->SaveAs("plots/cumulantRatios_IndiFitResults_Calib.pdf");
}

void plotIndiFitResults_Estimate()
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

	TCanvas *cnvs = new TCanvas("cnvs_IndiFitResults_Estimate","IndiFitResults_Estimate",400,350);

	TPad *pad = new TPad("pad_IndiFitResults_Estimate", "pad",0,0,1,1);
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

	TMultiGraph *mg1  = new TMultiGraph();

	mg1->Add(gr111,"pl3");
	// mg1->Add(gr112,"||");
	mg1->Add(gr113,"p");
	mg1->Add(gr114,"p");
	
	mg1->Add(gr121,"pl3");
	// mg1->Add(gr122,"||");
	mg1->Add(gr123,"p");
	mg1->Add(gr124,"p");

	mg1->Add(gr131,"pl3");
	// mg1->Add(gr132,"||");
	mg1->Add(gr133,"p");
	mg1->Add(gr134,"p");

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

	TLatex *txt1 = new TLatex(0.35,0.9,"Fit results (estimate)");
	txt1->SetNDC();
	txt1->SetTextSize(textSize);
	txt1->Draw();

	cnvs->SaveAs("plots/cumulantRatios_IndiFitResults_Estimate.pdf");
}

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

	double E_l = 5;
	double E_u = 250;
	int intvl = 100;
	double dE = pow(E_u/E_l,1.0/(intvl-1));
	double x;

	double energy_cleymans[intvl];

	double netLambda_sigma2byM_cleymans[intvl];
	double netProton_sigma2byM_cleymans[intvl];
	double netKaon_sigma2byM_cleymans[intvl];

	for (int i = 0; i < intvl; ++i)
	{
		x = E_l*pow(dE,i);
		energy_cleymans[i] = x;

		mub = 1.308/(1+0.273*x);
		muq = -0.0202/(1+0.126*x);
		mus = 0.214/(1+0.184*x);

		T = 0.166-0.139*mub*mub-0.053*pow(mub,4);

		netLambda_sigma2byM_cleymans[i] = netLambda_sigma2byM_hrg.getValue(T,mub,muq,mus);
		netProton_sigma2byM_cleymans[i] = netProton_sigma2byM_hrg.getValue(T,mub,muq,mus);
		netKaon_sigma2byM_cleymans[i] = netKaon_sigma2byM_hrg.getValue(T,mub,muq,mus);
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

	int Line_Cleymans = kDashed;
	double lineWidth_Cleymans = 2;

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

	auto *gr110 = new TGraph(intvl,energy_cleymans,netProton_sigma2byM_cleymans);
	gr110->SetLineStyle(Line_Cleymans);
	gr110->SetLineColorAlpha(Color_P,opaque);
	gr110->SetLineWidth(lineWidth_Cleymans);

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


	auto *gr120 = new TGraph(intvl,energy_cleymans,netKaon_sigma2byM_cleymans);
	gr120->SetLineStyle(Line_Cleymans);
	gr120->SetLineColorAlpha(Color_K,opaque);
	gr120->SetLineWidth(lineWidth_Cleymans);

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
	gr125->SetMarkerSize(markerSize);
	gr125->SetMarkerColorAlpha(Color_K,opaque);
	gr125->SetMarkerStyle(Marker_K_hrg);


	auto *gr130 = new TGraph(intvl,energy_cleymans,netLambda_sigma2byM_cleymans);
	gr130->SetLineStyle(Line_Cleymans);
	gr130->SetLineColorAlpha(Color_L,opaque);
	gr130->SetLineWidth(lineWidth_Cleymans);

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

	mg1->Add(gr110,"l");
	mg1->Add(gr111,"pl3");
	// mg1->Add(gr112,"||");
	mg1->Add(gr113,"p");
	mg1->Add(gr114,"p");
	mg1->Add(gr115,"p");
	
	mg1->Add(gr120,"l");
	mg1->Add(gr121,"pl3");
	// mg1->Add(gr122,"||");
	mg1->Add(gr123,"p");
	mg1->Add(gr124,"p");
	mg1->Add(gr125,"p");

	mg1->Add(gr130,"l");
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

	auto *grY = new TGraph();
	grY->SetMarkerSize(markerSize);
	grY->SetMarkerColorAlpha(kBlack,opaque);
	grY->SetLineStyle(Line_Cleymans);
	grY->SetLineWidth(lineWidth_Cleymans);

	TLegend *lgnd1 = new TLegend(0.15,0.73,0.56,0.95);
	lgnd1->SetNColumns(2);
	lgnd1->SetHeader("HRG Fits(data used)","C");
	lgnd1->AddEntry(grP,"net-p","p");
	lgnd1->AddEntry(grK,"net-k","p");
	lgnd1->AddEntry(grL,"net-#Lambda","p");
	lgnd1->AddEntry(grY,"yield","l");
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

int main()
{
	TApplication app {"app",0,0};
	// plotNetLambdaResults();
	// plotFitResults_Calib();
	// plotFitResults_Estimate();
	// plotIndiFitResults_Calib();
	// plotIndiFitResults_Estimate();
	plotIndiFitResults();
	
	cout << "\nPress Ctrl-C to exit..." << endl;

	app.Run();
	getchar();
	return 0;
}