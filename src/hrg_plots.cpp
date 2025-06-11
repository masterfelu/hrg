#include <iostream>
#include <cmath>
#include <fstream>

#include "Particles.h"
#include "ThermalFunctions.h"
#include "NumericalMinimization.h"

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

class Sigma2byM_Proton:public ThermalFunction
{

	Susceptibility1B *chi1B;
	Susceptibility2B *chi2B;

	public:

	Sigma2byM_Proton()
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
	~Sigma2byM_Proton()
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
		// std::cout << "chi2P/chi1P(" << T << ", " << mub << ", " << muq << ", " << mus << ") : " << num/denom << '\n';
		return num/denom;
	}
};

class Sigma2byM_Charge:public ThermalFunction
{

	Susceptibility1Q *chi1Q;
	Susceptibility2Q *chi2Q;

	public:

	Sigma2byM_Charge()
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
	~Sigma2byM_Charge()
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

class StrangenessNumberDensity:public ThermalFunction
{
	Susceptibility1S *chi1S;

	public:

	StrangenessNumberDensity()
	{
		ParticleSystem system;
		system.loadDefaultData();
		chi1S = new Susceptibility1S {system};
		chi1S->setSphericalUniformCoordinates();
		chi1S->setIntegralLimits(0,INFINITY);
		chi1S->setDecay(true);
	}
	~StrangenessNumberDensity()
	{
		delete chi1S;
		chi1S = nullptr;
	}
	double getValue(double T, double mub, double muq, double mus)
	{
		return chi1S->getValue(T,mub,muq,mus)*T*T*T;
	}
};

class ChargeToBaryonRatio:public ThermalFunction
{
	public:

	Susceptibility1B *chi1B;
	Susceptibility1Q *chi1Q;
	ChargeToBaryonRatio()
	{
		ParticleSystem system;
		system.loadDefaultData();

		chi1B = new Susceptibility1B {system};
		chi1B->setSphericalUniformCoordinates();
		chi1B->setIntegralLimits(0,INFINITY);
		chi1B->setDecay(true);

		chi1Q = new Susceptibility1Q {system};
		chi1Q->setSphericalUniformCoordinates();
		chi1Q->setIntegralLimits(0,INFINITY);
		chi1Q->setDecay(true);
	}
	~ChargeToBaryonRatio()
	{
		delete chi1B;
		delete chi1Q;
		chi1B = nullptr;
		chi1Q = nullptr;
	}
	double getValue(double T, double mub, double muq, double mus)
	{
		double num = chi1Q->getValue(T,mub,muq,mus);
		double denom = chi1B->getValue(T,mub,muq,mus);
		return num/denom;
	}
};

void numericalMinimization()
{
	Sigma2byM_Proton chiP;
	Sigma2byM_Charge chiQ;
	StrangenessNumberDensity nS;
	ChargeToBaryonRatio nQ_by_nB;
	cout << "-----------------------------------------------------------------" << '\n';
	cout << "        Determination of chemical freezout parameters using      " << '\n';
	cout << "             net-proton and net-charge fluctuations              " << '\n';
	cout << "-----------------------------------------------------------------" << '\n';
	ThermalMinimizer minimum {4};
	minimum.SetThermalFunctionEach(0,&chiP);
	minimum.SetThermalFunctionEach(1,&chiQ);
	minimum.SetThermalFunctionEach(2,&nS);
	minimum.SetThermalFunctionEach(3,&nQ_by_nB);

	minimum.SetFunctionValueEach(0,6.178635);
	minimum.SetFunctionValueEach(1,80.2318);
	minimum.SetFunctionValueEach(2,0);
	minimum.SetFunctionValueEach(3,0.4);

	minimum.setInitialParameters(0.145,0.02,-0.0007,0.005);

	minimum.minimize(100);
}

void numericalMinimizationCumulants()
{
	int n_expData;

	double *snn_exp;

	double *netP_sigma2byM_fit;

	double *netP_sigma2byM_0to5, *netP_sigma2byM_0to5_staterr, *netP_sigma2byM_0to5_syserr;
	double *netP_Ssigma_0to5, *netP_Ssigma_0to5_syserr, *netP_Ssigma_0to5_staterr;
	double *netP_ksigma2_0to5, *netP_ksigma2_0to5_syserr, *netP_ksigma2_0to5_staterr;

	double *netP_sigma2byM_70to80, *netP_sigma2byM_70to80_staterr, *netP_sigma2byM_70to80_syserr;
	double *netP_Ssigma_70to80, *netP_Ssigma_70to80_syserr, *netP_Ssigma_70to80_staterr;
	double *netP_ksigma2_70to80, *netP_ksigma2_70to80_syserr, *netP_ksigma2_70to80_staterr;

	double *netQ_sigma2byM_fit;

	double *netQ_sigma2byM_0to5, *netQ_sigma2byM_0to5_staterr, *netQ_sigma2byM_0to5_syserr;
	double *netQ_Ssigma_0to5, *netQ_Ssigma_0to5_syserr, *netQ_Ssigma_0to5_staterr;
	double *netQ_ksigma2_0to5, *netQ_ksigma2_0to5_syserr, *netQ_ksigma2_0to5_staterr;

	double *netQ_sigma2byM_70to80, *netQ_sigma2byM_70to80_staterr, *netQ_sigma2byM_70to80_syserr;
	double *netQ_Ssigma_70to80, *netQ_Ssigma_70to80_syserr, *netQ_Ssigma_70to80_staterr;
	double *netQ_ksigma2_70to80, *netQ_ksigma2_70to80_syserr, *netQ_ksigma2_70to80_staterr;

	ifstream inFileNetP;
	ifstream inFileNetQ;
	inFileNetP.open("data/cumulantRatios_netP.dat");
	inFileNetQ.open("data/momentRatios_netQ.dat");
	if (inFileNetP && inFileNetQ)
	{
		int l;
		inFileNetP >> l;
		n_expData = l;
		inFileNetQ >> l;
		n_expData = min(l,n_expData);
		l = n_expData;

		snn_exp = new double[l];

		netP_sigma2byM_fit = new double[l];

		netP_sigma2byM_0to5 = new double[l]; netP_sigma2byM_0to5_syserr = new double[l]; netP_sigma2byM_0to5_staterr = new double[l];
		netP_Ssigma_0to5 = new double[l]; netP_Ssigma_0to5_syserr = new double[l]; netP_Ssigma_0to5_staterr = new double[l];
		netP_ksigma2_0to5 = new double[l]; netP_ksigma2_0to5_syserr = new double[l]; netP_ksigma2_0to5_staterr = new double[l];

		netP_sigma2byM_70to80 = new double[l]; netP_sigma2byM_70to80_syserr = new double[l]; netP_sigma2byM_70to80_staterr = new double[l];
		netP_Ssigma_70to80 = new double[l]; netP_Ssigma_70to80_syserr = new double[l]; netP_Ssigma_70to80_staterr = new double[l];
		netP_ksigma2_70to80 = new double[l]; netP_ksigma2_70to80_syserr = new double[l]; netP_ksigma2_70to80_staterr = new double[l];

		netQ_sigma2byM_fit = new double[l];

		netQ_sigma2byM_0to5 = new double[l]; netQ_sigma2byM_0to5_syserr = new double[l]; netQ_sigma2byM_0to5_staterr = new double[l];
		netQ_Ssigma_0to5 = new double[l]; netQ_Ssigma_0to5_syserr = new double[l]; netQ_Ssigma_0to5_staterr = new double[l];
		netQ_ksigma2_0to5 = new double[l]; netQ_ksigma2_0to5_syserr = new double[l]; netQ_ksigma2_0to5_staterr = new double[l];

		netQ_sigma2byM_70to80 = new double[l]; netQ_sigma2byM_70to80_syserr = new double[l]; netQ_sigma2byM_70to80_staterr = new double[l];
		netQ_Ssigma_70to80 = new double[l]; netQ_Ssigma_70to80_syserr = new double[l]; netQ_Ssigma_70to80_staterr = new double[l];
		netQ_ksigma2_70to80 = new double[l]; netQ_ksigma2_70to80_syserr = new double[l]; netQ_ksigma2_70to80_staterr = new double[l];

		for (int i = 0; i < l; ++i)
		{
			inFileNetP >> snn_exp[i] >> netP_sigma2byM_0to5[i] >> netP_sigma2byM_0to5_staterr[i] >> netP_sigma2byM_0to5_syserr[i]
			 >> netP_Ssigma_0to5[i] >> netP_Ssigma_0to5_staterr[i] >> netP_Ssigma_0to5_syserr[i]
			 >> netP_ksigma2_0to5[i] >> netP_ksigma2_0to5_staterr[i] >> netP_ksigma2_0to5_syserr[i]
			 >> netP_sigma2byM_70to80[i] >> netP_sigma2byM_70to80_staterr[i] >> netP_sigma2byM_70to80_syserr[i]
			 >> netP_Ssigma_70to80[i] >> netP_Ssigma_70to80_staterr[i] >> netP_Ssigma_70to80_syserr[i]
			 >> netP_ksigma2_70to80[i] >> netP_ksigma2_70to80_staterr[i] >> netP_ksigma2_70to80_syserr[i];

			inFileNetQ >> snn_exp[i] >> netQ_sigma2byM_0to5[i] >> netQ_sigma2byM_0to5_staterr[i] >> netQ_sigma2byM_0to5_syserr[i]
			 >> netQ_Ssigma_0to5[i] >> netQ_Ssigma_0to5_staterr[i] >> netQ_Ssigma_0to5_syserr[i]
			 >> netQ_ksigma2_0to5[i] >> netQ_ksigma2_0to5_staterr[i] >> netQ_ksigma2_0to5_syserr[i]
			 >> netQ_sigma2byM_70to80[i] >> netQ_sigma2byM_70to80_staterr[i] >> netQ_sigma2byM_70to80_syserr[i]
			 >> netQ_Ssigma_70to80[i] >> netQ_Ssigma_70to80_staterr[i] >> netQ_Ssigma_70to80_syserr[i]
			 >> netQ_ksigma2_70to80[i] >> netQ_ksigma2_70to80_staterr[i] >> netQ_ksigma2_70to80_syserr[i];
		}
	}
	else{
		cout << "Unable to load cumulant/moment ratios...\n";
		n_expData = 0;
	}
	inFileNetP.close();
	inFileNetQ.close();

	double E_l = 5;
	double E_u = 250;
	int intvl = 100;
	double dE = pow(E_u/E_l,1.0/(intvl-1));

	double x;
	double T, mub, muq, mus;

	double snn[intvl];

	double netP_sigma2byM_hrg[intvl];
	double netQ_sigma2byM_hrg[intvl];

	Sigma2byM_Proton sig2byM_netP;
	Sigma2byM_Charge sig2byM_netQ;

	for (int i = 0; i < intvl; ++i)
	{
		x = E_l*pow(dE,i);
		snn[i] = x;

		mub = 1.308/(1+0.273*x);
		muq = -0.0202/(1+0.126*x);
		mus = 0.214/(1+0.184*x);

		T = 0.166-0.139*mub*mub-0.053*pow(mub,4);

		netP_sigma2byM_hrg[i] = sig2byM_netP.getValue(T,mub,muq,mus);
		netQ_sigma2byM_hrg[i] = sig2byM_netQ.getValue(T,mub,muq,mus);
	}
	
	StrangenessNumberDensity nS;
	ChargeToBaryonRatio nQ_by_nB;

	double *minPar, *minPar_err;

	double T_ch[n_expData];
	double mub_ch[n_expData];
	double muq_ch[n_expData];
	double mus_ch[n_expData]; 

	double T_ch_err[n_expData];
	double mub_ch_err[n_expData];
	double muq_ch_err[n_expData];
	double mus_ch_err[n_expData]; 

	cout << "-----------------------------------------------------------------" << '\n';
	cout << "        Determination of chemical freezout parameters using      " << '\n';
	cout << "             net-proton and net-charge fluctuations              " << '\n';
	cout << "-----------------------------------------------------------------" << '\n';
	ThermalMinimizer minimum {4};
	minimum.SetThermalFunctionEach(0,&sig2byM_netP);
	minimum.SetThermalFunctionEach(1,&sig2byM_netQ);
	minimum.SetThermalFunctionEach(2,&nS);
	minimum.SetThermalFunctionEach(3,&nQ_by_nB);

	ofstream outFileFreeze;
	outFileFreeze.open("data/freezout_netQP_fluc.dat");

	bool writeToFile = false;
	if (outFileFreeze)
	{
	 outFileFreeze << n_expData << '\n';
	 writeToFile = true;
	}

	for (int i = 0; i < n_expData; ++i)
	{
		x = snn_exp[i];
		cout << "                    sqrt(snn) = " << x << " GeV" << '\n';

		mub = 1.308/(1+0.273*x);
		muq = -0.0202/(1+0.126*x);
		mus = 0.214/(1+0.184*x);

		T = 0.166-0.139*mub*mub-0.053*pow(mub,4);

		minimum.SetFunctionValueEach(0,netP_sigma2byM_0to5[i]);
		minimum.SetFunctionValueEach(1,netQ_sigma2byM_0to5[i]);
		minimum.SetFunctionValueEach(2,0);
		minimum.SetFunctionValueEach(3,0.4);

		minimum.setInitialParameters(T,mub,muq,mus);
		minimum.minimize(1000);

		minPar = minimum.getMinimizedParameters();
		minPar_err = minimum.getMinimizedParameterErrors();

		T_ch[i] = minPar[0]; mub_ch[i] = minPar[1]; muq_ch[i] = minPar[2]; mus_ch[i] = minPar[3];
		T_ch_err[i] = minPar_err[0]; mub_ch_err[i] = minPar_err[1]; muq_ch_err[i] = minPar_err[2]; mus_ch_err[i] = minPar_err[3];

		if (writeToFile)
		{
			outFileFreeze << minPar[0] << '\t' << minPar_err[0] << '\t' << minPar[1] << '\t' << minPar_err[1]
			 << '\t' << minPar[2] << '\t' << minPar_err[2] << '\t' << minPar[3] << '\t' << minPar_err[3] << '\n';
		}
		
		netP_sigma2byM_fit[i] = sig2byM_netP.getValue(minPar[0],minPar[1],minPar[2],minPar[3]);
		netQ_sigma2byM_fit[i] = sig2byM_netQ.getValue(minPar[0],minPar[1],minPar[2],minPar[3]);
	}
	if (writeToFile)
	{
		outFileFreeze << "\n\\\\ T_ch\tT_ch_err\tmub_ch\tmub_ch_err\tmuq_ch\tmuq_ch_err\tmus_ch\tmus_ch_err\n" ;
		cout << "Freezout parameters saved to data/freezout_netQP_fluc.dat..." << '\n';
		outFileFreeze.close();
	}

	double topMargin = 0.01;
	double bottomMargin = 0.05;
	double rightMargin = 0.02;
	double leftMargin = 0.1;

	double textSize = 0.02;
	double markerSize = 2;
	double lineWidth = 3;
	double transparency = 0.75;

	TCanvas *cl1 = new TCanvas("cl_netP","netP",400,700);

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetEndErrorSize(10);

	double ySizeEach = (1 - topMargin - bottomMargin)/3;

	TPad *p11 = new TPad("p11", "p11",0,2*ySizeEach+bottomMargin,1,1);
	p11->SetLogx();
	p11->SetTickx();
	p11->SetBottomMargin(0);
	p11->SetTopMargin(topMargin/p11->GetHNDC());
	p11->SetRightMargin(rightMargin);
	p11->Draw();

	TPad *p12 = new TPad("p12", "p12",0,ySizeEach+bottomMargin,1,2*ySizeEach+bottomMargin);
	p12->SetLogx();
	p12->SetTickx();
	p12->SetTopMargin(0);
	p12->SetBottomMargin(0);
	p12->SetRightMargin(rightMargin);
	p12->Draw();

	TPad *p13 = new TPad("p13", "p13",0,0,1,ySizeEach+bottomMargin);
	p13->SetLogx();
	p13->SetTickx();
	p13->SetTopMargin(0);
	p13->SetRightMargin(rightMargin);
	p13->SetBottomMargin(bottomMargin/p13->GetHNDC());
	p13->Draw();

	p11->cd();

	auto *gr111 = new TGraph(intvl,snn,netP_sigma2byM_hrg);
	gr111->SetLineStyle(1);
	gr111->SetLineColorAlpha(kGreen+2,transparency);
	gr111->SetLineWidth(lineWidth);

	auto *gr112 = new TGraphErrors(n_expData,snn_exp,netP_sigma2byM_0to5,(double*)0,netP_sigma2byM_0to5_staterr);
	gr112->SetMarkerSize(markerSize);
	gr112->SetMarkerColorAlpha(kOrange-3,transparency);
	gr112->SetLineColorAlpha(kOrange-3,transparency);
	gr112->SetLineWidth(lineWidth);
	gr112->SetMarkerStyle(kFullCircle);
	auto *gr113 = new TGraphErrors(n_expData,snn_exp,netP_sigma2byM_0to5,(double*)0,netP_sigma2byM_0to5_syserr);
	gr113->SetLineWidth(lineWidth);
	gr113->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr114 = new TGraphErrors(n_expData,snn_exp,netP_sigma2byM_70to80,(double*)0,netP_sigma2byM_70to80_staterr);
	gr114->SetMarkerSize(markerSize);
	gr114->SetMarkerColorAlpha(kViolet-1,transparency);
	gr114->SetLineColorAlpha(kViolet-1,transparency);
	gr114->SetLineWidth(lineWidth);
	gr114->SetMarkerStyle(kFullSquare);
	auto *gr115 = new TGraphErrors(n_expData,snn_exp,netP_sigma2byM_70to80,(double*)0,netP_sigma2byM_70to80_syserr);
	gr115->SetLineWidth(lineWidth);
	gr115->SetLineColorAlpha(kViolet-1,transparency);
	auto *gr116 = new TGraph(n_expData,snn_exp,netP_sigma2byM_fit);
	gr116->SetMarkerSize(markerSize);
	gr116->SetMarkerColorAlpha(kBlue,transparency);
	gr116->SetMarkerStyle(kFullTriangleUp);

	TMultiGraph  *mg11  = new TMultiGraph();

	mg11->Add(gr111,"l");
	mg11->Add(gr112,"pz");
	mg11->Add(gr113,"||");
	mg11->Add(gr114,"pz");
	mg11->Add(gr115,"||");
	mg11->Add(gr116,"p");

	// mg11->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg11->GetYaxis()->SetTitle("#sigma^{2}/M");
	// mg11->GetXaxis()->SetRangeUser(2,250);
	mg11->GetYaxis()->SetRangeUser(0.15,9.25);
	mg11->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg11->GetXaxis()->CenterTitle(true);
	mg11->GetYaxis()->CenterTitle(true);
	mg11->GetXaxis()->SetNoExponent();
	mg11->GetXaxis()->SetLabelSize(textSize/p11->GetHNDC());
	mg11->GetXaxis()->SetTitleSize(textSize/p11->GetHNDC());
	mg11->GetXaxis()->SetTitleOffset(0.75);
	mg11->GetYaxis()->SetLabelSize(textSize/p11->GetHNDC());
	mg11->GetYaxis()->SetTitleSize(textSize/p11->GetHNDC());
	mg11->GetYaxis()->SetTickLength(0.02);
	mg11->GetYaxis()->SetTitleOffset(0.65);
	mg11->Draw("a");

	TLatex *text11 = new TLatex(0.4,0.75,"#splitline{net-proton results}{and calculations}");
	text11->SetNDC();
	text11->SetTextSize(textSize/p11->GetHNDC());
	text11->Draw();

	p12->cd();

	auto *gr121 = new TGraphErrors(n_expData,snn_exp,netP_Ssigma_0to5,(double*)0,netP_Ssigma_0to5_staterr);
	gr121->SetMarkerSize(markerSize);
	gr121->SetMarkerColorAlpha(kOrange-3,transparency);
	gr121->SetLineColorAlpha(kOrange-3,transparency);
	gr121->SetLineWidth(lineWidth);
	gr121->SetMarkerStyle(kFullCircle);
	auto *gr122 = new TGraphErrors(n_expData,snn_exp,netP_Ssigma_0to5,(double*)0,netP_Ssigma_0to5_syserr);
	gr122->SetLineWidth(lineWidth);
	gr122->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr123 = new TGraphErrors(n_expData,snn_exp,netP_Ssigma_70to80,(double*)0,netP_Ssigma_70to80_staterr);
	gr123->SetMarkerSize(markerSize);
	gr123->SetMarkerColorAlpha(kViolet-1,transparency);
	gr123->SetLineColorAlpha(kViolet-1,transparency);
	gr123->SetLineWidth(lineWidth);
	gr123->SetMarkerStyle(kFullSquare);
	auto *gr124 = new TGraphErrors(n_expData,snn_exp,netP_Ssigma_70to80,(double*)0,netP_Ssigma_70to80_syserr);
	gr124->SetLineWidth(lineWidth);
	gr124->SetLineColorAlpha(kViolet-1,transparency);

	auto *gr125 = new TGraph(intvl,snn,netP_sigma2byM_hrg);
	gr125->SetLineStyle(1);
	gr125->SetLineColorAlpha(kGreen+2,transparency);
	gr125->SetLineWidth(lineWidth);

	TMultiGraph  *mg12  = new TMultiGraph();

	mg12->Add(gr121,"pz");
	mg12->Add(gr122,"||");
	mg12->Add(gr123,"pz");
	mg12->Add(gr124,"||");
	mg12->Add(gr125,"l");

	// mg12->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg12->GetYaxis()->SetTitle("S#sigma");
	// mg12->GetXaxis()->SetLimits(5,250);
	mg12->GetYaxis()->SetRangeUser(-0.05,1.05);
	mg12->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg12->GetXaxis()->CenterTitle(true);
	mg12->GetYaxis()->CenterTitle(true);
	mg12->GetXaxis()->SetNoExponent();
	mg12->GetXaxis()->SetLabelSize(textSize/p12->GetHNDC());
	mg12->GetXaxis()->SetTitleSize(textSize/p12->GetHNDC());
	mg12->GetXaxis()->SetTitleOffset(0.75);
	mg12->GetYaxis()->SetLabelSize(textSize/p12->GetHNDC());
	mg12->GetYaxis()->SetTitleSize(textSize/p12->GetHNDC());
	mg12->GetYaxis()->SetTickLength(0.02);
	mg12->GetYaxis()->SetTitleOffset(0.65);
	mg12->Draw("a");

	TLatex *latex1 = new TLatex(0.6,0.75,"#splitline{0.4 < p_{T} < 0.8 (GeV/c)}{#left|y#right| #leq 0.5}");
	latex1->SetNDC();
	latex1->SetTextSize(textSize/p12->GetHNDC());
	latex1->Draw();

	TLegend *legend1 = new TLegend(0.12,0.01,0.4,0.35);
	legend1->AddEntry(gr112,"STAR Au-Au 0 - 5%","p");
	legend1->AddEntry(gr114,"STAR Au-Au 70 - 80%","p");
	legend1->AddEntry(gr111,"HRG + decay + cuts","l");
	legend1->AddEntry(gr116,"HRG + decay + cuts (fit)","p");
	legend1->SetTextSize(textSize/p12->GetHNDC());
	legend1->SetBorderSize(0);
	legend1->SetFillStyle(0);
	legend1->Draw();

	p13->cd();

	auto *gr131 = new TGraphErrors(n_expData,snn_exp,netP_ksigma2_0to5,(double*)0,netP_ksigma2_0to5_staterr);
	gr131->SetMarkerSize(markerSize);
	gr131->SetMarkerColorAlpha(kOrange-3,transparency);
	gr131->SetLineColorAlpha(kOrange-3,transparency);
	gr131->SetLineWidth(lineWidth);
	gr131->SetMarkerStyle(kFullCircle);
	auto *gr132 = new TGraphErrors(n_expData,snn_exp,netP_ksigma2_0to5,(double*)0,netP_ksigma2_0to5_syserr);
	gr132->SetLineWidth(lineWidth);
	gr132->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr133 = new TGraphErrors(n_expData,snn_exp,netP_ksigma2_70to80,(double*)0,netP_ksigma2_70to80_staterr);
	gr133->SetMarkerSize(markerSize);
	gr133->SetMarkerColorAlpha(kViolet-1,transparency);
	gr133->SetLineColorAlpha(kViolet-1,transparency);
	gr133->SetLineWidth(lineWidth);
	gr133->SetMarkerStyle(kFullSquare);
	auto *gr134 = new TGraphErrors(n_expData,snn_exp,netP_ksigma2_70to80,(double*)0,netP_ksigma2_70to80_syserr);
	gr134->SetLineWidth(lineWidth);
	gr134->SetLineColorAlpha(kViolet-1,transparency);

	auto *gr135 = new TGraph(intvl,snn,netP_sigma2byM_hrg);
	gr135->SetLineStyle(1);
	gr135->SetLineColorAlpha(kGreen+2,transparency);
	gr135->SetLineWidth(lineWidth);

	TMultiGraph  *mg13  = new TMultiGraph();

	mg13->Add(gr131,"pz");
	mg13->Add(gr132,"||");
	mg13->Add(gr133,"pz");
	mg13->Add(gr134,"||");
	mg13->Add(gr135,"l");

	mg13->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg13->GetYaxis()->SetTitle("k#sigma^{2}");
	mg13->GetYaxis()->SetRangeUser(0.25,1.15);
	mg13->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg13->GetXaxis()->CenterTitle(true);
	mg13->GetYaxis()->CenterTitle(true);
	mg13->GetXaxis()->SetNoExponent();
	mg13->GetXaxis()->SetLabelSize(textSize/p13->GetHNDC());
	mg13->GetXaxis()->SetTitleSize(textSize/p13->GetHNDC());
	mg13->GetXaxis()->SetTitleOffset(0.75);
	mg13->GetYaxis()->SetLabelSize(textSize/p13->GetHNDC());
	mg13->GetYaxis()->SetTitleSize(textSize/p13->GetHNDC());
	mg13->GetYaxis()->SetTickLength(0.02);
	mg13->GetYaxis()->SetTitleOffset(0.65);
	mg13->Draw("a");

	TText *text12 = new TText(0.2,0.9,"STAR Collab., PRL 112, 032302 (2014)");
	text12->SetNDC();
	text12->SetTextSize(textSize/p13->GetHNDC());
	text12->Draw();
	
	TCanvas *cl2 = new TCanvas("cl_netQ","netQ",400,700);

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetEndErrorSize(10);
	
	ySizeEach = (1 - topMargin - bottomMargin)/3;

	TPad *p21 = new TPad("p21", "p21",0,2*ySizeEach+bottomMargin,1,1);
	p21->SetLogx();
	p21->SetTickx();
	p21->SetBottomMargin(0);
	p21->SetTopMargin(topMargin/p21->GetHNDC());
	p21->SetRightMargin(rightMargin);
	p21->Draw();

	TPad *p22 = new TPad("p22", "p22",0,ySizeEach+bottomMargin,1,2*ySizeEach+bottomMargin);
	p22->SetLogx();
	p22->SetTickx();
	p22->SetTopMargin(0);
	p22->SetBottomMargin(0);
	p22->SetRightMargin(rightMargin);
	p22->Draw();

	TPad *p23 = new TPad("p23", "p23",0,0,1,ySizeEach+bottomMargin);
	p23->SetLogx();
	p23->SetTickx();
	p23->SetTopMargin(0);
	p23->SetRightMargin(rightMargin);
	p23->SetBottomMargin(bottomMargin/p23->GetHNDC());
	p23->Draw();

	p21->cd();

	auto *gr211 = new TGraph(intvl,snn,netQ_sigma2byM_hrg);
	gr211->SetLineStyle(1);
	gr211->SetLineColorAlpha(kGreen+2,transparency);
	gr211->SetLineWidth(lineWidth);

	auto *gr212 = new TGraphErrors(n_expData,snn_exp,netQ_sigma2byM_0to5,(double*)0,netQ_sigma2byM_0to5_staterr);
	gr212->SetMarkerSize(markerSize);
	gr212->SetMarkerColorAlpha(kOrange-3,transparency);
	gr212->SetLineColorAlpha(kOrange-3,transparency);
	gr212->SetLineWidth(lineWidth);
	gr212->SetMarkerStyle(kFullCircle);
	auto *gr213 = new TGraphErrors(n_expData,snn_exp,netQ_sigma2byM_0to5,(double*)0,netQ_sigma2byM_0to5_syserr);
	gr213->SetLineWidth(lineWidth);
	gr213->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr214 = new TGraphErrors(n_expData,snn_exp,netQ_sigma2byM_70to80,(double*)0,netQ_sigma2byM_70to80_staterr);
	gr214->SetMarkerSize(markerSize);
	gr214->SetMarkerColorAlpha(kViolet-1,transparency);
	gr214->SetLineColorAlpha(kViolet-1,transparency);
	gr214->SetLineWidth(lineWidth);
	gr214->SetMarkerStyle(kFullSquare);
	auto *gr215 = new TGraphErrors(n_expData,snn_exp,netQ_sigma2byM_70to80,(double*)0,netQ_sigma2byM_70to80_syserr);
	gr215->SetLineWidth(lineWidth);
	gr215->SetLineColorAlpha(kViolet-1,transparency);
	auto *gr216 = new TGraph(n_expData,snn_exp,netQ_sigma2byM_fit);
	gr216->SetMarkerSize(markerSize);
	gr216->SetMarkerColorAlpha(kBlue,transparency);
	gr216->SetMarkerStyle(kFullTriangleUp);

	TMultiGraph  *mg21  = new TMultiGraph();

	mg21->Add(gr211,"l");
	mg21->Add(gr212,"pz");
	mg21->Add(gr213,"||");
	mg21->Add(gr214,"pz");
	mg21->Add(gr215,"||");
	mg21->Add(gr216,"p");

	// mg21->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg21->GetYaxis()->SetTitle("#sigma^{2}/M");
	// mg21->GetXaxis()->SetRangeUser(2,250);
	mg21->GetYaxis()->SetRangeUser(-7.5,112.5);
	mg21->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg21->GetXaxis()->CenterTitle(true);
	mg21->GetYaxis()->CenterTitle(true);
	mg21->GetXaxis()->SetNoExponent();
	mg21->GetXaxis()->SetLabelSize(textSize/p21->GetHNDC());
	mg21->GetXaxis()->SetTitleSize(textSize/p21->GetHNDC());
	mg21->GetXaxis()->SetTitleOffset(0.75);
	mg21->GetYaxis()->SetLabelSize(textSize/p21->GetHNDC());
	mg21->GetYaxis()->SetTitleSize(textSize/p21->GetHNDC());
	mg21->GetYaxis()->SetTickLength(0.02);
	mg21->GetYaxis()->SetTitleOffset(0.65);
	mg21->Draw("a");

	TLatex *text221 = new TLatex(0.4,0.75,"#splitline{net-charge results}{and calculations}");
	text221->SetNDC();
	text221->SetTextSize(textSize/p21->GetHNDC());
	text221->Draw();

	p22->cd();

	auto *gr221 = new TGraphErrors(n_expData,snn_exp,netQ_Ssigma_0to5,(double*)0,netQ_Ssigma_0to5_staterr);
	gr221->SetMarkerSize(markerSize);
	gr221->SetMarkerColorAlpha(kOrange-3,transparency);
	gr221->SetLineColorAlpha(kOrange-3,transparency);
	gr221->SetLineWidth(lineWidth);
	gr221->SetMarkerStyle(kFullCircle);
	auto *gr222 = new TGraphErrors(n_expData,snn_exp,netQ_Ssigma_0to5,(double*)0,netQ_Ssigma_0to5_syserr);
	gr222->SetLineWidth(lineWidth);
	gr222->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr223 = new TGraphErrors(n_expData,snn_exp,netQ_Ssigma_70to80,(double*)0,netQ_Ssigma_70to80_staterr);
	gr223->SetMarkerSize(markerSize);
	gr223->SetMarkerColorAlpha(kViolet-1,transparency);
	gr223->SetLineColorAlpha(kViolet-1,transparency);
	gr223->SetLineWidth(lineWidth);
	gr223->SetMarkerStyle(kFullSquare);
	auto *gr224 = new TGraphErrors(n_expData,snn_exp,netQ_Ssigma_70to80,(double*)0,netQ_Ssigma_70to80_syserr);
	gr224->SetLineWidth(lineWidth);
	gr224->SetLineColorAlpha(kViolet-1,transparency);

	auto *gr225 = new TGraph(intvl,snn,netQ_sigma2byM_hrg);
	gr225->SetLineStyle(1);
	gr225->SetLineColorAlpha(kGreen+2,transparency);
	gr225->SetLineWidth(lineWidth);

	TMultiGraph  *mg22  = new TMultiGraph();

	mg22->Add(gr221,"pz");
	mg22->Add(gr222,"||");
	mg22->Add(gr223,"pz");
	mg22->Add(gr224,"||");
	mg22->Add(gr225,"l");

	// mg22->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg22->GetYaxis()->SetTitle("S#sigma");
	// mg22->GetXaxis()->SetLimits(5,250);
	mg22->GetYaxis()->SetRangeUser(-0.095,0.75);
	mg22->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg22->GetXaxis()->CenterTitle(true);
	mg22->GetYaxis()->CenterTitle(true);
	mg22->GetXaxis()->SetNoExponent();
	mg22->GetXaxis()->SetLabelSize(textSize/p22->GetHNDC());
	mg22->GetXaxis()->SetTitleSize(textSize/p22->GetHNDC());
	mg22->GetXaxis()->SetTitleOffset(0.75);
	mg22->GetYaxis()->SetLabelSize(textSize/p22->GetHNDC());
	mg22->GetYaxis()->SetTitleSize(textSize/p22->GetHNDC());
	mg22->GetYaxis()->SetTickLength(0.02);
	mg22->GetYaxis()->SetTitleOffset(0.65);
	mg22->Draw("a");

	TLegend *legend2 = new TLegend(0.5,0.6,0.85,0.85);
	legend2->AddEntry(gr212,"STAR Au-Au 0 - 5%","p");
	legend2->AddEntry(gr214,"STAR Au-Au 70 - 80%","p");
	legend2->AddEntry(gr211,"HRG + decay + cuts","l");
	legend2->AddEntry(gr216,"HRG + decay + cuts (fit)","p");
	legend2->SetTextSize(textSize/p22->GetHNDC());
	legend2->SetBorderSize(0);
	legend2->SetFillStyle(0);
	legend2->Draw();	

	p23->cd();

	auto *gr231 = new TGraphErrors(n_expData,snn_exp,netQ_ksigma2_0to5,(double*)0,netQ_ksigma2_0to5_staterr);
	gr231->SetMarkerSize(markerSize);
	gr231->SetMarkerColorAlpha(kOrange-3,transparency);
	gr231->SetLineColorAlpha(kOrange-3,transparency);
	gr231->SetLineWidth(lineWidth);
	gr231->SetMarkerStyle(kFullCircle);
	auto *gr232 = new TGraphErrors(n_expData,snn_exp,netQ_ksigma2_0to5,(double*)0,netQ_ksigma2_0to5_syserr);
	gr232->SetLineWidth(lineWidth);
	gr232->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr233 = new TGraphErrors(n_expData,snn_exp,netQ_ksigma2_70to80,(double*)0,netQ_ksigma2_70to80_staterr);
	gr233->SetMarkerSize(markerSize);
	gr233->SetMarkerColorAlpha(kViolet-1,transparency);
	gr233->SetLineColorAlpha(kViolet-1,transparency);
	gr233->SetLineWidth(lineWidth);
	gr233->SetMarkerStyle(kFullSquare);
	auto *gr234 = new TGraphErrors(n_expData,snn_exp,netQ_ksigma2_70to80,(double*)0,netQ_ksigma2_70to80_syserr);
	gr234->SetLineWidth(lineWidth);
	gr234->SetLineColorAlpha(kViolet-1,transparency);

	auto *gr235 = new TGraph(intvl,snn,netQ_sigma2byM_hrg);
	gr235->SetLineStyle(1);
	gr235->SetLineColorAlpha(kGreen+2,transparency);
	gr235->SetLineWidth(lineWidth);

	TMultiGraph  *mg23  = new TMultiGraph();

	mg23->Add(gr231,"pz");
	mg23->Add(gr232,"||");
	mg23->Add(gr233,"pz");
	mg23->Add(gr234,"||");
	mg23->Add(gr235,"l");

	mg23->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg23->GetYaxis()->SetTitle("k#sigma^{2}");
	mg23->GetYaxis()->SetRangeUser(-19.5,12.5);
	mg23->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg23->GetXaxis()->CenterTitle(true);
	mg23->GetYaxis()->CenterTitle(true);
	mg23->GetXaxis()->SetNoExponent();
	mg23->GetXaxis()->SetLabelSize(textSize/p23->GetHNDC());
	mg23->GetXaxis()->SetTitleSize(textSize/p23->GetHNDC());
	mg23->GetXaxis()->SetTitleOffset(0.75);
	mg23->GetYaxis()->SetLabelSize(textSize/p23->GetHNDC());
	mg23->GetYaxis()->SetTitleSize(textSize/p23->GetHNDC());
	mg23->GetYaxis()->SetTickLength(0.02);
	mg23->GetYaxis()->SetTitleOffset(0.65);
	mg23->Draw("a");

	TLatex *latex2 = new TLatex(0.5,0.35,"#splitline{0.2 < p_{T} < 2.0 (GeV/c)}{#left|#eta#right| #leq 0.5}");
	latex2->SetNDC();
	latex2->SetTextSize(textSize/p22->GetHNDC());
	latex2->Draw();

	TText *text22 = new TText(0.32,0.9,"STAR Collab., PRL 113, 092301 (2014)");
	text22->SetNDC();
	text22->SetTextSize(textSize/p23->GetHNDC());
	text22->Draw();


	TCanvas *cl3 = new TCanvas("cl_freeze_T_mub","freeze_T_mub",550,400);

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetEndErrorSize(10);

	TPad *p3 = new TPad("p3", "p3",0,0,1,1);
	p3->SetBottomMargin(0.12);
	p3->SetTopMargin(0.01);
	p3->SetRightMargin(0.01);
	p3->SetLeftMargin(0.11);
	p3->Draw();

	p3->cd();

	double mub_l = 0, mub_u = 1;
	double T_l = 0.1, T_u = 0.18;

	intvl = 100;

	double dT = (T_u-T_l)/(intvl-1);
	double dmub = (mub_u-mub_l)/(intvl-1);

	double T_ch_old[intvl];
	double mub_ch_old[intvl];

	for (int i = 0; i < intvl; ++i)
	{
		x = mub_l+i*dmub;
		mub_ch_old[i] = x;
		T_ch_old[i] = 0.166-0.139*x*x-0.053*pow(x,4);
	}

	auto *gr31 = new TGraphErrors(n_expData,mub_ch,T_ch,mub_ch_err,T_ch_err);
	gr31->SetMarkerSize(markerSize);
	gr31->SetMarkerColorAlpha(kOrange-3,transparency);
	gr31->SetLineColorAlpha(kOrange-3,transparency);
	gr31->SetLineWidth(lineWidth);
	gr31->SetMarkerStyle(kFullCircle);

	auto *gr32 = new TGraph(intvl,mub_ch_old,T_ch_old);
	gr32->SetLineColorAlpha(kGreen+2,transparency);
	gr32->SetLineWidth(lineWidth);

	TMultiGraph  *mg3  = new TMultiGraph();

	mg3->Add(gr31,"pz");
	mg3->Add(gr32,"l");

	textSize = 0.05;

	mg3->GetXaxis()->SetTitle("#mu_{B} (GeV)");
	mg3->GetYaxis()->SetTitle("T (GeV)");
	mg3->GetXaxis()->SetRangeUser(0.03,0.67);
	mg3->GetYaxis()->SetRangeUser(0.007,0.183);
	mg3->GetXaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg3->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg3->GetXaxis()->CenterTitle(true);
	mg3->GetYaxis()->CenterTitle(true);

	mg3->GetXaxis()->SetLabelSize(textSize);
	mg3->GetXaxis()->SetTitleSize(textSize);
	mg3->GetXaxis()->SetTickLength(0.02);
	mg3->GetXaxis()->SetTitleOffset(1);

	mg3->GetYaxis()->SetLabelSize(textSize);
	mg3->GetYaxis()->SetTitleSize(textSize);
	mg3->GetYaxis()->SetTickLength(0.02);
	mg3->GetYaxis()->SetTitleOffset(1);
	mg3->Draw("a");

	TLegend *legend3 = new TLegend(0.15,0.2,0.6,0.35);
	legend3->AddEntry(gr31,"Chem. freezout (using fluc.)","p");
	legend3->AddEntry(gr32,"Chem. freezout (using yield)","l");
	legend3->SetTextSize(textSize);
	legend3->SetBorderSize(0);
	legend3->SetFillStyle(0);
	legend3->Draw();
}

void generateCumulantRatios_netP()
{
	int n_expData;

	double *snn_exp;

	double *C1_0to5;
	double *C1_0to5_staterr;
	double *C1_0to5_syserr;

	double *C2_0to5;
	double *C2_0to5_staterr;
	double *C2_0to5_syserr;

	double *C3_0to5;
	double *C3_0to5_staterr;
	double *C3_0to5_syserr;

	double *C4_0to5;
	double *C4_0to5_staterr;
	double *C4_0to5_syserr;

	double *C1_70to80;
	double *C1_70to80_staterr;
	double *C1_70to80_syserr;

	double *C2_70to80;
	double *C2_70to80_staterr;
	double *C2_70to80_syserr;

	double *C3_70to80;
	double *C3_70to80_staterr;
	double *C3_70to80_syserr;

	double *C4_70to80;
	double *C4_70to80_staterr;
	double *C4_70to80_syserr;

	double *sigma2byM_0to5;
	double *sigma2byM_0to5_staterr;
	double *sigma2byM_0to5_syserr;

	double *Ssigma_0to5;
	double *Ssigma_0to5_syserr;
	double *Ssigma_0to5_staterr;

	double *ksigma2_0to5;
	double *ksigma2_0to5_syserr;
	double *ksigma2_0to5_staterr;

	double *sigma2byM_70to80;
	double *sigma2byM_70to80_staterr;
	double *sigma2byM_70to80_syserr;

	double *Ssigma_70to80;
	double *Ssigma_70to80_syserr;
	double *Ssigma_70to80_staterr;

	double *ksigma2_70to80;
	double *ksigma2_70to80_syserr;
	double *ksigma2_70to80_staterr;

	ifstream inFile;
	inFile.open("data/cumulants_netP.dat");
	ofstream outFile;
	outFile.open("data/cumulantRatios_netP.dat");
	if (inFile && outFile)
	{
		int l;
		inFile >> l;
		outFile << l << '\n';
		n_expData = l;

		snn_exp = new double[l];

		C1_0to5 = new double[l];
		C1_0to5_staterr = new double[l];
		C1_0to5_syserr = new double[l];

		C2_0to5 = new double[l];
		C2_0to5_staterr = new double[l];
		C2_0to5_syserr = new double[l];

		C3_0to5 = new double[l];
		C3_0to5_staterr = new double[l];
		C3_0to5_syserr = new double[l];

		C4_0to5 = new double[l];
		C4_0to5_staterr = new double[l];
		C4_0to5_syserr = new double[l];

		C1_70to80 = new double[l];
		C1_70to80_staterr = new double[l];
		C1_70to80_syserr = new double[l];

		C2_70to80 = new double[l];
		C2_70to80_staterr = new double[l];
		C2_70to80_syserr = new double[l];

		C3_70to80 = new double[l];
		C3_70to80_staterr = new double[l];
		C3_70to80_syserr = new double[l];

		C4_70to80 = new double[l];
		C4_70to80_staterr = new double[l];
		C4_70to80_syserr = new double[l];

		sigma2byM_0to5 = new double[l];
		sigma2byM_0to5_syserr = new double[l];
		sigma2byM_0to5_staterr = new double[l];

		Ssigma_0to5 = new double[l];
		Ssigma_0to5_syserr = new double[l];
		Ssigma_0to5_staterr = new double[l];

		ksigma2_0to5 = new double[l];
		ksigma2_0to5_syserr = new double[l];
		ksigma2_0to5_staterr = new double[l];

		sigma2byM_70to80 = new double[l];
		sigma2byM_70to80_syserr = new double[l];
		sigma2byM_70to80_staterr = new double[l];

		Ssigma_70to80 = new double[l];
		Ssigma_70to80_syserr = new double[l];
		Ssigma_70to80_staterr = new double[l];

		ksigma2_70to80 = new double[l];
		ksigma2_70to80_syserr = new double[l];
		ksigma2_70to80_staterr = new double[l];

		for (int i = 0; i < l; ++i)
		{
			inFile >> snn_exp[i] >> C1_0to5[i] >> C1_0to5_staterr[i] >> C1_0to5_syserr[i]
				>> C2_0to5[i] >> C2_0to5_staterr[i] >> C2_0to5_syserr[i]
				>> C3_0to5[i] >> C3_0to5_staterr[i] >> C3_0to5_syserr[i]
				>> C4_0to5[i] >> C4_0to5_staterr[i] >> C4_0to5_syserr[i]
				>> C1_70to80[i] >> C1_70to80_staterr[i] >> C1_70to80_syserr[i]
				>> C2_70to80[i] >> C2_70to80_staterr[i] >> C2_70to80_syserr[i]
				>> C3_70to80[i] >> C3_70to80_staterr[i] >> C3_70to80_syserr[i]
				>> C4_70to80[i] >> C4_70to80_staterr[i] >> C4_70to80_syserr[i];
			// cout << snn_exp[i] << '\t' << C1_0to5[i] << '\t' << C1_0to5_staterr[i] << '\t' << C1_0to5_syserr[i]
			// 	<< '\t' << C2_0to5[i] << '\t' << C2_0to5_staterr[i] << '\t' << C2_0to5_syserr[i]
			// 	<< '\t' << C3_0to5[i] << '\t' << C3_0to5_staterr[i] << '\t' << C3_0to5_syserr[i]
			// 	<< '\t' << C4_0to5[i] << '\t' << C4_0to5_staterr[i] << '\t' << C4_0to5_syserr[i]
			// 	<< '\t' << C1_70to80[i] << '\t' << C1_70to80_staterr[i] << '\t' << C1_70to80_syserr[i]
			// 	<< '\t' << C2_70to80[i] << '\t' << C2_70to80_staterr[i] << '\t' << C2_70to80_syserr[i]
			// 	<< '\t' << C3_70to80[i] << '\t' << C3_70to80_staterr[i] << '\t' << C3_70to80_syserr[i]
			// 	<< '\t' << C4_70to80[i] << '\t' << C4_70to80_staterr[i] << '\t' << C4_70to80_syserr[i] << '\n';

			sigma2byM_0to5[i] = C2_0to5[i]/C1_0to5[i];
			sigma2byM_0to5_staterr[i] = sigma2byM_0to5[i]*sqrt(pow(C1_0to5_staterr[i]/C1_0to5[i],2)+pow(C2_0to5_staterr[i]/C2_0to5[i],2));
			sigma2byM_0to5_syserr[i] = sigma2byM_0to5[i]*sqrt(pow(C1_0to5_syserr[i]/C1_0to5[i],2)+pow(C2_0to5_syserr[i]/C2_0to5[i],2));

			Ssigma_0to5[i] = C3_0to5[i]/C2_0to5[i];
			Ssigma_0to5_staterr[i] = Ssigma_0to5[i]*sqrt(pow(C3_0to5_staterr[i]/C3_0to5[i],2)+pow(C2_0to5_staterr[i]/C2_0to5[i],2));
			Ssigma_0to5_syserr[i] = Ssigma_0to5[i]*sqrt(pow(C3_0to5_syserr[i]/C3_0to5[i],2)+pow(C2_0to5_syserr[i]/C2_0to5[i],2));

			ksigma2_0to5[i] = C4_0to5[i]/C2_0to5[i];
			ksigma2_0to5_staterr[i] = ksigma2_0to5[i]*sqrt(pow(C4_0to5_staterr[i]/C4_0to5[i],2)+pow(C2_0to5_staterr[i]/C2_0to5[i],2));
			ksigma2_0to5_syserr[i] = ksigma2_0to5[i]*sqrt(pow(C4_0to5_syserr[i]/C4_0to5[i],2)+pow(C2_0to5_syserr[i]/C2_0to5[i],2));

			sigma2byM_70to80[i] = C2_70to80[i]/C1_70to80[i];
			sigma2byM_70to80_staterr[i] = sigma2byM_70to80[i]*sqrt(pow(C1_70to80_staterr[i]/C1_70to80[i],2)+pow(C2_70to80_staterr[i]/C2_70to80[i],2));
			sigma2byM_70to80_syserr[i] = sigma2byM_70to80[i]*sqrt(pow(C1_70to80_syserr[i]/C1_70to80[i],2)+pow(C2_70to80_syserr[i]/C2_70to80[i],2));

			Ssigma_70to80[i] = C3_70to80[i]/C2_70to80[i];
			Ssigma_70to80_staterr[i] = Ssigma_70to80[i]*sqrt(pow(C3_70to80_staterr[i]/C3_70to80[i],2)+pow(C2_70to80_staterr[i]/C2_70to80[i],2));
			Ssigma_70to80_syserr[i] = Ssigma_70to80[i]*sqrt(pow(C3_70to80_syserr[i]/C3_70to80[i],2)+pow(C2_70to80_syserr[i]/C2_70to80[i],2));

			ksigma2_70to80[i] = C4_70to80[i]/C2_70to80[i];
			ksigma2_70to80_staterr[i] = ksigma2_70to80[i]*sqrt(pow(C4_70to80_staterr[i]/C4_70to80[i],2)+pow(C2_70to80_staterr[i]/C2_70to80[i],2));
			ksigma2_70to80_syserr[i] = ksigma2_70to80[i]*sqrt(pow(C4_70to80_syserr[i]/C4_70to80[i],2)+pow(C2_70to80_syserr[i]/C2_70to80[i],2));

			outFile << snn_exp[i] << '\t' << sigma2byM_0to5[i] << '\t' << sigma2byM_0to5_staterr[i] << '\t' << sigma2byM_0to5_syserr[i]
			 << '\t' << Ssigma_0to5[i] << '\t' << Ssigma_0to5_staterr[i] << '\t' << Ssigma_0to5_syserr[i]
			 << '\t' << ksigma2_0to5[i] << '\t' << ksigma2_0to5_staterr[i] << '\t' << ksigma2_0to5_syserr[i]
			 << '\t' << sigma2byM_70to80[i] << '\t' << sigma2byM_70to80_staterr[i] << '\t' << sigma2byM_70to80_syserr[i]
			 << '\t' << Ssigma_70to80[i] << '\t' << Ssigma_70to80_staterr[i] << '\t' << Ssigma_70to80_syserr[i]
			 << '\t' << ksigma2_70to80[i] << '\t' << ksigma2_70to80_staterr[i] << '\t' << ksigma2_70to80_syserr[i] << '\n';
		}
		outFile << "\n\\\\ snn_exp\tsigma2byM_0to5\tsigma2byM_0to5_staterr\tsigma2byM_0to5_syserr\tSsigma_0to5\tSsigma_0to5_staterr\tSsigma_0to5_syserr\tksigma2_0to5\tksigma2_0to5_staterr\tksigma2_0to5_syserr\tsigma2byM_70to80\tsigma2byM_70to80_staterr\tsigma2byM_70to80_syserr\tSsigma_70to80\tSsigma_70to80_staterr\tSsigma_70to80_syserr\tksigma2_70to80\tksigma2_70to80_staterr\tksigma2_70to80_syserr\n";
		cout << "net-proton cumulant ratios saved in data/cumulantRatios_netP.dat...\n" ;
	}
	else{
		cout << "Unable to load input/output file...\n" ;
	}
	inFile.close();
	outFile.close();

	delete[] snn_exp;

	delete[] C1_0to5; delete[] C1_0to5_staterr; delete[] C1_0to5_syserr;
	delete[] C2_0to5; delete[] C2_0to5_staterr; delete[] C2_0to5_syserr;
	delete[] C3_0to5; delete[] C3_0to5_staterr; delete[] C3_0to5_syserr;
	delete[] C4_0to5; delete[] C4_0to5_staterr; delete[] C4_0to5_syserr;

	delete[] C1_70to80; delete[] C1_70to80_staterr; delete[] C1_70to80_syserr;
	delete[] C2_70to80; delete[] C2_70to80_staterr; delete[] C2_70to80_syserr;
	delete[] C3_70to80; delete[] C3_70to80_staterr; delete[] C3_70to80_syserr;
	delete[] C4_70to80; delete[] C4_70to80_staterr; delete[] C4_70to80_syserr;

	delete[] sigma2byM_0to5; delete[] sigma2byM_0to5_staterr; delete[] sigma2byM_0to5_syserr;
	delete[] Ssigma_0to5; delete[] Ssigma_0to5_syserr; delete[] Ssigma_0to5_staterr;
	delete[] ksigma2_0to5; delete[] ksigma2_0to5_syserr; delete[] ksigma2_0to5_staterr;

	delete[] sigma2byM_70to80; delete[] sigma2byM_70to80_staterr; delete[] sigma2byM_70to80_syserr;
	delete[] Ssigma_70to80; delete[] Ssigma_70to80_syserr; delete[] Ssigma_70to80_staterr;
	delete[] ksigma2_70to80; delete[] ksigma2_70to80_syserr; delete[] ksigma2_70to80_staterr;
}

void plotSimilarToSTARnetProton(){

	double pT_min = 0.4;		//kinematic cuts
	double pT_max = 0.8;
	double y_min = -0.5;		//rapidity
	double y_max = 0.5;

	double E_l = 5;
	double E_u = 250;
	int intvl = 100;
	double dE = pow(E_u/E_l,1.0/(intvl-1));

	double x,y;
	double num, denom;
	double T, mub, muq, mus;

	ParticleSystem system {};
	system.loadDefaultData();

	Susceptibility1 chi1 {system};
	Susceptibility2 chi2 {system};

	chi1.setDecay(true);
	chi2.setDecay(true);

	chi1.setDetectorRapidityCoordinates();
	chi2.setDetectorRapidityCoordinates();

	chi1.setIntegralLimits(pT_min, pT_max, y_min, y_max);
	chi2.setIntegralLimits(pT_min, pT_max, y_min, y_max);

	double snn[intvl];
	double sigma2byM_p_hrg[intvl];

	int n_expData;

	double *snn_exp;

	double *sigma2byM_0to5;
	double *sigma2byM_0to5_staterr;
	double *sigma2byM_0to5_syserr;

	double *Ssigma_0to5;
	double *Ssigma_0to5_syserr;
	double *Ssigma_0to5_staterr;

	double *ksigma2_0to5;
	double *ksigma2_0to5_syserr;
	double *ksigma2_0to5_staterr;

	double *sigma2byM_70to80;
	double *sigma2byM_70to80_staterr;
	double *sigma2byM_70to80_syserr;

	double *Ssigma_70to80;
	double *Ssigma_70to80_syserr;
	double *Ssigma_70to80_staterr;

	double *ksigma2_70to80;
	double *ksigma2_70to80_syserr;
	double *ksigma2_70to80_staterr;

	for (int i = 0; i < intvl; ++i)
	{
		x = E_l*pow(dE,i);
		snn[i] = x;

		mub = 1.308/(1+0.273*x);
		muq = -0.0202/(1+0.126*x);
		mus = 0.214/(1+0.184*x);

		T = 0.166-0.139*mub*mub-0.053*pow(mub,4);

		num = chi2.getValueEach(23,T,mub,muq,mus)+chi2.getValueEach(24,T,mub,muq,mus);
		denom = chi1.getValueEach(23,T,mub,muq,mus)-chi1.getValueEach(24,T,mub,muq,mus);
		y = num/denom;
		sigma2byM_p_hrg[i] = y;
	}

	ifstream inFile;
	inFile.open("data/cumulantRatios_netP.dat");
	if (inFile)
	{
		int l;
		inFile >> l;
		n_expData = l;

		snn_exp = new double[l];

		sigma2byM_0to5 = new double[l];
		sigma2byM_0to5_syserr = new double[l];
		sigma2byM_0to5_staterr = new double[l];

		Ssigma_0to5 = new double[l];
		Ssigma_0to5_syserr = new double[l];
		Ssigma_0to5_staterr = new double[l];

		ksigma2_0to5 = new double[l];
		ksigma2_0to5_syserr = new double[l];
		ksigma2_0to5_staterr = new double[l];

		sigma2byM_70to80 = new double[l];
		sigma2byM_70to80_syserr = new double[l];
		sigma2byM_70to80_staterr = new double[l];

		Ssigma_70to80 = new double[l];
		Ssigma_70to80_syserr = new double[l];
		Ssigma_70to80_staterr = new double[l];

		ksigma2_70to80 = new double[l];
		ksigma2_70to80_syserr = new double[l];
		ksigma2_70to80_staterr = new double[l];

		for (int i = 0; i < l; ++i)
		{
			inFile >> snn_exp[i] >> sigma2byM_0to5[i] >> sigma2byM_0to5_staterr[i] >> sigma2byM_0to5_syserr[i]
			 >> Ssigma_0to5[i] >> Ssigma_0to5_staterr[i] >> Ssigma_0to5_syserr[i]
			 >> ksigma2_0to5[i] >> ksigma2_0to5_staterr[i] >> ksigma2_0to5_syserr[i]
			 >> sigma2byM_70to80[i] >> sigma2byM_70to80_staterr[i] >> sigma2byM_70to80_syserr[i]
			 >> Ssigma_70to80[i] >> Ssigma_70to80_staterr[i] >> Ssigma_70to80_syserr[i]
			 >> ksigma2_70to80[i] >> ksigma2_70to80_staterr[i] >> ksigma2_70to80_syserr[i];
		}
	}
	else{
		cout << "Unable to load net-proton cumulant ratios...\n";
		n_expData = 0;
	}
	inFile.close();

	double topMargin = 0.01;
	double bottomMargin = 0.05;
	double rightMargin = 0.02;
	double leftMargin = 0.1;

	double textSize = 0.02;
	double markerSize = 2;
	double lineWidth = 3;
	double transparency = 0.75;

	TCanvas *cl1 = new TCanvas("cl_netP","netP",400,700);

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetEndErrorSize(10);

	double ySizeEach = (1 - topMargin - bottomMargin)/3;

	TPad *p1 = new TPad("p1", "p1",0,2*ySizeEach+bottomMargin,1,1);
	p1->SetLogx();
	p1->SetTickx();
	p1->SetBottomMargin(0);
	p1->SetTopMargin(topMargin/p1->GetHNDC());
	p1->SetRightMargin(rightMargin);
	p1->Draw();

	TPad *p2 = new TPad("p2", "p2",0,ySizeEach+bottomMargin,1,2*ySizeEach+bottomMargin);
	p2->SetLogx();
	p2->SetTickx();
	p2->SetTopMargin(0);
	p2->SetBottomMargin(0);
	p2->SetRightMargin(rightMargin);
	p2->Draw();

	TPad *p3 = new TPad("p3", "p3",0,0,1,ySizeEach+bottomMargin);
	p3->SetLogx();
	p3->SetTickx();
	p3->SetTopMargin(0);
	p3->SetRightMargin(rightMargin);
	p3->SetBottomMargin(bottomMargin/p3->GetHNDC());
	p3->Draw();

	p1->cd();

	auto *gr11 = new TGraphErrors(n_expData,snn_exp,sigma2byM_0to5,(double*)0,sigma2byM_0to5_staterr);
	gr11->SetMarkerSize(markerSize);
	gr11->SetMarkerColorAlpha(kOrange-3,transparency);
	gr11->SetLineColorAlpha(kOrange-3,transparency);
	gr11->SetLineWidth(lineWidth);
	gr11->SetMarkerStyle(kFullCircle);
	auto *gr12 = new TGraphErrors(n_expData,snn_exp,sigma2byM_0to5,(double*)0,sigma2byM_0to5_syserr);
	gr12->SetLineWidth(lineWidth);
	gr12->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr13 = new TGraphErrors(n_expData,snn_exp,sigma2byM_70to80,(double*)0,sigma2byM_70to80_staterr);
	gr13->SetMarkerSize(markerSize);
	gr13->SetMarkerColorAlpha(kViolet-1,transparency);
	gr13->SetLineColorAlpha(kViolet-1,transparency);
	gr13->SetLineWidth(lineWidth);
	gr13->SetMarkerStyle(kFullSquare);
	auto *gr14 = new TGraphErrors(n_expData,snn_exp,sigma2byM_70to80,(double*)0,sigma2byM_70to80_syserr);
	gr14->SetLineWidth(lineWidth);
	gr14->SetLineColorAlpha(kViolet-1,transparency);

	auto *gr15 = new TGraph(intvl,snn,sigma2byM_p_hrg);
	gr15->SetLineStyle(1);
	gr15->SetLineColorAlpha(kGreen+2,transparency);
	gr15->SetLineWidth(lineWidth);

	TMultiGraph  *mg1  = new TMultiGraph();

	mg1->Add(gr11,"pz");
	mg1->Add(gr12,"||");
	mg1->Add(gr13,"pz");
	mg1->Add(gr14,"||");
	mg1->Add(gr15,"l");

	// mg1->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg1->GetYaxis()->SetTitle("#sigma^{2}/M");
	// mg1->GetXaxis()->SetRangeUser(2,250);
	mg1->GetYaxis()->SetRangeUser(0.15,9.25);
	mg1->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg1->GetXaxis()->CenterTitle(true);
	mg1->GetYaxis()->CenterTitle(true);
	mg1->GetXaxis()->SetNoExponent();
	mg1->GetXaxis()->SetLabelSize(textSize/p1->GetHNDC());
	mg1->GetXaxis()->SetTitleSize(textSize/p1->GetHNDC());
	mg1->GetXaxis()->SetTitleOffset(0.75);
	mg1->GetYaxis()->SetLabelSize(textSize/p1->GetHNDC());
	mg1->GetYaxis()->SetTitleSize(textSize/p1->GetHNDC());
	mg1->GetYaxis()->SetTickLength(0.02);
	mg1->GetYaxis()->SetTitleOffset(0.65);
	mg1->Draw("a");

	TLatex *text1 = new TLatex(0.4,0.75,"#splitline{net-proton results}{and calculations}");
	text1->SetNDC();
	text1->SetTextSize(textSize/p1->GetHNDC());
	text1->Draw();

	p2->cd();

	auto *gr21 = new TGraphErrors(n_expData,snn_exp,Ssigma_0to5,(double*)0,Ssigma_0to5_staterr);
	gr21->SetMarkerSize(markerSize);
	gr21->SetMarkerColorAlpha(kOrange-3,transparency);
	gr21->SetLineColorAlpha(kOrange-3,transparency);
	gr21->SetLineWidth(lineWidth);
	gr21->SetMarkerStyle(kFullCircle);
	auto *gr22 = new TGraphErrors(n_expData,snn_exp,Ssigma_0to5,(double*)0,Ssigma_0to5_syserr);
	gr22->SetLineWidth(lineWidth);
	gr22->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr23 = new TGraphErrors(n_expData,snn_exp,Ssigma_70to80,(double*)0,Ssigma_70to80_staterr);
	gr23->SetMarkerSize(markerSize);
	gr23->SetMarkerColorAlpha(kViolet-1,transparency);
	gr23->SetLineColorAlpha(kViolet-1,transparency);
	gr23->SetLineWidth(lineWidth);
	gr23->SetMarkerStyle(kFullSquare);
	auto *gr24 = new TGraphErrors(n_expData,snn_exp,Ssigma_70to80,(double*)0,Ssigma_70to80_syserr);
	gr24->SetLineWidth(lineWidth);
	gr24->SetLineColorAlpha(kViolet-1,transparency);

	auto *gr25 = new TGraph(intvl,snn,sigma2byM_p_hrg);
	gr25->SetLineStyle(1);
	gr25->SetLineColorAlpha(kGreen+2,transparency);
	gr25->SetLineWidth(lineWidth);

	TMultiGraph  *mg2  = new TMultiGraph();

	mg2->Add(gr21,"pz");
	mg2->Add(gr22,"||");
	mg2->Add(gr23,"pz");
	mg2->Add(gr24,"||");
	mg2->Add(gr25,"l");

	// mg2->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg2->GetYaxis()->SetTitle("S#sigma");
	// mg2->GetXaxis()->SetLimits(5,250);
	mg2->GetYaxis()->SetRangeUser(-0.05,1.05);
	mg2->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg2->GetXaxis()->CenterTitle(true);
	mg2->GetYaxis()->CenterTitle(true);
	mg2->GetXaxis()->SetNoExponent();
	mg2->GetXaxis()->SetLabelSize(textSize/p2->GetHNDC());
	mg2->GetXaxis()->SetTitleSize(textSize/p2->GetHNDC());
	mg2->GetXaxis()->SetTitleOffset(0.75);
	mg2->GetYaxis()->SetLabelSize(textSize/p2->GetHNDC());
	mg2->GetYaxis()->SetTitleSize(textSize/p2->GetHNDC());
	mg2->GetYaxis()->SetTickLength(0.02);
	mg2->GetYaxis()->SetTitleOffset(0.65);
	mg2->Draw("a");

	TLatex *latex = new TLatex(0.6,0.75,"#splitline{0.4 < p_{T} < 0.8 (GeV/c)}{#left|y#right| #leq 0.5}");
	latex->SetNDC();
	latex->SetTextSize(textSize/p2->GetHNDC());
	latex->Draw();

	TLegend *legend = new TLegend(0.12,0.01,0.4,0.35);
	legend->AddEntry(gr11,"STAR Au-Au 0 - 5%","p");
	legend->AddEntry(gr13,"STAR Au-Au 70 - 80%","p");
	legend->AddEntry(gr15,"HRG + decay + cuts","l");
	legend->SetTextSize(textSize/p2->GetHNDC());
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	legend->Draw();

	p3->cd();

	auto *gr31 = new TGraphErrors(n_expData,snn_exp,ksigma2_0to5,(double*)0,ksigma2_0to5_staterr);
	gr31->SetMarkerSize(markerSize);
	gr31->SetMarkerColorAlpha(kOrange-3,transparency);
	gr31->SetLineColorAlpha(kOrange-3,transparency);
	gr31->SetLineWidth(lineWidth);
	gr31->SetMarkerStyle(kFullCircle);
	auto *gr32 = new TGraphErrors(n_expData,snn_exp,ksigma2_0to5,(double*)0,ksigma2_0to5_syserr);
	gr32->SetLineWidth(lineWidth);
	gr32->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr33 = new TGraphErrors(n_expData,snn_exp,ksigma2_70to80,(double*)0,ksigma2_70to80_staterr);
	gr33->SetMarkerSize(markerSize);
	gr33->SetMarkerColorAlpha(kViolet-1,transparency);
	gr33->SetLineColorAlpha(kViolet-1,transparency);
	gr33->SetLineWidth(lineWidth);
	gr33->SetMarkerStyle(kFullSquare);
	auto *gr34 = new TGraphErrors(n_expData,snn_exp,ksigma2_70to80,(double*)0,ksigma2_70to80_syserr);
	gr34->SetLineWidth(lineWidth);
	gr34->SetLineColorAlpha(kViolet-1,transparency);

	auto *gr35 = new TGraph(intvl,snn,sigma2byM_p_hrg);
	gr35->SetLineStyle(1);
	gr35->SetLineColorAlpha(kGreen+2,transparency);
	gr35->SetLineWidth(lineWidth);

	TMultiGraph  *mg3  = new TMultiGraph();

	mg3->Add(gr31,"pz");
	mg3->Add(gr32,"||");
	mg3->Add(gr33,"pz");
	mg3->Add(gr34,"||");
	mg3->Add(gr35,"l");

	mg3->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg3->GetYaxis()->SetTitle("k#sigma^{2}");
	mg3->GetYaxis()->SetRangeUser(0.25,1.15);
	mg3->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg3->GetXaxis()->CenterTitle(true);
	mg3->GetYaxis()->CenterTitle(true);
	mg3->GetXaxis()->SetNoExponent();
	mg3->GetXaxis()->SetLabelSize(textSize/p3->GetHNDC());
	mg3->GetXaxis()->SetTitleSize(textSize/p3->GetHNDC());
	mg3->GetXaxis()->SetTitleOffset(0.75);
	mg3->GetYaxis()->SetLabelSize(textSize/p3->GetHNDC());
	mg3->GetYaxis()->SetTitleSize(textSize/p3->GetHNDC());
	mg3->GetYaxis()->SetTickLength(0.02);
	mg3->GetYaxis()->SetTitleOffset(0.65);
	mg3->Draw("a");

	TText *text2 = new TText(0.2,0.9,"STAR Collab., PRL 112, 032302 (2014)");
	text2->SetNDC();
	text2->SetTextSize(textSize/p3->GetHNDC());
	text2->Draw();

	delete[] sigma2byM_0to5; delete[] sigma2byM_0to5_staterr; delete[] sigma2byM_0to5_syserr;
	delete[] Ssigma_0to5; delete[] Ssigma_0to5_syserr; delete[] Ssigma_0to5_staterr;
	delete[] ksigma2_0to5; delete[] ksigma2_0to5_syserr; delete[] ksigma2_0to5_staterr;

	delete[] sigma2byM_70to80; delete[] sigma2byM_70to80_staterr; delete[] sigma2byM_70to80_syserr;
	delete[] Ssigma_70to80; delete[] Ssigma_70to80_syserr; delete[] Ssigma_70to80_staterr;
	delete[] ksigma2_70to80; delete[] ksigma2_70to80_syserr; delete[] ksigma2_70to80_staterr;
}

void plotSimilarToSTARnetCharge(){

	double pT_min = 0.2;		//kinematic cuts
	double pT_max = 2.0;
	double eta_min = -0.5;		// pseudorapidity
	double eta_max = 0.5;

	double E_l = 5;
	double E_u = 250;
	int intvl = 100;
	double dE = pow(E_u/E_l,1.0/(intvl-1));

	double x,y;
	double num, denom;
	double T, mub, muq, mus;

	ParticleSystem system {};
	system.loadDefaultData();

	Susceptibility1Q chi1Q {system};
	Susceptibility2Q chi2Q {system};

	chi1Q.setDecay(true);
	chi2Q.setDecay(true);

	chi1Q.setDetectorPseudorapidityCoordinates();
	chi2Q.setDetectorPseudorapidityCoordinates();

	chi1Q.setIntegralLimits(pT_min, pT_max, eta_min, eta_max);
	chi2Q.setIntegralLimits(pT_min, pT_max, eta_min, eta_max);

	Susceptibility1 chi1 {system};
	Susceptibility2 chi2 {system};

	chi1.setDecay(true);
	chi2.setDecay(true);

	chi1.setDetectorPseudorapidityCoordinates();
	chi2.setDetectorPseudorapidityCoordinates();

	chi1.setIntegralLimits(pT_min, pT_max, eta_min, eta_max);
	chi2.setIntegralLimits(pT_min, pT_max, eta_min, eta_max);

	double snn[intvl];
	double sigma2byM_p_hrg[intvl];

	int n_expData;

	double *snn_exp;

	double *sigma2byM_0to5;
	double *sigma2byM_0to5_staterr;
	double *sigma2byM_0to5_syserr;

	double *Ssigma_0to5;
	double *Ssigma_0to5_syserr;
	double *Ssigma_0to5_staterr;

	double *ksigma2_0to5;
	double *ksigma2_0to5_syserr;
	double *ksigma2_0to5_staterr;

	double *sigma2byM_70to80;
	double *sigma2byM_70to80_staterr;
	double *sigma2byM_70to80_syserr;

	double *Ssigma_70to80;
	double *Ssigma_70to80_syserr;
	double *Ssigma_70to80_staterr;

	double *ksigma2_70to80;
	double *ksigma2_70to80_syserr;
	double *ksigma2_70to80_staterr;

	int species = system.getSpecies();
	SingleParticle *particle;

	for (int i = 0; i < intvl; ++i)
	{
		x = E_l*pow(dE,i);
		snn[i] = x;

		mub = 1.308/(1+0.273*x);
		muq = -0.0202/(1+0.126*x);
		mus = 0.214/(1+0.184*x);

		T = 0.166-0.139*mub*mub-0.053*pow(mub,4);
		
		num = 0; denom = 0;

		for (int i = 0; i < species; ++i)
		{
			particle = &system.getParticle(i);
			if(particle->isStable)
			{
				num += pow(particle->charge,2)*chi2.getValueEach(i,T,mub,muq,mus);
				denom += particle->charge*chi1.getValueEach(i,T,mub,muq,mus);
			}
		}

		// num = chi2Q.getValue(T,mub,muq,mus);
		// denom = chi1Q.getValue(T,mub,muq,mus);
		y = num/denom;
		sigma2byM_p_hrg[i] = y;
	}

	ifstream inFile;
	inFile.open("data/momentRatios_netQ.dat");
	if (inFile)
	{
		int l;
		inFile >> l;
		n_expData = l;

		snn_exp = new double[l];

		sigma2byM_0to5 = new double[l];
		sigma2byM_0to5_syserr = new double[l];
		sigma2byM_0to5_staterr = new double[l];

		Ssigma_0to5 = new double[l];
		Ssigma_0to5_syserr = new double[l];
		Ssigma_0to5_staterr = new double[l];

		ksigma2_0to5 = new double[l];
		ksigma2_0to5_syserr = new double[l];
		ksigma2_0to5_staterr = new double[l];

		sigma2byM_70to80 = new double[l];
		sigma2byM_70to80_syserr = new double[l];
		sigma2byM_70to80_staterr = new double[l];

		Ssigma_70to80 = new double[l];
		Ssigma_70to80_syserr = new double[l];
		Ssigma_70to80_staterr = new double[l];

		ksigma2_70to80 = new double[l];
		ksigma2_70to80_syserr = new double[l];
		ksigma2_70to80_staterr = new double[l];

		for (int i = 0; i < l; ++i)
		{
			inFile >> snn_exp[i] >> sigma2byM_0to5[i] >> sigma2byM_0to5_staterr[i] >> sigma2byM_0to5_syserr[i]
			 >> Ssigma_0to5[i] >> Ssigma_0to5_staterr[i] >> Ssigma_0to5_syserr[i]
			 >> ksigma2_0to5[i] >> ksigma2_0to5_staterr[i] >> ksigma2_0to5_syserr[i]
			 >> sigma2byM_70to80[i] >> sigma2byM_70to80_staterr[i] >> sigma2byM_70to80_syserr[i]
			 >> Ssigma_70to80[i] >> Ssigma_70to80_staterr[i] >> Ssigma_70to80_syserr[i]
			 >> ksigma2_70to80[i] >> ksigma2_70to80_staterr[i] >> ksigma2_70to80_syserr[i];
		}
	}
	else{
		cout << "Unable to load net-charge moment ratios...\n";
		n_expData = 0;
	}
	inFile.close();

	double topMargin = 0.01;
	double bottomMargin = 0.05;
	double rightMargin = 0.02;
	double leftMargin = 0.1;

	double textSize = 0.02;
	double markerSize = 2;
	double lineWidth = 3;
	double transparency = 0.75;

	TCanvas *cl1 = new TCanvas("cl_netQ","netQ",400,700);

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetEndErrorSize(10);

	double ySizeEach = (1 - topMargin - bottomMargin)/3;

	TPad *p1 = new TPad("p1", "p1",0,2*ySizeEach+bottomMargin,1,1);
	p1->SetLogx();
	p1->SetTickx();
	p1->SetBottomMargin(0);
	p1->SetTopMargin(topMargin/p1->GetHNDC());
	p1->SetRightMargin(rightMargin);
	p1->Draw();

	TPad *p2 = new TPad("p2", "p2",0,ySizeEach+bottomMargin,1,2*ySizeEach+bottomMargin);
	p2->SetLogx();
	p2->SetTickx();
	p2->SetTopMargin(0);
	p2->SetBottomMargin(0);
	p2->SetRightMargin(rightMargin);
	p2->Draw();

	TPad *p3 = new TPad("p3", "p3",0,0,1,ySizeEach+bottomMargin);
	p3->SetLogx();
	p3->SetTickx();
	p3->SetTopMargin(0);
	p3->SetRightMargin(rightMargin);
	p3->SetBottomMargin(bottomMargin/p3->GetHNDC());
	p3->Draw();

	p1->cd();

	auto *gr11 = new TGraphErrors(n_expData,snn_exp,sigma2byM_0to5,(double*)0,sigma2byM_0to5_staterr);
	gr11->SetMarkerSize(markerSize);
	gr11->SetMarkerColorAlpha(kOrange-3,transparency);
	gr11->SetLineColorAlpha(kOrange-3,transparency);
	gr11->SetLineWidth(lineWidth);
	gr11->SetMarkerStyle(kFullCircle);
	auto *gr12 = new TGraphErrors(n_expData,snn_exp,sigma2byM_0to5,(double*)0,sigma2byM_0to5_syserr);
	gr12->SetLineWidth(lineWidth);
	gr12->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr13 = new TGraphErrors(n_expData,snn_exp,sigma2byM_70to80,(double*)0,sigma2byM_70to80_staterr);
	gr13->SetMarkerSize(markerSize);
	gr13->SetMarkerColorAlpha(kViolet-1,transparency);
	gr13->SetLineColorAlpha(kViolet-1,transparency);
	gr13->SetLineWidth(lineWidth);
	gr13->SetMarkerStyle(kFullSquare);
	auto *gr14 = new TGraphErrors(n_expData,snn_exp,sigma2byM_70to80,(double*)0,sigma2byM_70to80_syserr);
	gr14->SetLineWidth(lineWidth);
	gr14->SetLineColorAlpha(kViolet-1,transparency);

	auto *gr15 = new TGraph(intvl,snn,sigma2byM_p_hrg);
	gr15->SetLineStyle(1);
	gr15->SetLineColorAlpha(kGreen+2,transparency);
	gr15->SetLineWidth(lineWidth);

	TMultiGraph  *mg1  = new TMultiGraph();

	mg1->Add(gr11,"pz");
	mg1->Add(gr12,"||");
	mg1->Add(gr13,"pz");
	mg1->Add(gr14,"||");
	mg1->Add(gr15,"l");

	// mg1->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg1->GetYaxis()->SetTitle("#sigma^{2}/M");
	// mg1->GetXaxis()->SetRangeUser(2,250);
	mg1->GetYaxis()->SetRangeUser(-7.5,112.5);
	mg1->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg1->GetXaxis()->CenterTitle(true);
	mg1->GetYaxis()->CenterTitle(true);
	mg1->GetXaxis()->SetNoExponent();
	mg1->GetXaxis()->SetLabelSize(textSize/p1->GetHNDC());
	mg1->GetXaxis()->SetTitleSize(textSize/p1->GetHNDC());
	mg1->GetXaxis()->SetTitleOffset(0.75);
	mg1->GetYaxis()->SetLabelSize(textSize/p1->GetHNDC());
	mg1->GetYaxis()->SetTitleSize(textSize/p1->GetHNDC());
	mg1->GetYaxis()->SetTickLength(0.02);
	mg1->GetYaxis()->SetTitleOffset(0.65);
	mg1->Draw("a");

	TLatex *text1 = new TLatex(0.4,0.75,"#splitline{net-charge results}{and calculations}");
	text1->SetNDC();
	text1->SetTextSize(textSize/p1->GetHNDC());
	text1->Draw();

	p2->cd();

	auto *gr21 = new TGraphErrors(n_expData,snn_exp,Ssigma_0to5,(double*)0,Ssigma_0to5_staterr);
	gr21->SetMarkerSize(markerSize);
	gr21->SetMarkerColorAlpha(kOrange-3,transparency);
	gr21->SetLineColorAlpha(kOrange-3,transparency);
	gr21->SetLineWidth(lineWidth);
	gr21->SetMarkerStyle(kFullCircle);
	auto *gr22 = new TGraphErrors(n_expData,snn_exp,Ssigma_0to5,(double*)0,Ssigma_0to5_syserr);
	gr22->SetLineWidth(lineWidth);
	gr22->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr23 = new TGraphErrors(n_expData,snn_exp,Ssigma_70to80,(double*)0,Ssigma_70to80_staterr);
	gr23->SetMarkerSize(markerSize);
	gr23->SetMarkerColorAlpha(kViolet-1,transparency);
	gr23->SetLineColorAlpha(kViolet-1,transparency);
	gr23->SetLineWidth(lineWidth);
	gr23->SetMarkerStyle(kFullSquare);
	auto *gr24 = new TGraphErrors(n_expData,snn_exp,Ssigma_70to80,(double*)0,Ssigma_70to80_syserr);
	gr24->SetLineWidth(lineWidth);
	gr24->SetLineColorAlpha(kViolet-1,transparency);

	auto *gr25 = new TGraph(intvl,snn,sigma2byM_p_hrg);
	gr25->SetLineStyle(1);
	gr25->SetLineColorAlpha(kGreen+2,transparency);
	gr25->SetLineWidth(lineWidth);

	TMultiGraph  *mg2  = new TMultiGraph();

	mg2->Add(gr21,"pz");
	mg2->Add(gr22,"||");
	mg2->Add(gr23,"pz");
	mg2->Add(gr24,"||");
	mg2->Add(gr25,"l");

	// mg2->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg2->GetYaxis()->SetTitle("S#sigma");
	// mg2->GetXaxis()->SetLimits(5,250);
	mg2->GetYaxis()->SetRangeUser(-0.095,0.75);
	mg2->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg2->GetXaxis()->CenterTitle(true);
	mg2->GetYaxis()->CenterTitle(true);
	mg2->GetXaxis()->SetNoExponent();
	mg2->GetXaxis()->SetLabelSize(textSize/p2->GetHNDC());
	mg2->GetXaxis()->SetTitleSize(textSize/p2->GetHNDC());
	mg2->GetXaxis()->SetTitleOffset(0.75);
	mg2->GetYaxis()->SetLabelSize(textSize/p2->GetHNDC());
	mg2->GetYaxis()->SetTitleSize(textSize/p2->GetHNDC());
	mg2->GetYaxis()->SetTickLength(0.02);
	mg2->GetYaxis()->SetTitleOffset(0.65);
	mg2->Draw("a");

	TLegend *legend = new TLegend(0.5,0.6,0.85,0.85);
	legend->AddEntry(gr11,"STAR Au-Au 0 - 5%","p");
	legend->AddEntry(gr13,"STAR Au-Au 70 - 80%","p");
	legend->AddEntry(gr15,"HRG + decay + cuts","l");
	legend->SetTextSize(textSize/p2->GetHNDC());
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	legend->Draw();

	p3->cd();

	auto *gr31 = new TGraphErrors(n_expData,snn_exp,ksigma2_0to5,(double*)0,ksigma2_0to5_staterr);
	gr31->SetMarkerSize(markerSize);
	gr31->SetMarkerColorAlpha(kOrange-3,transparency);
	gr31->SetLineColorAlpha(kOrange-3,transparency);
	gr31->SetLineWidth(lineWidth);
	gr31->SetMarkerStyle(kFullCircle);
	auto *gr32 = new TGraphErrors(n_expData,snn_exp,ksigma2_0to5,(double*)0,ksigma2_0to5_syserr);
	gr32->SetLineWidth(lineWidth);
	gr32->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr33 = new TGraphErrors(n_expData,snn_exp,ksigma2_70to80,(double*)0,ksigma2_70to80_staterr);
	gr33->SetMarkerSize(markerSize);
	gr33->SetMarkerColorAlpha(kViolet-1,transparency);
	gr33->SetLineColorAlpha(kViolet-1,transparency);
	gr33->SetLineWidth(lineWidth);
	gr33->SetMarkerStyle(kFullSquare);
	auto *gr34 = new TGraphErrors(n_expData,snn_exp,ksigma2_70to80,(double*)0,ksigma2_70to80_syserr);
	gr34->SetLineWidth(lineWidth);
	gr34->SetLineColorAlpha(kViolet-1,transparency);

	auto *gr35 = new TGraph(intvl,snn,sigma2byM_p_hrg);
	gr35->SetLineStyle(1);
	gr35->SetLineColorAlpha(kGreen+2,transparency);
	gr35->SetLineWidth(lineWidth);

	TMultiGraph  *mg3  = new TMultiGraph();

	mg3->Add(gr31,"pz");
	mg3->Add(gr32,"||");
	mg3->Add(gr33,"pz");
	mg3->Add(gr34,"||");
	mg3->Add(gr35,"l");

	mg3->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg3->GetYaxis()->SetTitle("k#sigma^{2}");
	mg3->GetYaxis()->SetRangeUser(-19.5,12.5);
	mg3->GetYaxis()->SetNdivisions(6, 2, 0, kTRUE);
	mg3->GetXaxis()->CenterTitle(true);
	mg3->GetYaxis()->CenterTitle(true);
	mg3->GetXaxis()->SetNoExponent();
	mg3->GetXaxis()->SetLabelSize(textSize/p3->GetHNDC());
	mg3->GetXaxis()->SetTitleSize(textSize/p3->GetHNDC());
	mg3->GetXaxis()->SetTitleOffset(0.75);
	mg3->GetYaxis()->SetLabelSize(textSize/p3->GetHNDC());
	mg3->GetYaxis()->SetTitleSize(textSize/p3->GetHNDC());
	mg3->GetYaxis()->SetTickLength(0.02);
	mg3->GetYaxis()->SetTitleOffset(0.65);
	mg3->Draw("a");

	TText *text2 = new TText(0.32,0.9,"STAR Collab., PRL 113, 092301 (2014)");
	text2->SetNDC();
	text2->SetTextSize(textSize/p3->GetHNDC());
	text2->Draw();

	TLatex *latex = new TLatex(0.5,0.35,"#splitline{0.2 < p_{T} < 2.0 (GeV/c)}{#left|#eta#right| #leq 0.5}");
	latex->SetNDC();
	latex->SetTextSize(textSize/p2->GetHNDC());
	latex->Draw();

	delete[] sigma2byM_0to5; delete[] sigma2byM_0to5_staterr; delete[] sigma2byM_0to5_syserr;
	delete[] Ssigma_0to5; delete[] Ssigma_0to5_syserr; delete[] Ssigma_0to5_staterr;
	delete[] ksigma2_0to5; delete[] ksigma2_0to5_syserr; delete[] ksigma2_0to5_staterr;

	delete[] sigma2byM_70to80; delete[] sigma2byM_70to80_staterr; delete[] sigma2byM_70to80_syserr;
	delete[] Ssigma_70to80; delete[] Ssigma_70to80_syserr; delete[] Ssigma_70to80_staterr;
	delete[] ksigma2_70to80; delete[] ksigma2_70to80_syserr; delete[] ksigma2_70to80_staterr;
}

int main()
{
	// numericalMinimization();
	TApplication app {"app", nullptr, nullptr};
	numericalMinimizationCumulants();

	cout << "\nPress Ctrl-C to exit..." << endl;

	// generateCumulantRatios_netP();
	// plotSimilarToSTARnetProton();
	// plotSimilarToSTARnetCharge();

	app.Run();
	getchar();
	return 0;
}