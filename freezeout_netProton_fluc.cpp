#include <iostream>
#include <cmath>
#include <fstream>

#include "Particles.h"
#include "ThermalFunctions.h"
#include "ThermalMinimizer.h"

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

class STAR_netProton_Ssigma:public ThermalFunction
{

	Susceptibility3B *chi3B;
	Susceptibility2B *chi2B;

	public:

	STAR_netProton_Ssigma()
	{
		ParticleSystem system;
		system.loadDefaultData();

		chi3B = new Susceptibility3B {system};		
		chi3B->setDetectorRapidityCoordinates();
		chi3B->setIntegralLimits(0.4,0.8,-0.5,0.5);
		chi3B->setDecay(true);

		chi2B = new Susceptibility2B {system};
		chi2B->setDecay(true);
		chi2B->setIntegralLimits(0.4,0.8,-0.5,0.5);
		chi2B->setDetectorRapidityCoordinates();
	}
	~STAR_netProton_Ssigma()
	{
		delete chi3B;
		delete chi2B;
		chi3B = nullptr;
		chi2B = nullptr;
	}

	double getValue(double T, double mub, double muq, double mus)
	{
		double denom = chi2B->getValueEach(23,T,mub,muq,mus)+chi2B->getValueEach(24,T,mub,muq,mus);
		double num = chi3B->getValueEach(23,T,mub,muq,mus)+chi3B->getValueEach(24,T,mub,muq,mus);
		return num/denom;
	}
};

class STAR_netProton_ksigma2:public ThermalFunction
{

	Susceptibility4B *chi4B;
	Susceptibility2B *chi2B;

	public:

	STAR_netProton_ksigma2()
	{
		ParticleSystem system;
		system.loadDefaultData();

		chi4B = new Susceptibility4B {system};		
		chi4B->setDetectorRapidityCoordinates();
		chi4B->setIntegralLimits(0.4,0.8,-0.5,0.5);
		chi4B->setDecay(true);

		chi2B = new Susceptibility2B {system};
		chi2B->setDecay(true);
		chi2B->setIntegralLimits(0.4,0.8,-0.5,0.5);
		chi2B->setDetectorRapidityCoordinates();
	}
	~STAR_netProton_ksigma2()
	{
		delete chi4B;
		delete chi2B;
		chi4B = nullptr;
		chi2B = nullptr;
	}

	double getValue(double T, double mub, double muq, double mus)
	{
		double denom = chi2B->getValueEach(23,T,mub,muq,mus)+chi2B->getValueEach(24,T,mub,muq,mus);
		double num = chi4B->getValueEach(23,T,mub,muq,mus)+chi4B->getValueEach(24,T,mub,muq,mus);
		return num/denom;
	}
};

class STAR_netLambda_Sigma2byM:public ThermalFunction
{

	Susceptibility1S *chi1S;
	Susceptibility2S *chi2S;

	public:

	STAR_netLambda_Sigma2byM()
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
	~STAR_netLambda_Sigma2byM()
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

void findChemFreezeoutParam_netProton_fluc()
{
	int n_expData;

	double *energy_exp, *netProton_energy_exp, *netLambda_energy_exp;


	double *netProton_sigma2byM_exp_0to5, *netProton_sigma2byM_exp_0to5_staterr, *netProton_sigma2byM_exp_0to5_syserr;
	double *netProton_Ssigma_exp_0to5, *netProton_Ssigma_exp_0to5_syserr, *netProton_Ssigma_exp_0to5_staterr;
	double *netProton_ksigma2_exp_0to5, *netProton_ksigma2_exp_0to5_syserr, *netProton_ksigma2_exp_0to5_staterr;

	double *netProton_sigma2byM_exp_70to80, *netProton_sigma2byM_exp_70to80_staterr, *netProton_sigma2byM_exp_70to80_syserr;
	double *netProton_Ssigma_exp_70to80, *netProton_Ssigma_exp_70to80_syserr, *netProton_Ssigma_exp_70to80_staterr;
	double *netProton_ksigma2_exp_70to80, *netProton_ksigma2_exp_70to80_syserr, *netProton_ksigma2_exp_70to80_staterr;


	double *netLambda_sigma2byM_exp_0to5, *netLambda_sigma2byM_exp_0to5_staterr, *netLambda_sigma2byM_exp_0to5_syserr;
	double *netLambda_Ssigma_exp_0to5, *netLambda_Ssigma_exp_0to5_syserr, *netLambda_Ssigma_exp_0to5_staterr;

	double *netLambda_sigma2byM_exp_70to80, *netLambda_sigma2byM_exp_70to80_staterr, *netLambda_sigma2byM_exp_70to80_syserr;
	double *netLambda_Ssigma_exp_70to80, *netLambda_Ssigma_exp_70to80_syserr, *netLambda_Ssigma_exp_70to80_staterr;

	int l, l1, l2;

	ifstream inFilenetProton;
	ifstream inFileNetLambda;
	inFilenetProton.open("data/cumulantRatios_netProton_freezeout_netProtonLambda_fluc.dat");
	inFileNetLambda.open("data/cumulantRatios_netLambda_freezeout_netProtonLambda_fluc.dat");
	if (inFilenetProton && inFileNetLambda)
	{
		inFilenetProton >> l1;
		inFileNetLambda >> l2;

		if (l1 != l2)
		{
			cout << "Data size mismatch...\n";
		}

		n_expData = min(l1,l2);
		l = n_expData;

		energy_exp = new double[l]; netProton_energy_exp = new double[l]; netLambda_energy_exp = new double[l]; 


		netProton_sigma2byM_exp_0to5 = new double[l]; netProton_sigma2byM_exp_0to5_syserr = new double[l]; netProton_sigma2byM_exp_0to5_staterr = new double[l];
		netProton_Ssigma_exp_0to5 = new double[l]; netProton_Ssigma_exp_0to5_syserr = new double[l]; netProton_Ssigma_exp_0to5_staterr = new double[l];
		netProton_ksigma2_exp_0to5 = new double[l]; netProton_ksigma2_exp_0to5_syserr = new double[l]; netProton_ksigma2_exp_0to5_staterr = new double[l];

		netProton_sigma2byM_exp_70to80 = new double[l]; netProton_sigma2byM_exp_70to80_syserr = new double[l]; netProton_sigma2byM_exp_70to80_staterr = new double[l];
		netProton_Ssigma_exp_70to80 = new double[l]; netProton_Ssigma_exp_70to80_syserr = new double[l]; netProton_Ssigma_exp_70to80_staterr = new double[l];
		netProton_ksigma2_exp_70to80 = new double[l]; netProton_ksigma2_exp_70to80_syserr = new double[l]; netProton_ksigma2_exp_70to80_staterr = new double[l];


		netLambda_sigma2byM_exp_0to5 = new double[l]; netLambda_sigma2byM_exp_0to5_syserr = new double[l]; netLambda_sigma2byM_exp_0to5_staterr = new double[l];
		netLambda_Ssigma_exp_0to5 = new double[l]; netLambda_Ssigma_exp_0to5_syserr = new double[l]; netLambda_Ssigma_exp_0to5_staterr = new double[l];

		netLambda_sigma2byM_exp_70to80 = new double[l]; netLambda_sigma2byM_exp_70to80_syserr = new double[l]; netLambda_sigma2byM_exp_70to80_staterr = new double[l];
		netLambda_Ssigma_exp_70to80 = new double[l]; netLambda_Ssigma_exp_70to80_syserr = new double[l]; netLambda_Ssigma_exp_70to80_staterr = new double[l];

		for (int i = 0; i < l; ++i)
		{
			inFilenetProton >> netProton_energy_exp[i] >> netProton_sigma2byM_exp_0to5[i] >> netProton_sigma2byM_exp_0to5_staterr[i] >> netProton_sigma2byM_exp_0to5_syserr[i]
			 >> netProton_Ssigma_exp_0to5[i] >> netProton_Ssigma_exp_0to5_staterr[i] >> netProton_Ssigma_exp_0to5_syserr[i]
			 >> netProton_ksigma2_exp_0to5[i] >> netProton_ksigma2_exp_0to5_staterr[i] >> netProton_ksigma2_exp_0to5_syserr[i]
			 >> netProton_sigma2byM_exp_70to80[i] >> netProton_sigma2byM_exp_70to80_staterr[i] >> netProton_sigma2byM_exp_70to80_syserr[i]
			 >> netProton_Ssigma_exp_70to80[i] >> netProton_Ssigma_exp_70to80_staterr[i] >> netProton_Ssigma_exp_70to80_syserr[i]
			 >> netProton_ksigma2_exp_70to80[i] >> netProton_ksigma2_exp_70to80_staterr[i] >> netProton_ksigma2_exp_70to80_syserr[i];

			inFileNetLambda >> netLambda_energy_exp[i] >> netLambda_sigma2byM_exp_0to5[i] >> netLambda_sigma2byM_exp_0to5_staterr[i] >> netLambda_sigma2byM_exp_0to5_syserr[i]
			 >> netLambda_Ssigma_exp_0to5[i] >> netLambda_Ssigma_exp_0to5_staterr[i] >> netLambda_Ssigma_exp_0to5_syserr[i]
			 >> netLambda_sigma2byM_exp_70to80[i] >> netLambda_sigma2byM_exp_70to80_staterr[i] >> netLambda_sigma2byM_exp_70to80_syserr[i]
			 >> netLambda_Ssigma_exp_70to80[i] >> netLambda_Ssigma_exp_70to80_staterr[i] >> netLambda_Ssigma_exp_70to80_syserr[i];

			 if (netProton_energy_exp[i] != netLambda_energy_exp[i])
			 {
			 	cout << "Energy value mismatch at line " << i+2 << '\n';
			 }
			 energy_exp[i] = (netProton_energy_exp[i]+netLambda_energy_exp[i])/2.0;
		}
	}
	else{
		cout << "Unable to load cumulant/moment ratios...\n";
		n_expData = 0;
	}
	inFilenetProton.close();
	inFileNetLambda.close();

	double E_l = 5;
	double E_u = 250;
	int intvl = 100;
	double dE = pow(E_u/E_l,1.0/(intvl-1));

	double x;
	double T, mub, muq, mus;

	double energy[intvl];

	double netProton_sigma2byM_hrg[intvl], netProton_Ssigma_hrg[intvl], netProton_ksigma2_hrg[intvl];
	double netLambda_sigma2byM_hrg[intvl], netLambda_Ssigma_hrg[intvl];

	STAR_netProton_sigma2byM netProton_sig2byM;
	STAR_netProton_Ssigma netProton_Ssigma;
	STAR_netProton_ksigma2 netProton_ksigma2;

	STAR_netLambda_Sigma2byM netLambda_sig2byM;
	STAR_netLambda_Ssigma netLambda_Ssigma;

	for (int i = 0; i < intvl; ++i)
	{
		x = E_l*pow(dE,i);
		energy[i] = x;

		mub = 1.308/(1+0.273*x);
		muq = -0.0202/(1+0.126*x);
		mus = 0.214/(1+0.184*x);

		T = 0.166-0.139*mub*mub-0.053*pow(mub,4);

		netProton_sigma2byM_hrg[i] = netProton_sig2byM.getValue(T,mub,muq,mus);
		netProton_Ssigma_hrg[i] = netProton_Ssigma.getValue(T,mub,muq,mus);
		netProton_ksigma2_hrg[i] = netProton_ksigma2.getValue(T,mub,muq,mus);

		netLambda_sigma2byM_hrg[i] = netLambda_sig2byM.getValue(T,mub,muq,mus);
		netLambda_Ssigma_hrg[i] = netLambda_Ssigma.getValue(T,mub,muq,mus);
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

	double netProton_sigma2byM_fit[n_expData], netProton_Ssigma_fit[n_expData], netProton_ksigma2_fit[n_expData];
	double netLambda_sigma2byM_fit[n_expData], netLambda_Ssigma_fit[n_expData];

	cout << "-----------------------------------------------------------------" << '\n';
	cout << "        Determination of chemical freezout parameters using      " << '\n';
	cout << "                      net-proton fluctuation                     " << '\n';
	cout << "-----------------------------------------------------------------" << '\n';
	ThermalMinimizer minimum {3};
	minimum.SetThermalFunctionEach(0,&netProton_sig2byM);
	// minimum.SetThermalFunctionEach(1,&netLambda_sig2byM);
	minimum.SetThermalFunctionEach(1,&nS);
	minimum.SetThermalFunctionEach(2,&nQ_by_nB);

	ofstream outFileFreeze;
	outFileFreeze.open("data/freezeout_netProton_fluc.dat");

	bool writeToFile = false;
	if (outFileFreeze)
	{
	 outFileFreeze << n_expData << '\t' << 9 << '\n';
	 writeToFile = true;
	}

	double err;

	for (int i = 0; i < n_expData; ++i)
	{
		x = energy_exp[i];
		cout << "                    sqrt(snn) = " << x << " GeV" << '\n';

		mub = 1.308/(1+0.273*x);
		muq = -0.0202/(1+0.126*x);
		mus = 0.214/(1+0.184*x);

		T = 0.166-0.139*mub*mub-0.053*pow(mub,4);

		err = sqrt(pow(netProton_sigma2byM_exp_0to5_staterr[i],2)+pow(netProton_sigma2byM_exp_0to5_syserr[i],2));
		minimum.SetFunctionValueEach(0,netProton_sigma2byM_exp_0to5[i],err);
		err = sqrt(pow(netLambda_sigma2byM_exp_0to5_staterr[i],2)+pow(netLambda_sigma2byM_exp_0to5_syserr[i],2));
		// minimum.SetFunctionValueEach(1,netLambda_sigma2byM_exp_0to5[i],err);
		minimum.SetFunctionValueEach(1,0,1e-3);
		minimum.SetFunctionValueEach(2,0.4,1e-3);

		minimum.setInitialParameters(T,mub,muq,mus);
		minimum.minimize(1000);

		minPar = minimum.getMinimizedParameters();
		minPar_err = minimum.getMinimizedParameterErrors();

		T_ch[i] = minPar[0]; mub_ch[i] = minPar[1]; muq_ch[i] = minPar[2]; mus_ch[i] = minPar[3];
		T_ch_err[i] = minPar_err[0]; mub_ch_err[i] = minPar_err[1]; muq_ch_err[i] = minPar_err[2]; mus_ch_err[i] = minPar_err[3];

		if (writeToFile)
		{
			outFileFreeze << x << '\t' << minPar[0] << '\t' << minPar_err[0] << '\t' << minPar[1] << '\t' << minPar_err[1]
			 << '\t' << minPar[2] << '\t' << minPar_err[2] << '\t' << minPar[3] << '\t' << minPar_err[3] << '\n';
		}
		
		netProton_sigma2byM_fit[i] = netProton_sig2byM.getValue(minPar[0],minPar[1],minPar[2],minPar[3]);
		netProton_Ssigma_fit[i] = netProton_Ssigma.getValue(minPar[0],minPar[1],minPar[2],minPar[3]);
		netProton_ksigma2_fit[i] = netProton_ksigma2.getValue(minPar[0],minPar[1],minPar[2],minPar[3]);

		netLambda_sigma2byM_fit[i] = netLambda_sig2byM.getValue(minPar[0],minPar[1],minPar[2],minPar[3]);
		netLambda_Ssigma_fit[i] = netLambda_Ssigma.getValue(minPar[0],minPar[1],minPar[2],minPar[3]);
	}

	if (writeToFile)
	{
		outFileFreeze << "\n\\\\ energy\tT_ch\tT_ch_err\tmub_ch\tmub_ch_err\tmuq_ch\tmuq_ch_err\tmus_ch\tmus_ch_err\n" ;
		cout << "Freeze out parameters saved to data/freezeout_netProton_fluc.dat..." << '\n';
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

	TCanvas *cl1 = new TCanvas("cl1","netProton_fluc",400,700);

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

	auto *gr111 = new TGraph(intvl,energy,netProton_sigma2byM_hrg);
	gr111->SetLineStyle(1);
	gr111->SetLineColorAlpha(kGreen+2,transparency);
	gr111->SetLineWidth(lineWidth);

	auto *gr112 = new TGraphErrors(n_expData,energy_exp,netProton_sigma2byM_exp_0to5,(double*)0,netProton_sigma2byM_exp_0to5_staterr);
	gr112->SetMarkerSize(markerSize);
	gr112->SetMarkerColorAlpha(kOrange-3,transparency);
	gr112->SetLineColorAlpha(kOrange-3,transparency);
	gr112->SetLineWidth(lineWidth);
	gr112->SetMarkerStyle(kFullCircle);
	auto *gr113 = new TGraphErrors(n_expData,energy_exp,netProton_sigma2byM_exp_0to5,(double*)0,netProton_sigma2byM_exp_0to5_syserr);
	gr113->SetLineWidth(lineWidth);
	gr113->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr114 = new TGraphErrors(n_expData,energy_exp,netProton_sigma2byM_exp_70to80,(double*)0,netProton_sigma2byM_exp_70to80_staterr);
	gr114->SetMarkerSize(markerSize);
	gr114->SetMarkerColorAlpha(kViolet-1,transparency);
	gr114->SetLineColorAlpha(kViolet-1,transparency);
	gr114->SetLineWidth(lineWidth);
	gr114->SetMarkerStyle(kFullSquare);
	auto *gr115 = new TGraphErrors(n_expData,energy_exp,netProton_sigma2byM_exp_70to80,(double*)0,netProton_sigma2byM_exp_70to80_syserr);
	gr115->SetLineWidth(lineWidth);
	gr115->SetLineColorAlpha(kViolet-1,transparency);
	auto *gr116 = new TGraph(n_expData,energy_exp,netProton_sigma2byM_fit);
	gr116->SetMarkerSize(markerSize);
	gr116->SetMarkerColorAlpha(kBlue,transparency);
	gr116->SetMarkerStyle(kFullTriangleUp);

	TMultiGraph *mg11  = new TMultiGraph();

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
	mg11->GetYaxis()->SetTitleOffset(0.75);
	mg11->Draw("a");

	TLatex *text11 = new TLatex(0.4,0.75,"#splitline{net-proton results}{and calculations}");
	text11->SetNDC();
	text11->SetTextSize(textSize/p11->GetHNDC());
	text11->Draw();

	p12->cd();

	auto *gr121 = new TGraphErrors(n_expData,energy_exp,netProton_Ssigma_exp_0to5,(double*)0,netProton_Ssigma_exp_0to5_staterr);
	gr121->SetMarkerSize(markerSize);
	gr121->SetMarkerColorAlpha(kOrange-3,transparency);
	gr121->SetLineColorAlpha(kOrange-3,transparency);
	gr121->SetLineWidth(lineWidth);
	gr121->SetMarkerStyle(kFullCircle);
	auto *gr122 = new TGraphErrors(n_expData,energy_exp,netProton_Ssigma_exp_0to5,(double*)0,netProton_Ssigma_exp_0to5_syserr);
	gr122->SetLineWidth(lineWidth);
	gr122->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr123 = new TGraphErrors(n_expData,energy_exp,netProton_Ssigma_exp_70to80,(double*)0,netProton_Ssigma_exp_70to80_staterr);
	gr123->SetMarkerSize(markerSize);
	gr123->SetMarkerColorAlpha(kViolet-1,transparency);
	gr123->SetLineColorAlpha(kViolet-1,transparency);
	gr123->SetLineWidth(lineWidth);
	gr123->SetMarkerStyle(kFullSquare);
	auto *gr124 = new TGraphErrors(n_expData,energy_exp,netProton_Ssigma_exp_70to80,(double*)0,netProton_Ssigma_exp_70to80_syserr);
	gr124->SetLineWidth(lineWidth);
	gr124->SetLineColorAlpha(kViolet-1,transparency);

	auto *gr125 = new TGraph(intvl,energy,netProton_Ssigma_hrg);
	gr125->SetLineStyle(1);
	gr125->SetLineColorAlpha(kGreen+2,transparency);
	gr125->SetLineWidth(lineWidth);
	auto *gr126 = new TGraph(n_expData,energy_exp,netProton_Ssigma_fit);
	gr126->SetMarkerSize(markerSize);
	gr126->SetMarkerColorAlpha(kBlue,transparency);
	gr126->SetMarkerStyle(kFullTriangleUp);

	TMultiGraph *mg12  = new TMultiGraph();

	mg12->Add(gr121,"pz");
	mg12->Add(gr122,"||");
	mg12->Add(gr123,"pz");
	mg12->Add(gr124,"||");
	mg12->Add(gr125,"l");
	mg12->Add(gr126,"p");

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
	mg12->GetYaxis()->SetTitleOffset(0.75);
	mg12->Draw("a");

	TLatex *latex1 = new TLatex(0.6,0.75,"#splitline{0.4 < p_{T} < 0.8 (GeV/c)}{#left|y#right| #leq 0.5}");
	latex1->SetNDC();
	latex1->SetTextSize(textSize/p12->GetHNDC());
	latex1->Draw();

	TLegend *legend1 = new TLegend(0.12,0.02,0.4,0.35);
	legend1->AddEntry(gr112,"STAR Au-Au 0 - 5%","p");
	legend1->AddEntry(gr114,"STAR Au-Au 70 - 80%","p");
	legend1->AddEntry(gr111,"HRG + decay + cuts","l");
	legend1->AddEntry(gr116,"HRG + decay + cuts (fit)","p");
	legend1->SetTextSize(textSize/p12->GetHNDC());
	legend1->SetBorderSize(0);
	legend1->SetFillStyle(0);
	legend1->Draw();

	p13->cd();

	auto *gr131 = new TGraphErrors(n_expData,energy_exp,netProton_ksigma2_exp_0to5,(double*)0,netProton_ksigma2_exp_0to5_staterr);
	gr131->SetMarkerSize(markerSize);
	gr131->SetMarkerColorAlpha(kOrange-3,transparency);
	gr131->SetLineColorAlpha(kOrange-3,transparency);
	gr131->SetLineWidth(lineWidth);
	gr131->SetMarkerStyle(kFullCircle);
	auto *gr132 = new TGraphErrors(n_expData,energy_exp,netProton_ksigma2_exp_0to5,(double*)0,netProton_ksigma2_exp_0to5_syserr);
	gr132->SetLineWidth(lineWidth);
	gr132->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr133 = new TGraphErrors(n_expData,energy_exp,netProton_ksigma2_exp_70to80,(double*)0,netProton_ksigma2_exp_70to80_staterr);
	gr133->SetMarkerSize(markerSize);
	gr133->SetMarkerColorAlpha(kViolet-1,transparency);
	gr133->SetLineColorAlpha(kViolet-1,transparency);
	gr133->SetLineWidth(lineWidth);
	gr133->SetMarkerStyle(kFullSquare);
	auto *gr134 = new TGraphErrors(n_expData,energy_exp,netProton_ksigma2_exp_70to80,(double*)0,netProton_ksigma2_exp_70to80_syserr);
	gr134->SetLineWidth(lineWidth);
	gr134->SetLineColorAlpha(kViolet-1,transparency);

	auto *gr135 = new TGraph(intvl,energy,netProton_ksigma2_hrg);
	gr135->SetLineStyle(1);
	gr135->SetLineColorAlpha(kGreen+2,transparency);
	gr135->SetLineWidth(lineWidth);
	auto *gr136 = new TGraph(n_expData,energy_exp,netProton_ksigma2_fit);
	gr136->SetMarkerSize(markerSize);
	gr136->SetMarkerColorAlpha(kBlue,transparency);
	gr136->SetMarkerStyle(kFullTriangleUp);

	TMultiGraph  *mg13  = new TMultiGraph();

	mg13->Add(gr131,"pz");
	mg13->Add(gr132,"||");
	mg13->Add(gr133,"pz");
	mg13->Add(gr134,"||");
	mg13->Add(gr135,"l");
	mg13->Add(gr136,"p");

	mg13->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg13->GetYaxis()->SetTitle("k#sigma^{2}");
	mg13->GetYaxis()->SetRangeUser(0.37,1.17);
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
	mg13->GetYaxis()->SetTitleOffset(0.75);
	mg13->Draw("a");

	TText *text12 = new TText(0.15,0.9,"STAR Collab., PRL 112, 032302 (2014)");
	text12->SetNDC();
	text12->SetTextSize(textSize/p13->GetHNDC());
	text12->Draw();

	cl1->SaveAs("plots/netProton_fluc.pdf");
	
	TCanvas *cl2 = new TCanvas("cl2","netLambda_fluc",400,600);

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetEndErrorSize(10);

	textSize = 0.025;
	
	ySizeEach = (1 - topMargin - bottomMargin)/2;

	TPad *p21 = new TPad("p21", "p21",0,ySizeEach+bottomMargin,1,1);
	p21->SetLogx();
	p21->SetTickx();
	p21->SetBottomMargin(0);
	p21->SetTopMargin(topMargin/p21->GetHNDC());
	p21->SetRightMargin(rightMargin);
	p21->Draw();

	TPad *p22 = new TPad("p22", "p22",0,0,1,ySizeEach+bottomMargin);
	p22->SetLogx();
	p22->SetTickx();
	p22->SetTopMargin(0);
	p22->SetBottomMargin(bottomMargin/p22->GetHNDC());
	p22->SetRightMargin(rightMargin);
	p22->Draw();

	p21->cd();

	auto *gr211 = new TGraph(intvl,energy,netLambda_sigma2byM_hrg);
	gr211->SetLineStyle(1);
	gr211->SetLineColorAlpha(kGreen+2,transparency);
	gr211->SetLineWidth(lineWidth);

	auto *gr212 = new TGraphErrors(n_expData,energy_exp,netLambda_sigma2byM_exp_0to5,(double*)0,netLambda_sigma2byM_exp_0to5_staterr);
	gr212->SetMarkerSize(markerSize);
	gr212->SetMarkerColorAlpha(kOrange-3,transparency);
	gr212->SetLineColorAlpha(kOrange-3,transparency);
	gr212->SetLineWidth(lineWidth);
	gr212->SetMarkerStyle(kFullCircle);
	auto *gr213 = new TGraphErrors(n_expData,energy_exp,netLambda_sigma2byM_exp_0to5,(double*)0,netLambda_sigma2byM_exp_0to5_syserr);
	gr213->SetLineWidth(lineWidth);
	gr213->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr214 = new TGraphErrors(n_expData,energy_exp,netLambda_sigma2byM_exp_70to80,(double*)0,netLambda_sigma2byM_exp_70to80_staterr);
	gr214->SetMarkerSize(markerSize);
	gr214->SetMarkerColorAlpha(kViolet-1,transparency);
	gr214->SetLineColorAlpha(kViolet-1,transparency);
	gr214->SetLineWidth(lineWidth);
	gr214->SetMarkerStyle(kFullSquare);
	auto *gr215 = new TGraphErrors(n_expData,energy_exp,netLambda_sigma2byM_exp_70to80,(double*)0,netLambda_sigma2byM_exp_70to80_syserr);
	gr215->SetLineWidth(lineWidth);
	gr215->SetLineColorAlpha(kViolet-1,transparency);
	auto *gr216 = new TGraph(n_expData,energy_exp,netLambda_sigma2byM_fit);
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
	mg21->GetYaxis()->SetRangeUser(0.3,12.7);
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
	mg21->GetYaxis()->SetTitleOffset(0.95);
	mg21->Draw("a");

	TLatex *text21 = new TLatex(0.4,0.75,"#splitline{net-lambda results}{and calculations}");
	text21->SetNDC();
	text21->SetTextSize(textSize/p21->GetHNDC());
	text21->Draw();

	p22->cd();

	auto *gr221 = new TGraphErrors(n_expData,energy_exp,netLambda_Ssigma_exp_0to5,(double*)0,netLambda_Ssigma_exp_0to5_staterr);
	gr221->SetMarkerSize(markerSize);
	gr221->SetMarkerColorAlpha(kOrange-3,transparency);
	gr221->SetLineColorAlpha(kOrange-3,transparency);
	gr221->SetLineWidth(lineWidth);
	gr221->SetMarkerStyle(kFullCircle);
	auto *gr222 = new TGraphErrors(n_expData,energy_exp,netLambda_Ssigma_exp_0to5,(double*)0,netLambda_Ssigma_exp_0to5_syserr);
	gr222->SetLineWidth(lineWidth);
	gr222->SetLineColorAlpha(kOrange-3,transparency);

	auto *gr223 = new TGraphErrors(n_expData,energy_exp,netLambda_Ssigma_exp_70to80,(double*)0,netLambda_Ssigma_exp_70to80_staterr);
	gr223->SetMarkerSize(markerSize);
	gr223->SetMarkerColorAlpha(kViolet-1,transparency);
	gr223->SetLineColorAlpha(kViolet-1,transparency);
	gr223->SetLineWidth(lineWidth);
	gr223->SetMarkerStyle(kFullSquare);
	auto *gr224 = new TGraphErrors(n_expData,energy_exp,netLambda_Ssigma_exp_70to80,(double*)0,netLambda_Ssigma_exp_70to80_syserr);
	gr224->SetLineWidth(lineWidth);
	gr224->SetLineColorAlpha(kViolet-1,transparency);

	auto *gr225 = new TGraph(intvl,energy,netLambda_Ssigma_hrg);
	gr225->SetLineStyle(1);
	gr225->SetLineColorAlpha(kGreen+2,transparency);
	gr225->SetLineWidth(lineWidth);
	auto *gr226 = new TGraph(n_expData,energy_exp,netLambda_Ssigma_fit);
	gr226->SetMarkerSize(markerSize);
	gr226->SetMarkerColorAlpha(kBlue,transparency);
	gr226->SetMarkerStyle(kFullTriangleUp);

	TMultiGraph  *mg22  = new TMultiGraph();

	mg22->Add(gr221,"pz");
	mg22->Add(gr222,"||");
	mg22->Add(gr223,"pz");
	mg22->Add(gr224,"||");
	mg22->Add(gr225,"l");
	mg22->Add(gr226,"p");

	mg22->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
	mg22->GetYaxis()->SetTitle("S#sigma");
	// mg22->GetXaxis()->SetLimits(5,250);
	mg22->GetYaxis()->SetRangeUser(-0.17,1.33);
	mg22->GetYaxis()->SetNdivisions(7, 2, 0, kTRUE);
	mg22->GetXaxis()->CenterTitle(true);
	mg22->GetYaxis()->CenterTitle(true);
	mg22->GetXaxis()->SetNoExponent();
	mg22->GetXaxis()->SetLabelSize(textSize/p22->GetHNDC());
	mg22->GetXaxis()->SetTitleSize(textSize/p22->GetHNDC());
	mg22->GetXaxis()->SetTitleOffset(0.75);
	mg22->GetYaxis()->SetLabelSize(textSize/p22->GetHNDC());
	mg22->GetYaxis()->SetTitleSize(textSize/p22->GetHNDC());
	mg22->GetYaxis()->SetTickLength(0.02);
	mg22->GetYaxis()->SetTitleOffset(0.95);
	mg22->Draw("a");

	TLegend *legend2 = new TLegend(0.15,0.13,0.4,0.4);
	legend2->AddEntry(gr212,"STAR Au-Au 0 - 5%","p");
	legend2->AddEntry(gr214,"STAR Au-Au 70 - 80%","p");
	legend2->AddEntry(gr211,"HRG + decay + cuts","l");
	legend2->AddEntry(gr216,"HRG + decay + cuts (fit)","p");
	legend2->SetTextSize(textSize/p22->GetHNDC());
	legend2->SetBorderSize(0);
	legend2->SetFillStyle(0);
	legend2->Draw();	

	TLatex *latex2 = new TLatex(0.6,0.65,"#splitline{0.9 < p_{T} < 2.0 (GeV/c)}{#left|y#right| #leq 0.5}");
	latex2->SetNDC();
	latex2->SetTextSize(textSize/p22->GetHNDC());
	latex2->Draw();

	TText *text22 = new TText(0.22,0.9,"STAR Collab., arXiv:2001.06419 (2020)");
	text22->SetNDC();
	text22->SetTextSize(textSize/p22->GetHNDC());
	text22->Draw();

	cl2->SaveAs("plots/netLambda_fluc.pdf");

	TCanvas *cl3 = new TCanvas("cl3","freezeout_netProton_fluc",550,400);

	gStyle->SetEndErrorSize(5);

	textSize = 0.05; lineWidth = 2;

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

	double T_ch_old[intvl], T_ch_old_err[intvl];
	double mub_ch_old[intvl];

	for (int i = 0; i < intvl; ++i)
	{
		x = mub_l+i*dmub;
		mub_ch_old[i] = x;
		T_ch_old[i] = 0.166-0.139*x*x-0.053*pow(x,4);
		T_ch_old_err[i] = sqrt(pow(0.002,2)+pow(0.016,2)*pow(x,4)+pow(0.021,2)*pow(x,8));
	}

	auto *gr31 = new TGraphErrors(intvl,mub_ch_old,T_ch_old,nullptr,T_ch_old_err);
	gr31->SetFillColorAlpha(kRed-9,transparency);
	gr31->SetFillStyle(1001);

	auto *gr32 = new TGraphErrors(n_expData,mub_ch,T_ch,mub_ch_err,T_ch_err);
	gr32->SetMarkerSize(markerSize);
	gr32->SetMarkerColorAlpha(kBlue,transparency);
	gr32->SetLineColorAlpha(kBlue,transparency);
	gr32->SetLineWidth(lineWidth);
	gr32->SetMarkerStyle(kFullFourTrianglesX);

	TMultiGraph  *mg3  = new TMultiGraph();

	mg3->Add(gr31,"3");
	mg3->Add(gr32,"p");

	mg3->GetXaxis()->SetTitle("#mu_{B} (GeV)");
	mg3->GetYaxis()->SetTitle("T (GeV)");
	mg3->GetXaxis()->SetRangeUser(0,0.43);
	mg3->GetYaxis()->SetRangeUser(0.083,0.187);
	mg3->GetXaxis()->SetNdivisions(8, 2, 0, kTRUE);
	mg3->GetYaxis()->SetNdivisions(7, 2, 0, kTRUE);
	mg3->GetXaxis()->CenterTitle(true);
	mg3->GetYaxis()->CenterTitle(true);

	mg3->GetXaxis()->SetLabelSize(textSize);
	mg3->GetXaxis()->SetTitleSize(textSize);
	mg3->GetXaxis()->SetTickLength(0.02);
	mg3->GetXaxis()->SetTitleOffset(1);

	mg3->GetYaxis()->SetLabelSize(textSize);
	mg3->GetYaxis()->SetTitleSize(textSize);
	mg3->GetYaxis()->SetTickLength(0.02);
	mg3->GetYaxis()->SetTitleOffset(1.1);
	mg3->Draw("a");

	TLegend *legend3 = new TLegend(0.15,0.15,0.6,0.35);
	legend3->AddEntry(gr32,"Chem. freezout (using fluc.)","p");
	legend3->AddEntry(gr31,"Chem. freezout (using yield)","f");
	legend3->SetTextSize(textSize);
	legend3->SetBorderSize(0);
	legend3->SetFillStyle(0);
	legend3->Draw();

	cl3->SaveAs("plots/freezeout_netProton_fluc.pdf");
}

int main()
{
	TApplication app {"app", nullptr, nullptr};
	findChemFreezeoutParam_netProton_fluc();

	cout << "\nPress Ctrl-C to exit..." << endl;

	app.Run();
	getchar();
	return 0;
}