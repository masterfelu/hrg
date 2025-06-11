#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

#include "TStyle.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TApplication.h"

#include "ThermalFunctions.h"
#include "Particles.h"

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

void plot_hadronProductionRatio()
{
	int binCount {14};	// no. of ratios

	/*---------------------------------------------------------------------------------------
	|									Data Loading Region									|
	----------------------------------------------------------------------------------------*/

	DataFileReader infile;

	infile.read("data/hadronProductionRatio.dat");
	double *particleRatio_exp {infile.getColumn(0)};
	double *particleRatio_exp_err {infile.getColumn(1)};
	infile.close();

	double energy {200};	// sqrt(snn) = 200 GeV

	infile.read("data/freezeout_netKaonLambda_fluc.dat");
	double *row_data {infile.getRow(4)};
	double T_ch_KDflucFit {row_data[1]}, T_ch_KDflucFit_err {row_data[2]};
	double mub_ch_KDflucFit {row_data[3]}, mub_ch_KDflucFit_err {row_data[4]};
	double muq_ch_KDflucFit {row_data[5]}, muq_ch_KDflucFit_err {row_data[6]};
	double mus_ch_KDflucFit {row_data[7]}, mus_ch_KDflucFit_err {row_data[8]};
	delete[] row_data;
	infile.close();


	infile.read("data/freezeout_netQP_fluc.dat");
	row_data  = infile.getRow(5);
	double T_ch_QPflucFit {row_data[1]}, T_ch_QPflucFit_err {row_data[2]};
	double mub_ch_QPflucFit {row_data[3]}, mub_ch_QPflucFit_err {row_data[4]};
	double muq_ch_QPflucFit {row_data[5]}, muq_ch_QPflucFit_err {row_data[6]};
	double mus_ch_QPflucFit {row_data[7]}, mus_ch_QPflucFit_err {row_data[8]};
	delete[] row_data;
	infile.close();

	/*---------------------------------------------------------------------------------------
	|								Yield calculation region								|
	----------------------------------------------------------------------------------------*/

	double mub_ch_yieldFit {1.308/(1+0.273*energy)};
	double muq_ch_yieldFit {-0.0202/(1+0.126*energy)};
	double mus_ch_yieldFit {0.214/(1+0.184*energy)};
	double T_ch_yieldFit {0.166-0.139*mub_ch_yieldFit*mub_ch_yieldFit-0.053*pow(mub_ch_yieldFit,4)};


	ParticleSystem hrg;
	hrg.loadDefaultData();

	hrg.setTemperature(T_ch_KDflucFit);
	hrg.setChemicalPotential(mub_ch_KDflucFit,muq_ch_KDflucFit,mus_ch_KDflucFit);

	Susceptibility1 n_KDflucFit {hrg};
	n_KDflucFit.setSphericalUniformCoordinates();
	n_KDflucFit.setIntegralLimits(0,INFINITY);
	n_KDflucFit.calculateValue();
	double ratio_KDflucFit[binCount];
	ratio_KDflucFit[0] = n_KDflucFit.getValueEach(34)/n_KDflucFit.getValueEach(33);
	ratio_KDflucFit[1] = n_KDflucFit.getValueEach(72)/n_KDflucFit.getValueEach(71);
	ratio_KDflucFit[2] = n_KDflucFit.getValueEach(191)/n_KDflucFit.getValueEach(190);
	ratio_KDflucFit[3] = n_KDflucFit.getValueEach(33)/n_KDflucFit.getValueEach(2);
	ratio_KDflucFit[4] = n_KDflucFit.getValueEach(34)/n_KDflucFit.getValueEach(2);
	ratio_KDflucFit[5] = n_KDflucFit.getValueEach(71)/n_KDflucFit.getValueEach(2);
	ratio_KDflucFit[6] = n_KDflucFit.getValueEach(72)/n_KDflucFit.getValueEach(2);
	ratio_KDflucFit[7] = n_KDflucFit.getValueEach(190)/n_KDflucFit.getValueEach(2);
	ratio_KDflucFit[8] = n_KDflucFit.getValueEach(191)/n_KDflucFit.getValueEach(2);
	ratio_KDflucFit[9] = n_KDflucFit.getValueEach(2)/n_KDflucFit.getValueEach(1);
	ratio_KDflucFit[10] = n_KDflucFit.getValueEach(24)/n_KDflucFit.getValueEach(23);
	ratio_KDflucFit[11] = n_KDflucFit.getValueEach(5)/n_KDflucFit.getValueEach(4);
	ratio_KDflucFit[12] = n_KDflucFit.getValueEach(4)/n_KDflucFit.getValueEach(2);
	ratio_KDflucFit[13] = n_KDflucFit.getValueEach(5)/n_KDflucFit.getValueEach(2);

	hrg.setTemperature(T_ch_QPflucFit);
	hrg.setChemicalPotential(mub_ch_QPflucFit,muq_ch_QPflucFit,mus_ch_QPflucFit);

	Susceptibility1 n_QPflucFit {hrg};
	n_QPflucFit.setSphericalUniformCoordinates();
	n_QPflucFit.setIntegralLimits(0,INFINITY);
	n_QPflucFit.calculateValue();
	double ratio_QPflucFit[binCount];
	ratio_QPflucFit[0] = n_QPflucFit.getValueEach(34)/n_QPflucFit.getValueEach(33);
	ratio_QPflucFit[1] = n_QPflucFit.getValueEach(72)/n_QPflucFit.getValueEach(71);
	ratio_QPflucFit[2] = n_QPflucFit.getValueEach(191)/n_QPflucFit.getValueEach(190);
	ratio_QPflucFit[3] = n_QPflucFit.getValueEach(33)/n_QPflucFit.getValueEach(2);
	ratio_QPflucFit[4] = n_QPflucFit.getValueEach(34)/n_QPflucFit.getValueEach(2);
	ratio_QPflucFit[5] = n_QPflucFit.getValueEach(71)/n_QPflucFit.getValueEach(2);
	ratio_QPflucFit[6] = n_QPflucFit.getValueEach(72)/n_QPflucFit.getValueEach(2);
	ratio_QPflucFit[7] = n_QPflucFit.getValueEach(190)/n_QPflucFit.getValueEach(2);
	ratio_QPflucFit[8] = n_QPflucFit.getValueEach(191)/n_QPflucFit.getValueEach(2);
	ratio_QPflucFit[9] = n_QPflucFit.getValueEach(2)/n_QPflucFit.getValueEach(1);
	ratio_QPflucFit[10] = n_QPflucFit.getValueEach(24)/n_QPflucFit.getValueEach(23);
	ratio_QPflucFit[11] = n_QPflucFit.getValueEach(5)/n_QPflucFit.getValueEach(4);
	ratio_QPflucFit[12] = n_QPflucFit.getValueEach(4)/n_QPflucFit.getValueEach(2);
	ratio_QPflucFit[13] = n_QPflucFit.getValueEach(5)/n_QPflucFit.getValueEach(2);

	hrg.setTemperature(T_ch_yieldFit);
	hrg.setChemicalPotential(mub_ch_yieldFit,muq_ch_yieldFit,mus_ch_yieldFit);

	Susceptibility1 n_yieldFit {hrg};
	n_yieldFit.setSphericalUniformCoordinates();
	n_yieldFit.setIntegralLimits(0,INFINITY);
	n_yieldFit.calculateValue();
	double ratio_yieldFit[binCount];
	ratio_yieldFit[0] = n_yieldFit.getValueEach(34)/n_yieldFit.getValueEach(33);
	ratio_yieldFit[1] = n_yieldFit.getValueEach(72)/n_yieldFit.getValueEach(71);
	ratio_yieldFit[2] = n_yieldFit.getValueEach(191)/n_yieldFit.getValueEach(190);
	ratio_yieldFit[3] = n_yieldFit.getValueEach(33)/n_yieldFit.getValueEach(2);
	ratio_yieldFit[4] = n_yieldFit.getValueEach(34)/n_yieldFit.getValueEach(2);
	ratio_yieldFit[5] = n_yieldFit.getValueEach(71)/n_yieldFit.getValueEach(2);
	ratio_yieldFit[6] = n_yieldFit.getValueEach(72)/n_yieldFit.getValueEach(2);
	ratio_yieldFit[7] = n_yieldFit.getValueEach(190)/n_yieldFit.getValueEach(2);
	ratio_yieldFit[8] = n_yieldFit.getValueEach(191)/n_yieldFit.getValueEach(2);
	ratio_yieldFit[9] = n_yieldFit.getValueEach(2)/n_yieldFit.getValueEach(1);
	ratio_yieldFit[10] = n_yieldFit.getValueEach(24)/n_yieldFit.getValueEach(23);
	ratio_yieldFit[11] = n_yieldFit.getValueEach(5)/n_yieldFit.getValueEach(4);
	ratio_yieldFit[12] = n_yieldFit.getValueEach(4)/n_yieldFit.getValueEach(2);
	ratio_yieldFit[13] = n_yieldFit.getValueEach(5)/n_yieldFit.getValueEach(2);


	string binName[binCount];

	binName[0] = "#color[0]{#bar{|}}#bar{#Lambda}/#Lambda";
	binName[1] = "#color[0]{#bar{|}}#bar{#Xi}/#Xi";
	binName[2] = "#color[0]{#bar{|}}#bar{#Omega}/#Omega";
	binName[3] = "#color[0]{#bar{|}}#Lambda/#pi^{-}";
	binName[4] = "#color[0]{#bar{|}}#bar{#Lambda}/#pi^{-}";
	binName[5] = "#color[0]{#bar{|}}#Xi/#pi^{-}";
	binName[6] = "#color[0]{#bar{|}}#bar{#Xi}/#pi^{-}";
	binName[7] = "#color[0]{#bar{|}}#Omega/#pi^{-}";
	binName[8] = "#color[0]{#bar{|}}#bar{#Omega}/#pi^{-}";
	binName[9] = "#color[0]{#bar{|}}#pi^{-}/#pi^{+}";
	binName[10] = "#color[0]{#bar{|}}#bar{p}/p";
	binName[11] = "#color[0]{#bar{|}}K^{-}/K^{+}";
	binName[12] = "#color[0]{#bar{|}}K^{+}/#pi^{-}";
	binName[13] = "#color[0]{#bar{|}}K^{-}/#pi^{-}";

	/*---------------------------------------------------------------------------------------
	|									Data Ploting region									|
	----------------------------------------------------------------------------------------*/

	double topMargin = 0.01, bottomMargin = 0.1, leftMargin = 0.08, rightMargin = 0.005;
	double textSize = 0.06;
	double markerSize = 1.3;
	double lineWidth = 2;
	double tickLength = 0.02;
	double opacity = 1;

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	gStyle->SetEndErrorSize(7);

	TCanvas *canvas = new TCanvas("c1","Hadron Production Ratio",700,450);

	canvas->SetTopMargin(topMargin);
	canvas->SetBottomMargin(bottomMargin);
	canvas->SetLeftMargin(leftMargin);
	canvas->SetRightMargin(rightMargin);

	canvas->SetGridx();

	canvas->SetTicky();
	canvas->SetLogy();

	TH1D *haxis = new TH1D("hist_axis","axis",binCount,0,binCount);

	haxis->GetXaxis()->SetTickLength(0);
	haxis->GetXaxis()->SetLabelSize(textSize);

	haxis->GetYaxis()->SetTitle("hadron production ratio");
	haxis->GetYaxis()->CenterTitle(true);
	haxis->GetYaxis()->SetLabelSize(textSize/2);
	haxis->GetYaxis()->SetTitleSize(textSize/2);
	haxis->GetYaxis()->SetTickLength(0.01);
	haxis->GetYaxis()->SetRangeUser(8e-5,8);

	double binCenter[binCount];
	double xError[binCount];

	for (int i = 0; i < binCount; ++i)
	{
		xError[i] = 0.8*haxis->GetBinWidth(i+1)/2;
		binCenter[i] = haxis->GetBinCenter(i+1);
		haxis->GetXaxis()->SetBinLabel(i+1,binName[i].c_str());
	}

	haxis->Draw("axis");

	TGraphErrors *gr1 = new TGraphErrors(binCount,binCenter,particleRatio_exp,0,particleRatio_exp_err);
	gr1->SetMarkerStyle(kFullCircle);
	gr1->SetMarkerColorAlpha(kGray+3,opacity);
	gr1->SetMarkerSize(markerSize);
	gr1->SetLineStyle(kSolid);
	gr1->SetLineColorAlpha(kGray+3,opacity);
	gr1->SetLineWidth(lineWidth);

	delete[] particleRatio_exp; delete[] particleRatio_exp_err;

	TGraphErrors *gr2 = new TGraphErrors(binCount,binCenter,ratio_yieldFit,xError,0);
	gr2->SetMarkerSize(0);
	gr2->SetLineStyle(kSolid);
	gr2->SetLineColorAlpha(kBlue,opacity);
	gr2->SetLineWidth(lineWidth);

	TGraphErrors *gr3 = new TGraphErrors(binCount,binCenter,ratio_KDflucFit,xError,0);
	gr3->SetMarkerSize(0);
	gr3->SetLineStyle(6);
	gr3->SetLineColorAlpha(kRed,opacity);
	gr3->SetLineWidth(lineWidth);

	TGraphErrors *gr4 = new TGraphErrors(binCount,binCenter,ratio_QPflucFit,xError,0);
	gr4->SetMarkerSize(0);
	gr4->SetLineStyle(7);
	gr4->SetLineColorAlpha(kSpring-7,opacity);
	gr4->SetLineWidth(lineWidth);

	TMultiGraph *mg = new TMultiGraph();

	mg->Add(gr1,"p");
	mg->Add(gr2,"p");
	mg->Add(gr3,"p");
	mg->Add(gr4,"p");

	mg->Draw();

	gPad->RedrawAxis();
	gPad->RedrawAxis("G");

	const char *label1 = "STAR Au-Au #sqrt{s_{NN}} = 200 GeV";

	char *label2 = new char[100];
	sprintf(label2,"T_{ch} = %.1f MeV, #mu_{B} = %.1f MeV (Yield fit)",T_ch_yieldFit*1e3, mub_ch_yieldFit*1e3);

	char *label3 = new char[100];
	sprintf(label3,"T_{ch} = %.1f MeV, #mu_{B} = %.1f MeV (K#Lambda fluc. fit)",T_ch_KDflucFit*1e3, mub_ch_KDflucFit*1e3);

	char *label4 = new char[100];
	sprintf(label4,"T_{ch} = %.1f MeV, #mu_{B} = %.1f MeV (QP fluc. fit)",T_ch_QPflucFit*1e3, mub_ch_QPflucFit*1e3);

	TLegend *legend = new TLegend(0.11,0.12,0.5,0.3);
	legend->AddEntry(gr1,label1,"p");
	legend->AddEntry(gr2,label2,"l");
	legend->AddEntry(gr3,label3,"l");
	legend->AddEntry(gr4,label4,"l");
	// legend->SetMargin(0.2);
	legend->SetTextAlign(12);
	legend->Draw();

	canvas->SaveAs("plots/hadronProductionRatio.pdf");
}

int main()
{
	TApplication app {"app", nullptr, nullptr};
	plot_hadronProductionRatio();

	cout << "Press Ctrl-C to exit..." << endl;
	app.Run();
	return 0;
}