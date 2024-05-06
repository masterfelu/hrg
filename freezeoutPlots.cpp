#include <iostream>
#include <cmath>
#include <fstream>

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

void plotFreezeoutAll()
{
	double mub_l = 0, mub_u = 1;
	double T_l = 0.1, T_u = 0.18;

	int intvl = 100;

	double dT = (T_u-T_l)/(intvl-1);
	double dmub = (mub_u-mub_l)/(intvl-1);

	double T_fo_old[intvl], T_fo_old_err[intvl];
	double mub_fo_old[intvl];

	double x;

	for (int i = 0; i < intvl; ++i)
	{
		x = mub_l+i*dmub;
		mub_fo_old[i] = x;
		T_fo_old[i] = 0.166-0.139*x*x-0.053*pow(x,4);
		T_fo_old_err[i] = sqrt(pow(0.002,2)+pow(0.016,2)*pow(x,4)+pow(0.021,2)*pow(x,8));
	}

	int n_fo_PQ_alba, n_fo_PQ, n_fo_KL, n_fo_PL, n_fo_KQ;

	double *T_fo_PQ_alba, *T_fo_err_PQ_alba;
	double *mub_fo_PQ_alba, *mub_fo_err_PQ_alba;

	double *T_fo_PQ, *T_fo_err_PQ;
	double *mub_fo_PQ, *mub_fo_err_PQ;

	double *T_fo_KL, *T_fo_err_KL;
	double *mub_fo_KL, *mub_fo_err_KL;

	double *T_fo_PL, *T_fo_err_PL;
	double *mub_fo_PL, *mub_fo_err_PL;

	double *T_fo_KQ, *T_fo_err_KQ;
	double *mub_fo_KQ, *mub_fo_err_KQ;

	DataFileReader infile;

	infile.read("data/freezeout_netQP_fluc_alba.dat");
	if(infile.isDataLoaded())
	{
		n_fo_PQ_alba = infile.nRows();
		T_fo_PQ_alba = infile.getColumn(1); T_fo_err_PQ_alba = infile.getColumn(2);
		mub_fo_PQ_alba = infile.getColumn(3); mub_fo_err_PQ_alba = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netQP_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_PQ = infile.nRows();
		T_fo_PQ = infile.getColumn(1); T_fo_err_PQ = infile.getColumn(2);
		mub_fo_PQ = infile.getColumn(3); mub_fo_err_PQ = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netKaonLambda_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_KL = infile.nRows();
		T_fo_KL = infile.getColumn(1); T_fo_err_KL = infile.getColumn(2);
		mub_fo_KL = infile.getColumn(3); mub_fo_err_KL = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netProtonLambda_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_PL = infile.nRows();
		T_fo_PL = infile.getColumn(1); T_fo_err_PL = infile.getColumn(2);
		mub_fo_PL = infile.getColumn(3); mub_fo_err_PL = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netKaonCharge_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_KQ = infile.nRows();
		T_fo_KQ = infile.getColumn(1); T_fo_err_KQ = infile.getColumn(2);
		mub_fo_KQ = infile.getColumn(3); mub_fo_err_KQ = infile.getColumn(4);
		infile.close();
	}

	double topMargin = 0.01;
	double bottomMargin = 0.05;
	double rightMargin = 0.02;
	double leftMargin = 0.1;

	double textSize = 0.05;
	double markerSize = 2;
	double lineWidth = 2;
	double transparency = 0.75;

	gStyle->SetEndErrorSize(5);

	TCanvas *cv1 = new TCanvas("cv1","freezeout curve",550,400);

	TPad *pd1 = new TPad("pd1", "pd1",0,0,1,1);
	pd1->SetBottomMargin(0.12);
	pd1->SetTopMargin(0.01);
	pd1->SetRightMargin(0.01);
	pd1->SetLeftMargin(0.11);
	pd1->Draw();

	pd1->cd();

	auto *gr11 = new TGraphErrors(intvl,mub_fo_old,T_fo_old,nullptr,T_fo_old_err);
	gr11->SetFillColorAlpha(kRed-9,transparency);
	gr11->SetFillStyle(1001);

	auto *gr12 = new TGraphErrors(n_fo_PQ_alba,mub_fo_PQ_alba,T_fo_PQ_alba,mub_fo_err_PQ_alba,T_fo_err_PQ_alba);
	gr12->SetFillColorAlpha(kYellow+1,transparency);
	gr12->SetFillStyle(1001);

	auto *gr13 = new TGraphErrors(n_fo_KL,mub_fo_KL,T_fo_KL,mub_fo_err_KL,T_fo_err_KL);
	gr13->SetFillColorAlpha(kGreen+2,transparency);
	gr13->SetFillStyle(1001);
	gr13->SetMarkerSize(markerSize);
	gr13->SetMarkerColorAlpha(kGreen+2,transparency);
	gr13->SetLineColorAlpha(kGreen+2,transparency);
	gr13->SetLineWidth(lineWidth);
	gr13->SetMarkerStyle(kFullFourTrianglesX);

	auto *gr14 = new TGraphErrors(n_fo_PL,mub_fo_PL,T_fo_PL,mub_fo_err_PL,T_fo_err_PL);
	gr14->SetFillColorAlpha(kViolet-1,transparency);
	gr14->SetFillStyle(1001);
	gr14->SetMarkerSize(markerSize);
	gr14->SetMarkerColorAlpha(kViolet-1,transparency);
	gr14->SetLineColorAlpha(kViolet-1,transparency);
	gr14->SetLineWidth(lineWidth);
	gr14->SetMarkerStyle(kFullFourTrianglesX);

	auto *gr15 = new TGraphErrors(n_fo_KQ,mub_fo_KQ,T_fo_KQ,mub_fo_err_KQ,T_fo_err_KQ);
	gr15->SetFillColorAlpha(kMagenta,transparency);
	gr15->SetFillStyle(1001);
	gr15->SetMarkerSize(markerSize);
	gr15->SetMarkerColorAlpha(kMagenta,transparency);
	gr15->SetLineColorAlpha(kMagenta,transparency);
	gr15->SetLineWidth(lineWidth);
	gr15->SetMarkerStyle(kFullFourTrianglesX);

	auto *grtemp1 = new TGraphErrors(n_fo_PQ,mub_fo_PQ,T_fo_PQ,mub_fo_err_PQ,T_fo_err_PQ);
	grtemp1->SetMarkerSize(markerSize);
	grtemp1->SetMarkerColorAlpha(kBlue,transparency);
	grtemp1->SetLineColorAlpha(kBlue,transparency);
	grtemp1->SetLineWidth(lineWidth);
	grtemp1->SetMarkerStyle(kFullFourTrianglesX);

	TMultiGraph  *mg1  = new TMultiGraph();

	mg1->Add(gr11,"3");
	mg1->Add(gr12,"3");
	mg1->Add(gr13,"p");
	mg1->Add(gr14,"p");
	mg1->Add(gr15,"p");
	mg1->Add(grtemp1,"p");

	mg1->GetXaxis()->SetTitle("#mu_{B} (GeV)");
	mg1->GetYaxis()->SetTitle("T (GeV)");
	mg1->GetXaxis()->SetRangeUser(0,0.43);
	mg1->GetYaxis()->SetRangeUser(0.083,0.187);
	mg1->GetXaxis()->SetNdivisions(8, 2, 0, kTRUE);
	mg1->GetYaxis()->SetNdivisions(7, 2, 0, kTRUE);
	mg1->GetXaxis()->CenterTitle(true);
	mg1->GetYaxis()->CenterTitle(true);

	mg1->GetXaxis()->SetLabelSize(textSize);
	mg1->GetXaxis()->SetTitleSize(textSize);
	mg1->GetXaxis()->SetTickLength(0.02);
	mg1->GetXaxis()->SetTitleOffset(1);

	mg1->GetYaxis()->SetLabelSize(textSize);
	mg1->GetYaxis()->SetTitleSize(textSize);
	mg1->GetYaxis()->SetTickLength(0.02);
	mg1->GetYaxis()->SetTitleOffset(1.1);
	mg1->Draw("a");

	TLegend *lg1 = new TLegend(0.15,0.15,0.6,0.45);

	lg1->AddEntry(gr11,"Yield fit, Cleymans et al.","f");
	lg1->AddEntry(gr12,"pq fluc. fit, Alba et al.","f");
	lg1->AddEntry(gr13,"k#Lambda fluc. fit, this work","p");
	lg1->AddEntry(gr14,"p#Lambda fluc. fit, this work","p");
	lg1->AddEntry(gr15,"kq fluc. fit, this work","p");
	lg1->AddEntry(grtemp1,"pq fluc. fit, this work","p");

	lg1->SetTextSize(textSize);
	lg1->SetBorderSize(0);
	lg1->SetFillStyle(0);
	lg1->Draw();

	cv1->SaveAs("plots/freezeoutAll.pdf");
}

void plotFreezeoutMass()
{
	double mub_l = 0, mub_u = 1;
	double T_l = 0.1, T_u = 0.18;

	int intvl = 100;

	double dT = (T_u-T_l)/(intvl-1);
	double dmub = (mub_u-mub_l)/(intvl-1);

	double T_fo_old[intvl], T_fo_old_err[intvl];
	double mub_fo_old[intvl];

	double x;

	for (int i = 0; i < intvl; ++i)
	{
		x = mub_l+i*dmub;
		mub_fo_old[i] = x;
		T_fo_old[i] = 0.166-0.139*x*x-0.053*pow(x,4);
		T_fo_old_err[i] = sqrt(pow(0.002,2)+pow(0.016,2)*pow(x,4)+pow(0.021,2)*pow(x,8));
	}

	int n_fo_PQ_alba, n_fo_PQ, n_fo_KL, n_fo_PL, n_fo_KQ;

	double *energy_PQ_alba;
	double *T_fo_PQ_alba, *T_fo_err_PQ_alba;
	double *mub_fo_PQ_alba, *mub_fo_err_PQ_alba;

	double *energy_PQ;
	double *T_fo_PQ, *T_fo_err_PQ;
	double *mub_fo_PQ, *mub_fo_err_PQ;

	double *energy_KL;
	double *T_fo_KL, *T_fo_err_KL;
	double *mub_fo_KL, *mub_fo_err_KL;

	double *energy_PL;
	double *T_fo_PL, *T_fo_err_PL;
	double *mub_fo_PL, *mub_fo_err_PL;

	double *energy_KQ;
	double *T_fo_KQ, *T_fo_err_KQ;
	double *mub_fo_KQ, *mub_fo_err_KQ;

	DataFileReader infile;

	infile.read("data/freezeout_netQP_fluc_alba.dat");
	if(infile.isDataLoaded())
	{
		n_fo_PQ_alba = infile.nRows(); energy_PQ_alba = infile.getColumn(0);
		T_fo_PQ_alba = infile.getColumn(1); T_fo_err_PQ_alba = infile.getColumn(2);
		mub_fo_PQ_alba = infile.getColumn(3); mub_fo_err_PQ_alba = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netQP_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_PQ = infile.nRows();  energy_PQ = infile.getColumn(0);
		T_fo_PQ = infile.getColumn(1); T_fo_err_PQ = infile.getColumn(2);
		mub_fo_PQ = infile.getColumn(3); mub_fo_err_PQ = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netKaonLambda_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_KL = infile.nRows();
		T_fo_KL = infile.getColumn(1); T_fo_err_KL = infile.getColumn(2);
		mub_fo_KL = infile.getColumn(3); mub_fo_err_KL = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netProtonLambda_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_PL = infile.nRows();
		T_fo_PL = infile.getColumn(1); T_fo_err_PL = infile.getColumn(2);
		mub_fo_PL = infile.getColumn(3); mub_fo_err_PL = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netKaonCharge_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_KQ = infile.nRows();
		T_fo_KQ = infile.getColumn(1); T_fo_err_KQ = infile.getColumn(2);
		mub_fo_KQ = infile.getColumn(3); mub_fo_err_KQ = infile.getColumn(4);
		infile.close();
	}

	double topMargin = 0.01;
	double bottomMargin = 0.12;
	double rightMargin = 0.01;
	double leftMargin = 0.11;

	double textSize = 0.05;
	double markerSize = 2;
	double lineWidth = 2;
	double transparency = 0.75;

	gStyle->SetEndErrorSize(5);

	TCanvas *cv1 = new TCanvas("cv1_Mass","freezeout curve",550,400);

	TPad *pd1 = new TPad("pd1_Mass", "pd1",0,0,1,1);
	pd1->SetBottomMargin(bottomMargin);
	pd1->SetTopMargin(topMargin);
	pd1->SetRightMargin(rightMargin);
	pd1->SetLeftMargin(leftMargin);
	pd1->Draw();

	pd1->cd();

	auto *gr11 = new TGraphErrors(intvl,mub_fo_old,T_fo_old,nullptr,T_fo_old_err);
	gr11->SetFillColorAlpha(kRed-9,transparency);
	gr11->SetFillStyle(1001);

	auto *gr12 = new TGraphErrors(n_fo_PQ_alba,mub_fo_PQ_alba,T_fo_PQ_alba,mub_fo_err_PQ_alba,T_fo_err_PQ_alba);
	gr12->SetFillColorAlpha(kYellow+1,transparency);
	gr12->SetFillStyle(1001);

	auto *gr13 = new TGraphErrors(n_fo_KL,mub_fo_KL,T_fo_KL,mub_fo_err_KL,T_fo_err_KL);
	gr13->SetFillColorAlpha(kGreen+2,transparency);
	gr13->SetFillStyle(1001);
	gr13->SetMarkerSize(markerSize);
	gr13->SetMarkerColorAlpha(kGreen+2,transparency);
	gr13->SetLineColorAlpha(kGreen+2,transparency);
	gr13->SetLineWidth(lineWidth);
	gr13->SetMarkerStyle(kFullFourTrianglesX);

	auto *gr14 = new TGraphErrors(n_fo_PL,mub_fo_PL,T_fo_PL,mub_fo_err_PL,T_fo_err_PL);
	gr14->SetFillColorAlpha(kViolet-1,transparency);
	gr14->SetFillStyle(1001);
	gr14->SetMarkerSize(markerSize);
	gr14->SetMarkerColorAlpha(kViolet-1,transparency);
	gr14->SetLineColorAlpha(kViolet-1,transparency);
	gr14->SetLineWidth(lineWidth);
	gr14->SetMarkerStyle(kFullFourTrianglesX);

	auto *gr15 = new TGraphErrors(n_fo_KQ,mub_fo_KQ,T_fo_KQ,mub_fo_err_KQ,T_fo_err_KQ);
	gr15->SetFillColorAlpha(kMagenta,transparency);
	gr15->SetFillStyle(1001);
	gr15->SetMarkerSize(markerSize);
	gr15->SetMarkerColorAlpha(kMagenta,transparency);
	gr15->SetLineColorAlpha(kMagenta,transparency);
	gr15->SetLineWidth(lineWidth);
	gr15->SetMarkerStyle(kFullFourTrianglesX);

	auto *grtemp1 = new TGraphErrors(n_fo_PQ,mub_fo_PQ,T_fo_PQ,mub_fo_err_PQ,T_fo_err_PQ);
	grtemp1->SetMarkerSize(markerSize);
	grtemp1->SetMarkerColorAlpha(kBlue,transparency);
	grtemp1->SetLineColorAlpha(kBlue,transparency);
	grtemp1->SetLineWidth(lineWidth);
	grtemp1->SetMarkerStyle(kFullFourTrianglesX);

	TMultiGraph  *mg1  = new TMultiGraph();

	mg1->Add(gr11,"3");
	// mg1->Add(gr12,"3"); // Alba plots
	// mg1->Add(gr13,"p");
	mg1->Add(gr14,"pl");
	mg1->Add(gr15,"pl");
	// mg1->Add(grtemp1,"p");

	mg1->GetXaxis()->SetTitle("#mu_{B} (GeV)");
	mg1->GetYaxis()->SetTitle("T (GeV)");
	mg1->GetXaxis()->SetRangeUser(0,0.43);
	mg1->GetYaxis()->SetRangeUser(0.083,0.187);
	mg1->GetXaxis()->SetNdivisions(8, 2, 0, kTRUE);
	mg1->GetYaxis()->SetNdivisions(7, 2, 0, kTRUE);
	mg1->GetXaxis()->CenterTitle(true);
	mg1->GetYaxis()->CenterTitle(true);

	mg1->GetXaxis()->SetLabelSize(textSize);
	mg1->GetXaxis()->SetTitleSize(textSize);
	mg1->GetXaxis()->SetTickLength(0.02);
	mg1->GetXaxis()->SetTitleOffset(1);

	mg1->GetYaxis()->SetLabelSize(textSize);
	mg1->GetYaxis()->SetTitleSize(textSize);
	mg1->GetYaxis()->SetTickLength(0.02);
	mg1->GetYaxis()->SetTitleOffset(1.1);
	mg1->Draw("a");

	TLegend *lg1 = new TLegend(0.15,0.15,0.6,0.45);

	lg1->AddEntry(gr11,"Yield fit, Cleymans et al.","f");
	// lg1->AddEntry(gr12,"pq fluc. fit, Alba et al.","f");
	// lg1->AddEntry(gr13,"k#Lambda fluc. fit, this work","p");
	lg1->AddEntry(gr14,"p#Lambda fluc. fit, this work","pl");
	lg1->AddEntry(gr15,"kq fluc. fit, this work","pl");
	// lg1->AddEntry(grtemp1,"pq fluc. fit, this work","p");

	lg1->SetTextSize(textSize);
	lg1->SetBorderSize(0);
	lg1->SetFillStyle(0);
	lg1->Draw();

	cv1->SaveAs("plots/freezeoutMass.pdf");
}

void plotFreezeoutCharge()
{
	double mub_l = 0, mub_u = 1;
	double T_l = 0.1, T_u = 0.18;

	int intvl = 100;

	double dT = (T_u-T_l)/(intvl-1);
	double dmub = (mub_u-mub_l)/(intvl-1);

	double T_fo_old[intvl], T_fo_old_err[intvl];
	double mub_fo_old[intvl];

	double x;

	for (int i = 0; i < intvl; ++i)
	{
		x = mub_l+i*dmub;
		mub_fo_old[i] = x;
		T_fo_old[i] = 0.166-0.139*x*x-0.053*pow(x,4);
		T_fo_old_err[i] = sqrt(pow(0.002,2)+pow(0.016,2)*pow(x,4)+pow(0.021,2)*pow(x,8));
	}

	int n_fo_PQ_alba, n_fo_PQ, n_fo_KL, n_fo_PL, n_fo_KQ;

	double *T_fo_PQ_alba, *T_fo_err_PQ_alba;
	double *mub_fo_PQ_alba, *mub_fo_err_PQ_alba;

	double *T_fo_PQ, *T_fo_err_PQ;
	double *mub_fo_PQ, *mub_fo_err_PQ;

	double *T_fo_KL, *T_fo_err_KL;
	double *mub_fo_KL, *mub_fo_err_KL;

	double *T_fo_PL, *T_fo_err_PL;
	double *mub_fo_PL, *mub_fo_err_PL;

	double *T_fo_KQ, *T_fo_err_KQ;
	double *mub_fo_KQ, *mub_fo_err_KQ;

	DataFileReader infile;

	infile.read("data/freezeout_netQP_fluc_alba.dat");
	if(infile.isDataLoaded())
	{
		n_fo_PQ_alba = infile.nRows();
		T_fo_PQ_alba = infile.getColumn(1); T_fo_err_PQ_alba = infile.getColumn(2);
		mub_fo_PQ_alba = infile.getColumn(3); mub_fo_err_PQ_alba = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netQP_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_PQ = infile.nRows();
		T_fo_PQ = infile.getColumn(1); T_fo_err_PQ = infile.getColumn(2);
		mub_fo_PQ = infile.getColumn(3); mub_fo_err_PQ = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netKaonLambda_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_KL = infile.nRows();
		T_fo_KL = infile.getColumn(1); T_fo_err_KL = infile.getColumn(2);
		mub_fo_KL = infile.getColumn(3); mub_fo_err_KL = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netProtonLambda_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_PL = infile.nRows();
		T_fo_PL = infile.getColumn(1); T_fo_err_PL = infile.getColumn(2);
		mub_fo_PL = infile.getColumn(3); mub_fo_err_PL = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netKaonCharge_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_KQ = infile.nRows();
		T_fo_KQ = infile.getColumn(1); T_fo_err_KQ = infile.getColumn(2);
		mub_fo_KQ = infile.getColumn(3); mub_fo_err_KQ = infile.getColumn(4);
		infile.close();
	}

	double topMargin = 0.01;
	double bottomMargin = 0.05;
	double rightMargin = 0.02;
	double leftMargin = 0.1;

	double textSize = 0.05;
	double markerSize = 2;
	double lineWidth = 2;
	double transparency = 0.75;

	gStyle->SetEndErrorSize(5);

	TCanvas *cv1 = new TCanvas("cv1_Charge","freezeout curve",550,400);

	TPad *pd1 = new TPad("pd1_Charge", "pd1",0,0,1,1);
	pd1->SetBottomMargin(0.12);
	pd1->SetTopMargin(0.01);
	pd1->SetRightMargin(0.01);
	pd1->SetLeftMargin(0.11);
	pd1->Draw();

	pd1->cd();

	auto *gr11 = new TGraphErrors(intvl,mub_fo_old,T_fo_old,nullptr,T_fo_old_err);
	gr11->SetFillColorAlpha(kRed-9,transparency);
	gr11->SetFillStyle(1001);

	auto *gr12 = new TGraphErrors(n_fo_PQ_alba,mub_fo_PQ_alba,T_fo_PQ_alba,mub_fo_err_PQ_alba,T_fo_err_PQ_alba);
	gr12->SetFillColorAlpha(kYellow+1,transparency);
	gr12->SetFillStyle(1001);

	auto *gr13 = new TGraphErrors(n_fo_KL,mub_fo_KL,T_fo_KL,mub_fo_err_KL,T_fo_err_KL);
	gr13->SetFillColorAlpha(kGreen+2,transparency);
	gr13->SetFillStyle(1001);
	gr13->SetMarkerSize(markerSize);
	gr13->SetMarkerColorAlpha(kGreen+2,transparency);
	gr13->SetLineColorAlpha(kGreen+2,transparency);
	gr13->SetLineWidth(lineWidth);
	gr13->SetMarkerStyle(kFullFourTrianglesX);

	auto *gr14 = new TGraphErrors(n_fo_PL,mub_fo_PL,T_fo_PL,mub_fo_err_PL,T_fo_err_PL);
	gr14->SetFillColorAlpha(kViolet-1,transparency);
	gr14->SetFillStyle(1001);
	gr14->SetMarkerSize(markerSize);
	gr14->SetMarkerColorAlpha(kViolet-1,transparency);
	gr14->SetLineColorAlpha(kViolet-1,transparency);
	gr14->SetLineWidth(lineWidth);
	gr14->SetMarkerStyle(kFullFourTrianglesX);

	auto *gr15 = new TGraphErrors(n_fo_KQ,mub_fo_KQ,T_fo_KQ,mub_fo_err_KQ,T_fo_err_KQ);
	gr15->SetFillColorAlpha(kMagenta,transparency);
	gr15->SetFillStyle(1001);
	gr15->SetMarkerSize(markerSize);
	gr15->SetMarkerColorAlpha(kMagenta,transparency);
	gr15->SetLineColorAlpha(kMagenta,transparency);
	gr15->SetLineWidth(lineWidth);
	gr15->SetMarkerStyle(kFullFourTrianglesX);

	auto *grtemp1 = new TGraphErrors(n_fo_PQ,mub_fo_PQ,T_fo_PQ,mub_fo_err_PQ,T_fo_err_PQ);
	grtemp1->SetMarkerSize(markerSize);
	grtemp1->SetMarkerColorAlpha(kBlue,transparency);
	grtemp1->SetLineColorAlpha(kBlue,transparency);
	grtemp1->SetLineWidth(lineWidth);
	grtemp1->SetMarkerStyle(kFullFourTrianglesX);

	TMultiGraph  *mg1  = new TMultiGraph();

	mg1->Add(gr11,"3");
	// mg1->Add(gr12,"3"); // alba plot
	mg1->Add(gr13,"p");
	// mg1->Add(gr14,"pl");
	// mg1->Add(gr15,"pl");
	mg1->Add(grtemp1,"p");

	mg1->GetXaxis()->SetTitle("#mu_{B} (GeV)");
	mg1->GetYaxis()->SetTitle("T (GeV)");
	mg1->GetXaxis()->SetRangeUser(0,0.43);
	mg1->GetYaxis()->SetRangeUser(0.083,0.187);
	mg1->GetXaxis()->SetNdivisions(8, 2, 0, kTRUE);
	mg1->GetYaxis()->SetNdivisions(7, 2, 0, kTRUE);
	mg1->GetXaxis()->CenterTitle(true);
	mg1->GetYaxis()->CenterTitle(true);

	mg1->GetXaxis()->SetLabelSize(textSize);
	mg1->GetXaxis()->SetTitleSize(textSize);
	mg1->GetXaxis()->SetTickLength(0.02);
	mg1->GetXaxis()->SetTitleOffset(1);

	mg1->GetYaxis()->SetLabelSize(textSize);
	mg1->GetYaxis()->SetTitleSize(textSize);
	mg1->GetYaxis()->SetTickLength(0.02);
	mg1->GetYaxis()->SetTitleOffset(1.1);
	mg1->Draw("a");

	TLegend *lg1 = new TLegend(0.15,0.15,0.6,0.45);

	lg1->AddEntry(gr11,"Yield fit, Cleymans et al.","f");
	// lg1->AddEntry(gr12,"pq fluc. fit, Alba et al.","f");
	lg1->AddEntry(gr13,"k#Lambda fluc. fit, this work","p");
	// lg1->AddEntry(gr14,"p#Lambda fluc. fit, this work","pl");
	// lg1->AddEntry(gr15,"kq fluc. fit, this work","pl");
	lg1->AddEntry(grtemp1,"pq fluc. fit, this work","p");

	lg1->SetTextSize(textSize);
	lg1->SetBorderSize(0);
	lg1->SetFillStyle(0);
	lg1->Draw();

	cv1->SaveAs("plots/freezeoutCharge.pdf");
}

void plotFreezeoutIndiFit()
{
	double mub_l = 0, mub_u = 1;
	double T_l = 0.1, T_u = 0.18;

	int intvl = 100;

	double dT = (T_u-T_l)/(intvl-1);
	double dmub = (mub_u-mub_l)/(intvl-1);

	double T_fo_old[intvl], T_fo_old_err[intvl];
	double mub_fo_old[intvl];

	double x;

	for (int i = 0; i < intvl; ++i)
	{
		x = mub_l+i*dmub;
		mub_fo_old[i] = x;
		T_fo_old[i] = 0.166-0.139*x*x-0.053*pow(x,4);
		T_fo_old_err[i] = sqrt(pow(0.002,2)+pow(0.016,2)*pow(x,4)+pow(0.021,2)*pow(x,8));
	}

	int n_fo_P, n_fo_K, n_fo_L;

	double *T_fo_P, *T_fo_err_P;
	double *mub_fo_P, *mub_fo_err_P;

	double *T_fo_K, *T_fo_err_K;
	double *mub_fo_K, *mub_fo_err_K;

	double *T_fo_L, *T_fo_err_L;
	double *mub_fo_L, *mub_fo_err_L;

	DataFileReader infile;

	infile.read("data/freezeout_netProton_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_P = infile.nRows();
		T_fo_P = infile.getColumn(1); T_fo_err_P = infile.getColumn(2);
		mub_fo_P = infile.getColumn(3); mub_fo_err_P = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netKaon_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_K = infile.nRows();
		T_fo_K = infile.getColumn(1); T_fo_err_K = infile.getColumn(2);
		mub_fo_K = infile.getColumn(3); mub_fo_err_K = infile.getColumn(4);
		infile.close();
	}

	infile.read("data/freezeout_netLambda_fluc.dat");
	if(infile.isDataLoaded())
	{
		n_fo_L = infile.nRows();
		T_fo_L = infile.getColumn(1); T_fo_err_L = infile.getColumn(2);
		mub_fo_L = infile.getColumn(3); mub_fo_err_L = infile.getColumn(4);
		infile.close();
	}

	double topMargin = 0.01;
	double bottomMargin = 0.05;
	double rightMargin = 0.02;
	double leftMargin = 0.1;

	double textSize = 0.05;
	double markerSize = 2.5;
	double lineWidth = 2;
	double opaque = 1, transparent = 0.75;

	gStyle->SetEndErrorSize(5);

	int Color_old = TColor::GetColor("#8a8a8a");
	int Color_P = TColor::GetColor("#e31e1e");
	int Color_K = TColor::GetColor("#11a600");
	int Color_L = TColor::GetColor("#0835c7");

	int Marker_P = kFullFourTrianglesX;
	int Marker_K = kFullCrossX;
	int Marker_L = kFullFourTrianglesPlus;	
	

	TCanvas *cv1 = new TCanvas("cv1_IndiFit","freezeout curve",550,400);

	TPad *pd1 = new TPad("pd1_IndiFit", "pd1",0,0,1,1);
	pd1->SetBottomMargin(0.12);
	pd1->SetTopMargin(0.01);
	pd1->SetRightMargin(0.01);
	pd1->SetLeftMargin(0.11);
	pd1->Draw();

	pd1->cd();

	auto *gr11 = new TGraphErrors(intvl,mub_fo_old,T_fo_old,nullptr,T_fo_old_err);
	gr11->SetFillColorAlpha(Color_old,transparent);
	gr11->SetFillStyle(1001);

	auto *gr12 = new TGraphErrors(n_fo_P,mub_fo_P,T_fo_P,mub_fo_err_P,T_fo_err_P);
	gr12->SetFillColorAlpha(Color_P,opaque);
	gr12->SetFillStyle(1001);
	gr12->SetMarkerSize(markerSize);
	gr12->SetMarkerColorAlpha(Color_P,opaque);
	gr12->SetLineColorAlpha(Color_P,opaque);
	gr12->SetLineWidth(lineWidth);
	gr12->SetMarkerStyle(Marker_P);

	auto *gr13 = new TGraphErrors(n_fo_K,mub_fo_K,T_fo_K,mub_fo_err_K,T_fo_err_K);
	gr13->SetFillColorAlpha(Color_K,opaque);
	gr13->SetFillStyle(1001);
	gr13->SetMarkerSize(markerSize);
	gr13->SetMarkerColorAlpha(Color_K,opaque);
	gr13->SetLineColorAlpha(Color_K,opaque);
	gr13->SetLineWidth(lineWidth);
	gr13->SetMarkerStyle(Marker_K);

	auto *gr14 = new TGraphErrors(n_fo_L,mub_fo_L,T_fo_L,mub_fo_err_L,T_fo_err_L);
	gr14->SetFillColorAlpha(Color_L,opaque);
	gr14->SetFillStyle(1001);
	gr14->SetMarkerSize(markerSize);
	gr14->SetMarkerColorAlpha(Color_L,opaque);
	gr14->SetLineColorAlpha(Color_L,opaque);
	gr14->SetLineWidth(lineWidth);
	gr14->SetMarkerStyle(Marker_L);

	TMultiGraph  *mg1  = new TMultiGraph();

	mg1->Add(gr11,"3");
	mg1->Add(gr12,"p");
	mg1->Add(gr13,"p");
	mg1->Add(gr14,"p");

	mg1->GetXaxis()->SetTitle("#mu_{B} (GeV)");
	mg1->GetYaxis()->SetTitle("T (GeV)");
	mg1->GetXaxis()->SetRangeUser(0,0.43);
	mg1->GetYaxis()->SetRangeUser(0.083,0.187);
	mg1->GetXaxis()->SetNdivisions(8, 2, 0, kTRUE);
	mg1->GetYaxis()->SetNdivisions(7, 2, 0, kTRUE);
	mg1->GetXaxis()->CenterTitle(true);
	mg1->GetYaxis()->CenterTitle(true);

	mg1->GetXaxis()->SetLabelSize(textSize);
	mg1->GetXaxis()->SetTitleSize(textSize);
	mg1->GetXaxis()->SetTickLength(0.02);
	mg1->GetXaxis()->SetTitleOffset(1);

	mg1->GetYaxis()->SetLabelSize(textSize);
	mg1->GetYaxis()->SetTitleSize(textSize);
	mg1->GetYaxis()->SetTickLength(0.02);
	mg1->GetYaxis()->SetTitleOffset(1.1);
	mg1->Draw("a");

	TLegend *lg1 = new TLegend(0.2,0.15,0.7,0.45);

	lg1->AddEntry(gr11,"Yield fit, Cleymans et al.","f");
	lg1->AddEntry(gr12,"proton fluc. fit","lp");
	lg1->AddEntry(gr13,"kaon fluc. fit","lp");
	lg1->AddEntry(gr14,"#Lambda fluc. fit","lp");

	lg1->SetTextSize(textSize);
	lg1->SetBorderSize(0);
	lg1->SetFillStyle(0);
	lg1->Draw();

	cv1->SaveAs("plots/freezeoutIndiFit.pdf");
}

int main()
{
	TApplication app {"app", nullptr, nullptr};
	// plotFreezeoutMass();
	// plotFreezeoutCharge();
	plotFreezeoutIndiFit();

	cout << "\nPress Ctrl-C to exit..." << endl;

	app.Run();
	getchar();
	return 0;
}