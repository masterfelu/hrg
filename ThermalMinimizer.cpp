#include <iostream>
#include <iomanip>
#include <cmath>

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

#include "ThermalMinimizer.h"
#include "ThermalFunctions.h"
#include "Particles.h"

ThermalMinimizer *minFuncObjGlobal;

ThermalMinimizer::ThermalMinimizer(int fc)
{
	isMinimized = false;
	funcCount = fc;
	funcVal = new double[fc];
	funcVal_err = new double[fc];
	tFunc =  new ThermalFunction*[fc];
}

ThermalMinimizer::~ThermalMinimizer()
{
	delete[] funcVal;
	delete[] funcVal_err;
	// delete[] tFunc;
	funcVal = nullptr;
	funcVal_err = nullptr;
	tFunc = nullptr;
}

void ThermalMinimizer::setInitialParameters(double T, double mub, double muq, double mus)
{
	T_init = T;
	mub_init = mub;
	muq_init = muq;
	mus_init = mus;
	isMinimized = false;
}

void ThermalMinimizer::SetFunctionValueEach(int id, double val, double val_err)
{
	if(id < 0 || id >= funcCount){ return; }
	funcVal[id] = val;
	funcVal_err[id] = val_err;
}

void ThermalMinimizer::SetThermalFunctionEach(int id, ThermalFunction *f)
{	
	if(id < 0 || id >= funcCount){ return; }
	// delete tFunc[id];
	tFunc[id] = f;
}

double ThermalMinimizer::minFunction(const double *x)
{
	double T = x[0];
	double mub = x[1];
	double muq = x[2];
	double mus = x[3];

	double sum = 0;

	for (int i = 0; i < funcCount; ++i)
	{
		sum += pow(tFunc[i]->getValue(T,mub,muq,mus)-funcVal[i],2)/pow(funcVal_err[i],2);	// for chi2 test
		// sum += pow(tFunc[i]->getValue(T,mub,muq,mus)-funcVal[i],2);							// for net-pQ paper results
	}
	return sum;					// for chi2 test
	// return sqrt(sum)/funcCount;		// for net-pQ paper results
}

void ThermalMinimizer::minimize(int maxFC, double tol, bool showLog)
{
	ROOT::Minuit2::Minuit2Minimizer minimum {ROOT::Minuit2::kMigrad};
	minimum.SetMaxFunctionCalls(maxFC);
	minimum.SetMaxIterations(maxFC);
	minimum.SetTolerance(tol);
	// minimum.SetErrorDef(1);	// standard dev. for chi2 distr. with n dof

	minFuncObjGlobal = this;
	ROOT::Math::Functor f (&ThermalMinFunction_wrapped,4);
	double step = 1e-2;

	minimum.SetFunction(f);

	minimum.SetVariable(0,"T",T_init, step);
	minimum.SetVariable(1,"mub",mub_init, step);
	minimum.SetVariable(2,"muq",muq_init, step);
	minimum.SetVariable(3,"mus",mus_init, step);

	minimum.SetVariableLimits(0,0.1,0.2);
	minimum.SetVariableLimits(1,-0.5,0.5);
	minimum.SetVariableLimits(2,-0.1,0.1);
	minimum.SetVariableLimits(3,-0.1,0.1);

	// do the minimization
	bool hasConverged = minimum.Minimize();

	const double *x = minimum.X();
	const double *xerr = minimum.Errors();
	
	T_final = x[0];
	mub_final = x[1];
	muq_final = x[2];
	mus_final = x[3];

	// minimum.GetMinosError(0,T_err_low,T_err_up);
	// minimum.GetMinosError(1,mub_err_low,mub_err_up);
	// minimum.GetMinosError(2,muq_err_low,muq_err_up);
	// minimum.GetMinosError(3,mus_err_low,mus_err_up);

	// T_err = (T_err_low+T_err_up)/2;
	// mub_err = (mub_err_low+mub_err_up)/2;
	// muq_err = (muq_err_low+muq_err_up)/2;
	// mus_err = (mus_err_low+mus_err_up)/2;

	T_err = xerr[0];
	mub_err = xerr[1];
	muq_err = xerr[2];
	mus_err = xerr[3];

	if (showLog)
	{
		std::cout << "-----------------------------------------------------------------" << '\n';
		std::cout << "-----------------------------------------------------------------" << '\n';
		std::cout << "        Determination of chemical freezout parameters" << '\n';
		std::cout << "   functions: " << funcCount << " , tol.: " << tol << " , max. calls(approx.): " << maxFC << '\n';
		std::cout << "              Number of function calls(each) : " << minimum.NCalls() << '\n';
		std::cout << "                     Minimization status : " << hasConverged << '\n';
		std::cout << "-----------------------------------------------------------------" << '\n';
		std::cout << " Par. name               Init. Val.    Final. Val.          Error" << '\n';
		std::cout << "-----------------------------------------------------------------" << '\n';
		std::cout << std::left << std::setw(20) << " Temp." << std::right << std::setw(15) << T_init << std::setw(15) << T_final << std::setw(15) << T_err << '\n';
		std::cout << std::left << std::setw(20) << " B chem. pot." << std::right << std::setw(15) << mub_init << std::setw(15) << mub_final << std::setw(15) << mub_err << '\n';
		std::cout << std::left << std::setw(20) << " Q chem. pot." << std::right << std::setw(15) << muq_init << std::setw(15) << muq_final << std::setw(15) << muq_err << '\n';
		std::cout << std::left << std::setw(20) << " S chem. pot." << std::right << std::setw(15) << mus_init << std::setw(15) << mus_final << std::setw(15) << mus_err << '\n';

		std::cout << "-----------------------------------------------------------------" << '\n';
		std::cout << "   #f     Value Obt.     Value Req.    Error prov.     Chi2 Error" << '\n';
		std::cout << "-----------------------------------------------------------------" << '\n';
		double val, chi2err;
		for (int i = 0; i < funcCount; ++i)
		{
			val = tFunc[i]->getValue(x[0],x[1],x[2],x[3]);
			chi2err = pow((val - funcVal[i])/funcVal_err[i],2);
			std::cout << std::setw(5) << i+1 << std::setw(15) << val
			 << std::setw(15) << funcVal[i] << std::setw(15) << funcVal_err[i] << std::setw(15) << chi2err << '\n';
		}
		std::cout << "-----------------------------------------------------------------" << '\n';
		std::cout << std::setw(50) << "Minimiser final value" <<  std::setw(15) << minimum.MinValue() << '\n';
		std::cout << "-----------------------------------------------------------------" << '\n';
		std::cout << "-----------------------------------------------------------------" << '\n';
	}

	isMinimized = true;
	minimum.Clear();
}

double* ThermalMinimizer::getMinimizedParameters()
{
	if(isMinimized)
	{
		double *x = new double[4];
		x[0] = T_final;
		x[1] = mub_final;
		x[2] = muq_final;
		x[3] = mus_final;
		return x;
	}
	return nullptr;
}

double* ThermalMinimizer::getMinimizedParameterErrors()
{
	if(isMinimized)
	{
		double *x = new double[4];
		x[0] = T_err;
		x[1] = mub_err;
		x[2] = muq_err;
		x[3] = mus_err;
		return x;
	}
	return nullptr;
}

double ThermalMinFunction_wrapped(const double *var)
{
	return minFuncObjGlobal->minFunction(var);
}
