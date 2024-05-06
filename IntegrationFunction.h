#ifndef INTEGRAL_FUNCTION_H
#define INTEGRAL_FUNCTION_H

#include "TF1.h"
#include "TF2.h"

static int n_gaussLegendreSamplingPoints = 0;
static double *gaussLegendreWeights = nullptr;
static double *gaussLegendreAbscissas = nullptr;
static bool gaussLegendreSamplingPointsGenerated = false;

void generateGaussLegendreSamplingPoint(int n = 20);

class IntegrationFunction
{
	private:

	double (*function)(double*,double*){ nullptr };

	int n_param { 0 };
	int n_var { 1 };

	TF1 function1D {};
	TF2 function2D {};

	double integralLowerLimit1 { 0 };
	double integralUpperLimit1 { 0 };
	double integralLowerLimit2 { 0 };
	double integralUpperLimit2 { 0 };

	// flags

	enum IntegrationType
	{
		SELF,
		ROOT,
		ROOT_FAST
	};

	IntegrationType integrationType { ROOT };

	bool useIntegralFast { false };
	bool functionInitialized { false };

	public:

	IntegrationFunction() = default;
	IntegrationFunction(int npar, int nvar = 1, bool integralFast = false);

	~IntegrationFunction(){}

	void setLimits(double a1, double b1, double a2 = 0, double b2 = 0);
	void setFunction(double (*f)(double*,double*));
	double integrate(double *param);
};

#endif