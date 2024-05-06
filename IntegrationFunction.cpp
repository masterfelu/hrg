#include "TF1.h"
#include "TF2.h"
#include "Math/Integrator.h"

#include "IntegrationFunction.h"

IntegrationFunction::IntegrationFunction(int npar, int nvar, bool integralFast)
{
	n_param = npar;
	n_var = nvar; 
	useIntegralFast = integralFast;
}

void generateGaussLegendreSamplingPoint(int n)
{
	TF1 f{};
	n_gaussLegendreSamplingPoints = n;
	gaussLegendreWeights = new double[n];
	gaussLegendreAbscissas = new double[n];
	f.CalcGaussLegendreSamplingPoints(n,gaussLegendreAbscissas,gaussLegendreWeights);
	gaussLegendreSamplingPointsGenerated = true;
}

void IntegrationFunction::setLimits(double a1, double b1, double a2, double b2)
{
	integralLowerLimit1 = a1;
	integralUpperLimit1 = b1;
	integralLowerLimit2 = a2;
	integralUpperLimit2 = b2;
}

void IntegrationFunction::setFunction(double (*f)(double*,double*))
{
	function = f;
	if (n_var == 1)
	{
		ROOT::Math::IntegratorOneDimOptions::SetDefaultIntegrator("Gauss");
		function1D = {"1D_Integral", f, integralLowerLimit1, integralUpperLimit1, n_param};
		functionInitialized = true;
	}
	else if (n_var == 2)
	{
		function2D = {"2D_Integral", f, integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2, n_param};
		functionInitialized = true;
	}
}

double IntegrationFunction::integrate(double *param){
	if(!functionInitialized){
		return 0;
	}
	if (n_var == 1)
	{
		function1D.SetParameters(param);
		if(useIntegralFast)
		{
			if(!gaussLegendreSamplingPointsGenerated){	generateGaussLegendreSamplingPoint();	}
			return function1D.IntegralFast(n_gaussLegendreSamplingPoints, gaussLegendreAbscissas, gaussLegendreWeights, integralLowerLimit1, integralUpperLimit1, param);
		}
		return function1D.Integral(integralLowerLimit1, integralUpperLimit1);
	}
	else if (n_var == 2)
	{
		function2D.SetParameters(param);
		return function2D.Integral(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
	}
	return 0;
}