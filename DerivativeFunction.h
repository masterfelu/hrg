#ifndef DERIVATIVE_FUNCTION_H
#define DERIVATIVE_FUNCTION_H

#include "TF1.h"

class DerivativeFunction
{
	private:

	double (*function)(double,double*){ nullptr };

	int n_param { 0 };

	TF1 function1D {};

	// flags

	enum DerivativeType
	{
		SELF,
		ROOT
	};

	DerivativeType derivativeType { ROOT };

	bool functionInitialized { false };

	public:

	DerivativeFunction() = default;
	DerivativeFunction(int npar, int nvar = 1, bool integralFast = false);

	~DerivativeFunction(){}

	void setLimits(double a1, double b1, double a2 = 0, double b2 = 0);
	void setFunction(double (*f)(double*,double*));
	double integrate(double *param);
};

#endif