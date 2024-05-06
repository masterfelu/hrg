#ifndef THERMAL_MINIMIZER_H
#define THERMAL_MINIMIZER_H

#include "ThermalFunctions.h"

class ThermalMinimizer
{
	private:

	int funcCount {0};
	double *funcVal {nullptr};
	double *funcVal_err {nullptr};

	double T_init {0.16};
	double mub_init {0.01};
	double muq_init {-0.001};
	double mus_init {0.001};

	double T_final {0};
	double mub_final {0};
	double muq_final {0};
	double mus_final {0};

	double T_err {0};
	double mub_err {0};
	double muq_err {0};
	double mus_err {0};

	double T_err_low {0};
	double mub_err_low {0};
	double muq_err_low {0};
	double mus_err_low {0};

	double T_err_up {0};
	double mub_err_up {0};
	double muq_err_up {0};
	double mus_err_up {0};

	ThermalFunction **tFunc {nullptr};

	bool isMinimized {false};

	public:
	
	ThermalMinimizer(int fc);
	~ThermalMinimizer();
	void setInitialParameters(double T, double mub, double muq, double mus);
	void SetFunctionValueEach(int id, double val, double val_err = 1);
	void SetThermalFunctionEach(int id, ThermalFunction *f);
	double minFunction(const double *x);
	
	void minimize(int maxFC = 1000, double tol = 0.01, bool showLog = true);
	double *getMinimizedParameters();
	double* getMinimizedParameterErrors();
};

double ThermalMinFunction_wrapped(const double *var);
extern ThermalMinimizer *minFuncObjGlobal;

#endif

