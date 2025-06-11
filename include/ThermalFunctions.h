#ifndef THERMAL_FUNCTIONS_H
#define THERMAL_FUNCTIONS_H

#include "Particles.h"

class ThermalFunction
{
	protected:

	ParticleSystem system {};

	int species { 0 };				// number of species in system
	double value { 0 };				// value of the function
	double *valueEach { nullptr };	// value for each particle (if needed)

	//integral limits

	double integralLowerLimit1 { 0 };	// Variable 1
	double integralUpperLimit1 { 0 };	// generally momentum (single integral) or transverse momentum (double integral)

	double integralLowerLimit2 { 0 };	// variable 2 (in case of double integral)
	double integralUpperLimit2 { 0 };	// generally rapidity/pseudorapidity

	short int integralDimension { 1 }; 	//1D integral (1) or 2D integral (2)

	// flags

	bool isSingleParticleSystem { false };	// whether the system only contains a single particle
	bool isQuantityCalculated { false };	// whether the necessary calculations are done or not
	bool useFastIntegral { false };			// whether use fast integral or not

	enum CoordinateSystem
	{	
		SPHERICAL_UNIFORM,
		DETECTOR_RAPIDITY,
		DETECTOR_PSEUDORAPIDITY
	};
	CoordinateSystem coordinateSystem { SPHERICAL_UNIFORM };

	public:

	ThermalFunction() = default;
	ThermalFunction(const ThermalFunction &copy);
	ThermalFunction& operator=(const ThermalFunction &copy) = delete;
	ThermalFunction(const ParticleSystem &sys);
	~ThermalFunction();

	void setIntegralLimits(double a1, double b1, double a2 = 0, double b2 = 0);

	ParticleSystem &getParticleSystem(){ return system; }
	
	void setSphericalUniformCoordinates();
	void setDetectorRapidityCoordinates();
	void setDetectorPseudorapidityCoordinates();

	bool isDetectorCoordinates()			{	return coordinateSystem != SPHERICAL_UNIFORM;		}
	bool isRapidityCoordinates()			{	return coordinateSystem == DETECTOR_RAPIDITY;		}
	bool isPseudorapidityCoordinates()		{	return coordinateSystem == DETECTOR_PSEUDORAPIDITY;	}

	virtual void calculateValue(){ isQuantityCalculated = true; }
	double getValue();
	double getValueEach(int index);
	virtual double getValue(double T, double mub, double muq, double mus){ return 0; }
	virtual double getValueEach(int idx, double T, double mub, double muq, double mus){ return 0; }
};

class Pressure : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	

	Pressure() = default;
	Pressure(ParticleSystem &sys);
	~Pressure();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility1 : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	

	bool enableDecay { true };

	Susceptibility1() = default;
	Susceptibility1(ParticleSystem &sys);
	~Susceptibility1();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility1B : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	
	int *baryonNumber { nullptr };

	bool enableDecay { true };

	Susceptibility1B() = default;
	Susceptibility1B(ParticleSystem &sys);
	~Susceptibility1B();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility1Q : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	
	int *charge { nullptr };

	bool enableDecay { true };

	Susceptibility1Q() = default;
	Susceptibility1Q(ParticleSystem &sys);
	~Susceptibility1Q();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility1S : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	
	int *strangeness { nullptr };

	bool enableDecay { true };

	Susceptibility1S() = default;
	Susceptibility1S(ParticleSystem &sys);
	~Susceptibility1S();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility2 : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	

	bool enableDecay { true };

	Susceptibility2() = default;
	Susceptibility2(ParticleSystem &sys);
	~Susceptibility2();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility2B : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	
	int *baryonNumber { nullptr };

	bool enableDecay { true };

	Susceptibility2B() = default;
	Susceptibility2B(ParticleSystem &sys);
	~Susceptibility2B();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility2Q : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	
	int *charge { nullptr };

	bool enableDecay { true };

	Susceptibility2Q() = default;
	Susceptibility2Q(ParticleSystem &sys);
	~Susceptibility2Q();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility2S : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	
	int *strangeness { nullptr };

	bool enableDecay { true };

	Susceptibility2S() = default;
	Susceptibility2S(ParticleSystem &sys);
	~Susceptibility2S();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility3 : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	

	bool enableDecay { true };

	Susceptibility3() = default;
	Susceptibility3(ParticleSystem &sys);
	~Susceptibility3();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility3B : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	
	int *baryonNumber { nullptr };

	bool enableDecay { true };

	Susceptibility3B() = default;
	Susceptibility3B(ParticleSystem &sys);
	~Susceptibility3B();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility3Q : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	
	int *charge { nullptr };

	bool enableDecay { true };

	Susceptibility3Q() = default;
	Susceptibility3Q(ParticleSystem &sys);
	~Susceptibility3Q();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility3S : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	
	int *strangeness { nullptr };

	bool enableDecay { true };

	Susceptibility3S() = default;
	Susceptibility3S(ParticleSystem &sys);
	~Susceptibility3S();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility4 : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	

	bool enableDecay { true };

	Susceptibility4() = default;
	Susceptibility4(ParticleSystem &sys);
	~Susceptibility4();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility4B : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	
	int *baryonNumber { nullptr };

	bool enableDecay { true };

	Susceptibility4B() = default;
	Susceptibility4B(ParticleSystem &sys);
	~Susceptibility4B();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility4Q : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	
	int *charge { nullptr };

	bool enableDecay { true };

	Susceptibility4Q() = default;
	Susceptibility4Q(ParticleSystem &sys);
	~Susceptibility4Q();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

class Susceptibility4S : public ThermalFunction
{
	private:

	public:

	double temperature { 0 };
	double *chemicalPotential { nullptr };
	float *mass { nullptr };
	double *spinDegeneracy { nullptr };	
	int *strangeness { nullptr };

	bool enableDecay { true };

	Susceptibility4S() = default;
	Susceptibility4S(ParticleSystem &sys);
	~Susceptibility4S();

	using ThermalFunction::getValue;
	using ThermalFunction::getValueEach;

	void calculateValue();
	double getValue(double T, double mub, double muq, double mus);
	double getValueEach(int idx, double T, double mub, double muq, double mus);
	double getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus);
	void setDecay(bool decayFlag) { enableDecay = decayFlag; }

	static double integrandSphericalUniformCoordinates(double *var, double *par);
	static double integrandDetectorRapidityCoordinates(double *var, double *par);
	static double integrandDetectorPseudorapidityCoordinates(double *var, double *par);
};

#endif