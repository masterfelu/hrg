#include <cmath>

#include "ThermalFunctions.h"
#include "Particles.h"
#include "IntegrationFunction.h"
#include "MiscellaneousFunctions.h"

const double PI = 3.1415926539;

// ------------------------------------------ThermalFunction------------------------------------------ //

ThermalFunction::ThermalFunction(const ParticleSystem &sys)
{
	system  = ParticleSystem { sys };
	species  = sys.getSpecies();
	valueEach = new double[species];
}

ThermalFunction::ThermalFunction(const ThermalFunction &copy):
ThermalFunction { copy.system }
{
	value = copy.value ;
	valueEach = new double[species];
	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = copy.valueEach[i];
	}
	integralLowerLimit1 = copy.integralLowerLimit1;
	integralUpperLimit1 = copy.integralUpperLimit1;
	integralLowerLimit2 = copy.integralLowerLimit2;
	integralUpperLimit2 = copy.integralUpperLimit2;
	integralDimension = copy.integralDimension;
	isQuantityCalculated = copy.isQuantityCalculated;
	isSingleParticleSystem = copy.isSingleParticleSystem;
}

ThermalFunction::~ThermalFunction()
{
	delete[] valueEach;
	valueEach = nullptr;
	// std::cout << "ThermalFunction destructor called..." << '\n';
}

void ThermalFunction::setIntegralLimits(double a1, double b1, double a2, double b2)
{
	integralLowerLimit1 = a1;
	integralUpperLimit1 = b1;

	integralLowerLimit2 = a2;
	integralUpperLimit2 = b2;
}

void ThermalFunction::setSphericalUniformCoordinates()
{
	integralDimension = 1;
	coordinateSystem = SPHERICAL_UNIFORM;
}

void ThermalFunction::setDetectorRapidityCoordinates()
{
	integralDimension = 2;
	coordinateSystem = DETECTOR_RAPIDITY;
}

void ThermalFunction::setDetectorPseudorapidityCoordinates()
{
	integralDimension = 2;
	coordinateSystem = DETECTOR_PSEUDORAPIDITY;
}

double ThermalFunction::getValue()
{
	if (isQuantityCalculated)
	{
		return value;
	}
	return 0;
}

double ThermalFunction::getValueEach(int idx)
{
	if (isQuantityCalculated)
	{
		return valueEach[idx];
	}
	return 0;
}

// ------------------------------------------Pressure------------------------------------------ //

Pressure::Pressure(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);

	}
}

Pressure::~Pressure()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
}

double Pressure::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};

	double val,E,pVol,func;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;
	func = sign*g*T*log(1+sign*exp(-(E-mu)/T));

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Pressure::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double y = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];

	double val,E,pVol,mT,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);
	func = sign*g*T*log(1+sign*exp(-(E-mu)/T));

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Pressure::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];

	double val,E,pVol,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);
	func = sign*g*T*log(1+sign*exp(-(E-mu)/T));

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Pressure::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[4];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 4, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = T;
			param[1] = system.getParticle(i).baryonNumber*mub;
			param[1] += system.getParticle(i).charge*muq;
			param[1] += system.getParticle(i).strangeness*mus;
			param[2] = system.getParticle(i).mass;
			param[3] = system.getParticle(i).spinDegeneracy;
			value += integral.integrate(param);
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = T;
			param[1] = system.getParticle(i).baryonNumber*mub;
			param[1] += system.getParticle(i).charge*muq;
			param[1] += system.getParticle(i).strangeness*mus;
			param[2] = system.getParticle(i).mass;
			param[3] = system.getParticle(i).spinDegeneracy;
			value += integral.integrate(param);
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = T;
			param[1] = system.getParticle(i).baryonNumber*mub;
			param[1] += system.getParticle(i).charge*muq;
			param[1] += system.getParticle(i).strangeness*mus;
			param[2] = system.getParticle(i).mass;
			param[3] = system.getParticle(i).spinDegeneracy;
			value += integral.integrate(param);
		}
	}
	return value;
}

void Pressure::calculateValue()
{
	value = 0;
	double param[4];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 4, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility1------------------------------------------ //

Susceptibility1::Susceptibility1(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);

	}
}

Susceptibility1::~Susceptibility1()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
}

double Susceptibility1::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = n/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility1::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double y = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = n/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility1::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = n/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility1::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[4];
	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 4, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility1::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility1::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility1::calculateValue()
{
	value = 0;
	double param[4];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 4, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{				
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility1B------------------------------------------ //

Susceptibility1B::Susceptibility1B(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	baryonNumber = new int[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);
		baryonNumber[i] = system.getParticle(i).baryonNumber;

	}
}

Susceptibility1B::~Susceptibility1B()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	delete[] baryonNumber;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
	baryonNumber = nullptr;
}

double Susceptibility1B::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double B {par[4]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = B*n/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility1B::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double y = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double B {par[4]};

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = B*n/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility1B::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double B {par[4]};

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = B*n/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility1B::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[5];
	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.baryonNumber;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.baryonNumber;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.baryonNumber;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility1B::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility1B::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility1B::calculateValue()
{
	value = 0;
	double param[5];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = baryonNumber[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = baryonNumber[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = baryonNumber[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility1Q------------------------------------------ //

Susceptibility1Q::Susceptibility1Q(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	charge = new int[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);
		charge[i] = system.getParticle(i).charge;

	}
}

Susceptibility1Q::~Susceptibility1Q()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	delete[] charge;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
	charge = nullptr;
}

double Susceptibility1Q::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double Q {par[4]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = Q*n/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility1Q::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double y = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double Q {par[4]};

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = Q*n/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility1Q::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double Q {par[4]};

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = Q*n/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility1Q::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[5];

	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.charge;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.charge;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.charge;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility1Q::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility1Q::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility1Q::calculateValue()
{
	value = 0;
	double param[5];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = charge[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = charge[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = charge[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility1S------------------------------------------ //

Susceptibility1S::Susceptibility1S(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	strangeness = new int[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);
		strangeness[i] = system.getParticle(i).strangeness;

	}
}

Susceptibility1S::~Susceptibility1S()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	delete[] strangeness;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
	strangeness = nullptr;
}

double Susceptibility1S::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double S {par[4]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = S*n/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility1S::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double y = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double S {par[4]};

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = S*n/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility1S::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double S {par[4]};

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = S*n/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility1S::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[5];

	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.strangeness;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.strangeness;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.strangeness;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility1S::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility1S::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility1S::calculateValue()
{
	value = 0;
	double param[5];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = strangeness[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = strangeness[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = strangeness[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility2------------------------------------------ //

Susceptibility2::Susceptibility2(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);

	}
}

Susceptibility2::~Susceptibility2()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
}

double Susceptibility2::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = n*(1-sign*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility2::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double y = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = n*(1-sign*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility2::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = n*(1-sign*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility2::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[4];
	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 4, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility2::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility2::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility2::calculateValue()
{
	value = 0;
	double param[4];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 4, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{				
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility2B------------------------------------------ //

Susceptibility2B::Susceptibility2B(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	baryonNumber = new int[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);
		baryonNumber[i] = system.getParticle(i).baryonNumber;

	}
}

Susceptibility2B::~Susceptibility2B()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	delete[] baryonNumber;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
}

double Susceptibility2B::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double B {par[4]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = B*B*n*(1-sign*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility2B::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double y = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double B {par[4]};

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = B*B*n*(1-sign*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility2B::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double B {par[4]};

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = B*B*n*(1-sign*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility2B::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[5];

	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.baryonNumber;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.baryonNumber;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.baryonNumber;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility2B::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility2B::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility2B::calculateValue()
{
	value = 0;
	double param[5];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = baryonNumber[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = baryonNumber[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = baryonNumber[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility2Q------------------------------------------ //

Susceptibility2Q::Susceptibility2Q(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	charge = new int[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);
		charge[i] = system.getParticle(i).charge;

	}
}

Susceptibility2Q::~Susceptibility2Q()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	delete[] charge;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
	charge = nullptr;
}

double Susceptibility2Q::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double Q {par[4]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = Q*Q*n*(1-sign*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility2Q::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double y = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double Q {par[4]};

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = Q*Q*n*(1-sign*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility2Q::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double Q {par[4]};

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = Q*Q*n*(1-sign*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility2Q::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[5];

	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.charge;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.charge;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.charge;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility2Q::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility2Q::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility2Q::calculateValue()
{
	value = 0;
	double param[5];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = charge[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = charge[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = charge[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility2S------------------------------------------ //

Susceptibility2S::Susceptibility2S(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	strangeness = new int[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);
		strangeness[i] = system.getParticle(i).strangeness;

	}
}

Susceptibility2S::~Susceptibility2S()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	delete[] strangeness;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
	strangeness = nullptr;
}

double Susceptibility2S::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double S {par[4]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = S*S*n*(1-sign*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility2S::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double y = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double S {par[4]};

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = S*S*n*(1-sign*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility2S::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double S {par[4]};

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = S*S*n*(1-sign*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility2S::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[5];

	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.strangeness;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.strangeness;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.strangeness;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility2S::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility2S::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility2S::calculateValue()
{
	value = 0;
	double param[5];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = strangeness[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = strangeness[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = strangeness[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility3------------------------------------------ //

Susceptibility3::Susceptibility3(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);

	}
}

Susceptibility3::~Susceptibility3()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
}

double Susceptibility3::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility3::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double y = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility3::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility3::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[4];
	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 4, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility3::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility3::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility3::calculateValue()
{
	value = 0;
	double param[4];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 4, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{				
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility3B------------------------------------------ //

Susceptibility3B::Susceptibility3B(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	baryonNumber = new int[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);
		baryonNumber[i] = system.getParticle(i).baryonNumber;

	}
}

Susceptibility3B::~Susceptibility3B()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	delete[] baryonNumber;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
}

double Susceptibility3B::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double B {par[4]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = B*B*B*n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility3B::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double y = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double B {par[4]};

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = B*B*B*n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility3B::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double B {par[4]};

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = B*B*B*n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility3B::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[5];

	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.baryonNumber;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.baryonNumber;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.baryonNumber;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility3B::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility3B::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility3B::calculateValue()
{
	value = 0;
	double param[5];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = baryonNumber[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = baryonNumber[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = baryonNumber[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility3Q------------------------------------------ //

Susceptibility3Q::Susceptibility3Q(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	charge = new int[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);
		charge[i] = system.getParticle(i).charge;

	}
}

Susceptibility3Q::~Susceptibility3Q()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	delete[] charge;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
	charge = nullptr;
}

double Susceptibility3Q::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double Q {par[4]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = Q*Q*Q*n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility3Q::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double y = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double Q {par[4]};

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = Q*Q*Q*n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility3Q::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];
	double Q {par[4]};

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = Q*Q*Q*n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility3Q::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[5];

	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.charge;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.charge;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.charge;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility3Q::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility3Q::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility3Q::calculateValue()
{
	value = 0;
	double param[5];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = charge[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = charge[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = charge[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility3S------------------------------------------ //

Susceptibility3S::Susceptibility3S(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	strangeness = new int[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);
		strangeness[i] = system.getParticle(i).strangeness;

	}
}

Susceptibility3S::~Susceptibility3S()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	delete[] strangeness;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
	strangeness = nullptr;
}

double Susceptibility3S::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double S {par[4]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = S*S*S*n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility3S::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT {var[0]};
	double y {var[1]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double S {par[4]};

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = S*S*S*n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility3S::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double S {par[4]};

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = S*S*S*n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility3S::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[5];

	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.strangeness;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.strangeness;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.strangeness;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility3S::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility3S::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility3S::calculateValue()
{
	value = 0;
	double param[5];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = strangeness[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = strangeness[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = strangeness[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility4------------------------------------------ //

Susceptibility4::Susceptibility4(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);

	}
}

Susceptibility4::~Susceptibility4()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
}

double Susceptibility4::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = n*(1-7*sign*n+12*n*n-6*sign*n*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility4::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double y = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = n*(1-7*sign*n+12*n*n-6*sign*n*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility4::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT = var[0];
	double eta = var[1];

	double T = par[0];
	double mu = par[1];
	double m = par[2];
	double g = par[3];

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = n*(1-7*sign*n+12*n*n-6*sign*n*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility4::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[4];
	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 4, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility4::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility4::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility4::calculateValue()
{
	value = 0;
	double param[4];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 4, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 4, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{				
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility4B------------------------------------------ //

Susceptibility4B::Susceptibility4B(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	baryonNumber = new int[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);
		baryonNumber[i] = system.getParticle(i).baryonNumber;

	}
}

Susceptibility4B::~Susceptibility4B()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	delete[] baryonNumber;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
}

double Susceptibility4B::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double B {par[4]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = B*B*B*B*n*(1-7*sign*n+12*n*n-6*sign*n*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility4B::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT {var[0]};
	double y {var[1]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double B {par[4]};

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = B*B*B*B*n*(1-7*sign*n+12*n*n-6*sign*n*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility4B::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT {var[0]};
	double eta {var[1]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double B {par[4]};

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = B*B*B*B*n*(1-7*sign*n+12*n*n-6*sign*n*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility4B::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[5];

	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.baryonNumber;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.baryonNumber;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.baryonNumber;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility4B::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility4B::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility4B::calculateValue()
{
	value = 0;
	double param[5];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = baryonNumber[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = baryonNumber[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = baryonNumber[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility4Q------------------------------------------ //

Susceptibility4Q::Susceptibility4Q(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	charge = new int[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);
		charge[i] = system.getParticle(i).charge;

	}
}

Susceptibility4Q::~Susceptibility4Q()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	delete[] charge;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
	charge = nullptr;
}

double Susceptibility4Q::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double Q {par[4]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = Q*Q*Q*n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility4Q::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT {var[0]};
	double y {var[1]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double Q {par[4]};

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = Q*Q*Q*Q*n*(1-7*sign*n+12*n*n-6*sign*n*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility4Q::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT {var[0]};
	double eta {var[1]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double Q {par[4]};

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = Q*Q*Q*Q*n*(1-7*sign*n+12*n*n-6*sign*n*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility4Q::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[5];

	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.charge;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.charge;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.charge;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility4Q::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility4Q::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility4Q::calculateValue()
{
	value = 0;
	double param[5];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = charge[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = charge[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = charge[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}

// ------------------------------------------Susceptibility4S------------------------------------------ //

Susceptibility4S::Susceptibility4S(ParticleSystem &sys)
	: ThermalFunction {sys}
{
	temperature = system.getTemperature();
	mass = new float[species];
	chemicalPotential = new double[species];
	spinDegeneracy = new double[species];
	strangeness = new int[species];
	for (int i = 0; i < species; ++i)
	{
		mass[i] = system.getParticle(i).mass;
		spinDegeneracy[i] = system.getParticle(i).spinDegeneracy;
		chemicalPotential[i] = system.getChemicalPotentialEach(i);
		strangeness[i] = system.getParticle(i).strangeness;

	}
}

Susceptibility4S::~Susceptibility4S()
{
	delete[] mass;
	delete[] chemicalPotential;
	delete[] spinDegeneracy;
	delete[] strangeness;
	mass = nullptr;
	chemicalPotential = nullptr;
	spinDegeneracy = nullptr;
	strangeness = nullptr;
}

double Susceptibility4S::integrandSphericalUniformCoordinates(double *var, double *par)
{
	double p {var[0]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double S {par[4]};

	double val,E,pVol,func,n;
	int sign;

	E = sqrt(m*m+p*p);
	sign = pow(-1,g);
	pVol = 4*PI*p*p;

	n = g/(exp((E-mu)/T)+sign);
	func = S*S*S*n*(1-3*sign*n+2*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility4S::integrandDetectorRapidityCoordinates(double *var, double *par)
{
	double pT {var[0]};
	double y {var[1]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double S {par[4]};

	double val,E,pVol,mT,n,func;
	int sign;

	mT = sqrt(m*m+pT*pT);
	E = mT*cosh(y);
	sign = pow(-1,g);
	pVol = 2*PI*pT*mT*cosh(y);

	n = g/(exp((E-mu)/T)+sign);
	func = S*S*S*S*n*(1-7*sign*n+12*n*n-6*sign*n*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility4S::integrandDetectorPseudorapidityCoordinates(double *var, double *par)
{
	double pT {var[0]};
	double eta {var[1]};

	double T {par[0]};
	double mu {par[1]};
	double m {par[2]};
	double g {par[3]};
	double S {par[4]};

	double val,E,pVol,n,func;
	int sign;

	E = sqrt(pT*cosh(eta)*pT*cosh(eta)+m*m);
	sign = pow(-1,g);
	pVol = 2*PI*pT*pT*cosh(eta);


	n = g/(exp((E-mu)/T)+sign);
	func = S*S*S*S*n*(1-7*sign*n+12*n*n-6*sign*n*n*n)/(T*T*T);

	val = pVol*func/(8*PI*PI*PI);
	return val;
}

double Susceptibility4S::getValueEachWithoutDecay(int idx, double T, double mub, double muq, double mus)
{
	double value = 0;
	double param[5];

	SingleParticle currParticle = system.getParticle(idx);

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.strangeness;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.strangeness;
		value = integral.integrate(param);
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		param[0] = T;
		param[1] = currParticle.baryonNumber*mub;
		param[1] += currParticle.charge*muq;
		param[1] += currParticle.strangeness*mus;
		param[2] = currParticle.mass;
		param[3] = currParticle.spinDegeneracy;
		param[4] = currParticle.strangeness;
		value = integral.integrate(param);
	}
	return value;
}

double Susceptibility4S::getValue(double T, double mub, double muq, double mus)
{
	double value = 0;
	double valueEach[species];

	for (int i = 0; i < species; ++i)
	{
		valueEach[i] = getValueEachWithoutDecay(i,T,mub,muq,mus);
		value += valueEach[i];
	}

	if(enableDecay)
	{
		double newValue = 0;
		double newValueEach[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

double Susceptibility4S::getValueEach(int idx, double T, double mub, double muq, double mus)
{
	double value = getValueEachWithoutDecay(idx,T,mub,muq,mus);
	if(enableDecay)
	{
		double newValue = value;
		double valueEach = 0;

		bool decayFlag {false};

		SingleParticle stableParticle, decayParticle;
		stableParticle = system.getParticle(idx);
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		if (stableParticle.isStable)
		{
			for (int i = 0; i < species; ++i)
			{
				decayParticle = system.getParticle(i);
				if (!decayParticle.isStable)
				{
					if(decayParticle.n_decayChannels > maxL){
						// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
						maxL = decayParticle.n_decayChannels;

						nr = new int[maxL];
						br = new float[maxL];
					}
					decayFlag = false;
					for (int j = 0; j < decayParticle.n_decayChannels; ++j)
					{
						br[j] = decayParticle.dChannels[j].branchingRatio;
						nr[j] = 0;
						for (int k = 0; k < decayParticle.dChannels[j].daughters; ++k)
						{
							if (idx == decayParticle.dChannels[j].daughterID[k])
							{
								nr[j] += 1;
								decayFlag = true;
							}
						}
					}
					nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
					// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
					// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
					// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);
					valueEach = 0;
					if (decayFlag)
					{
						valueEach = getValueEachWithoutDecay(i,T,mub,muq,mus);
					}
					newValue += valueEach*nr1;
				}
			}
		}
		delete[] nr;
		delete[] br;
		value = newValue;
	}
	return value;
}

void Susceptibility4S::calculateValue()
{
	value = 0;
	double param[5];

	if(coordinateSystem == SPHERICAL_UNIFORM)
	{
		IntegrationFunction integral { 5, 1, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1);
		integral.setFunction(integrandSphericalUniformCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = strangeness[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_RAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorRapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = strangeness[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	else if(coordinateSystem == DETECTOR_PSEUDORAPIDITY)
	{
		IntegrationFunction integral { 5, 2, useFastIntegral };
		integral.setLimits(integralLowerLimit1, integralUpperLimit1, integralLowerLimit2, integralUpperLimit2);
		integral.setFunction(integrandDetectorPseudorapidityCoordinates);
		for (int i = 0; i < species; ++i)
		{
			param[0] = temperature;
			param[1] = chemicalPotential[i];
			param[2] = mass[i];
			param[3] = spinDegeneracy[i];
			param[4] = strangeness[i];
			valueEach[i] = integral.integrate(param);
			value += valueEach[i];
		}
	}
	if(enableDecay)
	{
		double newValue = 0;
		double *newValueEach = new double[species];

		SingleParticle stableParticle, decayParticle;
		int maxL = 60;

		int *nr = new int[maxL];
		float *br = new float[maxL];

		double nr1, dnr2, dnr3, dnr4;
		for (int i = 0; i < species; ++i)
		{
			newValueEach[i] = valueEach[i];
			stableParticle = system.getParticle(i);
			if (stableParticle.isStable)
			{
				for (int j = 0; j < species; ++j)
				{
					decayParticle = system.getParticle(j);
					if (!decayParticle.isStable)
					{
						if(decayParticle.n_decayChannels > maxL){
							// cout << "Need bigger array of size " << particles[i]->decayChannels << endl;
							maxL = decayParticle.n_decayChannels;

							nr = new int[maxL];
							br = new float[maxL];
						}
						for (int k = 0; k < decayParticle.n_decayChannels; ++k)
						{
							br[k] = decayParticle.dChannels[k].branchingRatio;
							nr[k] = 0;
							for (int l = 0; l < decayParticle.dChannels[k].daughters; ++l)
							{
								if (i == decayParticle.dChannels[k].daughterID[l])
								{
									nr[k] += 1;
								}
							}
						}
						nr1 = moment(nr,br,decayParticle.n_decayChannels,1);
						// dnr2 = centralMoment(nr,br,particles[i]->decayChannels,2);
						// dnr3 = centralMoment(nr,br,particles[i]->decayChannels,3);
						// dnr4 = centralMoment(nr,br,particles[i]->decayChannels,4);

						newValueEach[i] += valueEach[j]*nr1;
					}
				}
				newValue += newValueEach[i];
			}
		}
		delete[] valueEach;
		valueEach = newValueEach;
		value = newValue;
		delete[] nr;
		delete[] br;
	}
	isQuantityCalculated = true;
}