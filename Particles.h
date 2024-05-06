#ifndef PARTICLES_H
#define PARTICLES_H

#include <string>
#include <iostream>

////////////////////////////////////////////
// This class contains all the attibutes  //
// necessary to define a single           //
// hadron and information about its       //
// decay.                                 //
////////////////////////////////////////////

class SingleParticle
{

	public:

	// Particle attributes

	//Basic attributes

	long PDGID { 0 };					// PDG ID
	int ID { -1 };						// Index of particle (starts from 0)
	std::string name { };				// Name
	float mass { 0 };					// mass(in GeV)
	short int spinDegeneracy { 0 };		// spin degeneracy

	// Conserved charges

	short int baryonNumber { 0 };		// Baryon number	
	short int charge { 0 };				// Charge
	short int strangeness { 0 };		// Strangeness

	// Extra attributes

	float radius { 0 };					// Hadron radius (in GeV^-1)

	// Flags

	bool isStable { false };			// Stable or not
	bool isAntiParticle { false };		// Whether it is anti-particle or not

	// Decay attributes(useful if not stable)

	int n_decayChannels { 0 };			// No. of decay Channels
	struct decayChannel
	{
		int daughters { 0 };			// No. of daughters in the channel
		int *daughterID { nullptr };	// IDs of each daughters in the channel
		float branchingRatio { 0 };		// Branching ratio of the channel	
		~decayChannel();
	}; 

	decayChannel *dChannels { nullptr };

	// default constructor

	SingleParticle(){}

	SingleParticle(int i_ID, long i_PDGID, std::string i_name, float i_m,
	 short int i_spinDegeneracy, short int i_B, short int i_Q, short int i_S,
	 bool i_isStable, bool i_isAntiParticle)
		: PDGID { i_PDGID },
		ID { i_ID },
		name { i_name },
		mass { i_m },
		spinDegeneracy { i_spinDegeneracy },
		baryonNumber { i_B },
		charge { i_Q },
		strangeness { i_S},
		isStable { i_isStable},
		isAntiParticle { i_isAntiParticle }
	{
	}

	// copy constructor

	SingleParticle(const SingleParticle &copy);
	SingleParticle& operator=(const SingleParticle &copy);

	// destructor

	~SingleParticle();

};

class ParticleSystem
{
	private:

	//particle attributes

	int species { 0 };
	SingleParticle *particles { nullptr };

	//thermal attributes of the system

	double temperature { 0 };					// in GeV units
	double baryonChemicalPotential { 0 };		// in GeV units
	double chargeChemicalPotential { 0 };		// in GeV units
	double strangenessChemicalPotential { 0 };	// in GeV units
	double *chemicalPotentialEach { nullptr };	// chemical potential of each species

	// flags

	bool PDGDataLoaded { false };
	bool decaysLoaded { false };

	public:

	ParticleSystem(){}

	// copy constructor

	ParticleSystem(const ParticleSystem &copy);
	ParticleSystem& operator=(const ParticleSystem &copy);

	ParticleSystem(double T, double mub, double muq, double mus);

	// destructor

	~ParticleSystem();

	int getSpecies() const{	return species;	}

	double getTemperature () const 					{	return temperature;	}
	double getBaryonChemicalPotential () const		{	return baryonChemicalPotential;	}
	double getChargeChemicalPotential () const		{	return chargeChemicalPotential;	}
	double getStrangenessChemicalPotential () const	{	return strangenessChemicalPotential;	}
	const double *getChemicalPotentialAll () const	{	return chemicalPotentialEach;			}
	double getChemicalPotentialEach (int ID);

	void setTemperature(double T);
	void setBaryonChemicalPotential(double mub);
	void setChargeChemicalPotential(double muq);
	void settStrangenessChemicalPotential(double mus);
	void setChemicalPotential(double mub, double muq, double mus);
	void setThermalParameters(double T, double mub, double muq, double mus);
	void calculateChemicalPotentialEach();

	int getIDfromPDGID(long PDGID);
	int getAntiParticleID(int ID);

	void addParticle(SingleParticle &p);
	SingleParticle &getParticle(int ID);

	void loadPDGTable(std::string fname = "data/pdg_hadrons.dat");
	void loadDecayChannels(std::string fname = "data/decays.dat");
	void loadDefaultData();
	void printDataToFile(std::string fname = "PDG_and_decays.txt");
};

#endif