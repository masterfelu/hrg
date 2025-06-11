#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>

#include "Particles.h"

SingleParticle::decayChannel::~decayChannel()
{
	delete[] daughterID;
	daughterID = nullptr;
	// std::cout << "decayChannel destructor called..." << '\n';
}

SingleParticle::SingleParticle(const SingleParticle &copy)
	: PDGID { copy.PDGID },
	ID { copy.ID },
	mass { copy.mass },
	spinDegeneracy { copy.spinDegeneracy },
	baryonNumber { copy.baryonNumber },
	charge { copy.charge },
	strangeness { copy.strangeness },
	isStable { copy.isStable },
	isAntiParticle { copy.isAntiParticle },
	n_decayChannels { copy.n_decayChannels },
	name { copy.name }
{
	if(n_decayChannels > 0)
	{
		dChannels = new decayChannel[n_decayChannels];
		for (int i = 0; i < n_decayChannels; ++i)
		{
			dChannels[i].daughters = copy.dChannels[i].daughters;
			dChannels[i].branchingRatio = copy.dChannels[i].branchingRatio;
			dChannels[i].daughterID = new int[dChannels[i].daughters];
			for (int j = 0; j < copy.dChannels[i].daughters; ++j)
			{

				dChannels[i].daughterID[j] = copy.dChannels[i].daughterID[j];
			}
		}
	}
}

SingleParticle& SingleParticle::operator=(const SingleParticle &copy)
{
	if(this != &copy)
	{
		PDGID = copy.PDGID;
		ID = copy.ID;
		mass = copy.mass;
		spinDegeneracy = copy.spinDegeneracy;
		baryonNumber = copy.baryonNumber;
		charge = copy.charge;
		strangeness = copy.strangeness;
		isStable = copy.isStable;
		isAntiParticle = copy.isAntiParticle;
		n_decayChannels = copy.n_decayChannels;
		name = copy.name;
		if(n_decayChannels > 0)
		{
			delete[] dChannels;
			dChannels = new decayChannel[n_decayChannels];
			for (int i = 0; i < n_decayChannels; ++i)
			{
				dChannels[i].daughters = copy.dChannels[i].daughters;
				dChannels[i].branchingRatio = copy.dChannels[i].branchingRatio;
				dChannels[i].daughterID = new int[dChannels[i].daughters];
				for (int j = 0; j < copy.dChannels[i].daughters; ++j)
				{

					dChannels[i].daughterID[j] = copy.dChannels[i].daughterID[j];
				}
			}
		}
	}
	return *this;
}

SingleParticle::~SingleParticle()
{
	delete[] dChannels;
	dChannels = nullptr;
	// std::cout << "SingleParticle destructor called..." << '\n';
}

ParticleSystem::ParticleSystem(const ParticleSystem &copy)
	: temperature { copy.temperature },
	baryonChemicalPotential { copy.baryonChemicalPotential },
	chargeChemicalPotential { copy.chargeChemicalPotential },
	strangenessChemicalPotential { copy.strangenessChemicalPotential },
	species { copy.species }, PDGDataLoaded { copy.PDGDataLoaded },
	decaysLoaded {copy.decaysLoaded }
{
	particles = new SingleParticle[species];
	chemicalPotentialEach = new double[species];
	for (int i = 0; i < species; ++i)
	{
		particles[i] = SingleParticle { copy.particles[i] };
		chemicalPotentialEach[i] = copy.chemicalPotentialEach[i];
	}
}

ParticleSystem& ParticleSystem::operator=(const ParticleSystem &copy)
{
	if(this != &copy)
	{
		temperature = copy.temperature;
		baryonChemicalPotential = copy.baryonChemicalPotential;
		chargeChemicalPotential = copy.chargeChemicalPotential;
		strangenessChemicalPotential = copy.strangenessChemicalPotential;
		species = copy.species;
		PDGDataLoaded = copy.PDGDataLoaded;
		decaysLoaded = copy.decaysLoaded;
		delete[] particles;
		particles = new SingleParticle[species];
		delete[] chemicalPotentialEach;
		chemicalPotentialEach = new double[species];
		for (int i = 0; i < species; ++i)
		{
			particles[i] = SingleParticle { copy.particles[i] };
			chemicalPotentialEach[i] = copy.chemicalPotentialEach[i];
		}
	}
	return *this;
}

ParticleSystem::ParticleSystem(double T, double mub, double muq, double mus)
	: temperature { T },
	baryonChemicalPotential { mub },
	chargeChemicalPotential { muq },
	strangenessChemicalPotential { mus}
{
}

ParticleSystem::~ParticleSystem()
{
	delete[] particles;
	delete[] chemicalPotentialEach;
	particles = nullptr;
	chemicalPotentialEach = nullptr;
	// std::cout << "ParticleSystem destructor called..." << '\n';
}

void ParticleSystem::addParticle(SingleParticle &p)
{
	SingleParticle *temp = new SingleParticle[species+1];
	double *tempMus = new double[species+1];
	for (int i = 0; i < species; ++i)
	{
		temp[i] = SingleParticle { particles[i] };
	}
	temp[species] = SingleParticle { p };
	delete[] particles;
	delete[] chemicalPotentialEach;
	particles = temp;
	chemicalPotentialEach = tempMus;
	species++;
	calculateChemicalPotentialEach();
}

SingleParticle &ParticleSystem::getParticle(int ID)
{
	if(ID < 0 || ID >= species)
	{
		SingleParticle* empty = new SingleParticle{}; 
	 	return *empty;
	}
	return particles[ID];
}

int ParticleSystem::getIDfromPDGID(long PDGID)
{
	for (int i = 0; i < species; ++i)
	{
		if(particles[i].PDGID == PDGID)
		{
			return particles[i].ID;
		}
	}
	return -1;
}

int ParticleSystem::getAntiParticleID(int ID)
{	
	if(ID < 0 || ID >= species){	return -1;	}
	int antiID = getIDfromPDGID(-particles[ID].PDGID);
	if(antiID != -1)
	{
		return antiID;
	}
	else{
		return ID;	// Particle is its own anti-particle
	}
}

void ParticleSystem::loadPDGTable(std::string fname)
{
	std::ifstream inFile;
	inFile.open(fname);
	long PDGID;			//PDG ID
	std::string name;		//Name
	float m;			//mass(in GeV)
	short int g;		//spin degeneracy

	short int B;		//Baryon number
	short int S;		//Strangeness
	short int Q;		//Charge
	bool isStable;		//Stable or not
	bool isAnti;		//Anti-particle or not

	if(inFile){
		inFile >> species;
		particles = new SingleParticle[species];
		for(int i = 0; i < species ; i++){
			inFile >> PDGID >> name >> isStable >> m >> g >> B >> Q >> S;
			if (PDGID < 0)	{	isAnti = true;	}
			else 			{	isAnti = false;	}
			particles[i] = SingleParticle { i, PDGID, name, m, g, B, Q, S, isStable, isAnti};
		}
	}
	calculateChemicalPotentialEach();
	PDGDataLoaded = true;
	inFile.close();
}

void ParticleSystem::calculateChemicalPotentialEach()
{
	delete[] chemicalPotentialEach;
	chemicalPotentialEach = new double[species];
	for (int i = 0; i < species; ++i)
	{
		chemicalPotentialEach[i] = particles[i].baryonNumber*baryonChemicalPotential;
		chemicalPotentialEach[i] += particles[i].charge*chargeChemicalPotential;
		chemicalPotentialEach[i] += particles[i].strangeness*strangenessChemicalPotential;
	}
}

void ParticleSystem::setTemperature(double T)
{
	temperature = T;
}

void ParticleSystem::setBaryonChemicalPotential(double mub)
{
	baryonChemicalPotential = mub;
	calculateChemicalPotentialEach();
}

void ParticleSystem::setChargeChemicalPotential(double muq)
{
	chargeChemicalPotential = muq;
	calculateChemicalPotentialEach();
}

void ParticleSystem::settStrangenessChemicalPotential(double mus)
{
	strangenessChemicalPotential = mus;
	calculateChemicalPotentialEach();
}

void ParticleSystem::setChemicalPotential(double mub, double muq, double mus)
{
	baryonChemicalPotential = mub;
	chargeChemicalPotential = muq;
	strangenessChemicalPotential = mus;
	calculateChemicalPotentialEach();
}

double ParticleSystem::getChemicalPotentialEach(int ID)
{
	if(ID < 0 || ID >= species){ return 0; }
	return chemicalPotentialEach[ID];
}

void ParticleSystem::setThermalParameters(double T, double mub, double muq, double mus)
{
	temperature = T;
	baryonChemicalPotential = mub;
	chargeChemicalPotential = muq;
	strangenessChemicalPotential = mus;
	calculateChemicalPotentialEach();
}

void ParticleSystem::loadDecayChannels(std::string fname){
	std::ifstream inFile;
	inFile.open(fname);
	if(inFile){
		long PDGID;

		int dSpecies, ID,  daughterID, n_decayChannels, daughters;

		int antiFlag;					//Raised if antiparticle exists
		int antiID, antidaughterID;		

		float BR;
		inFile >> dSpecies;
		for (int i = 0; i < dSpecies; ++i)
		{
			inFile >> PDGID;
			ID = getIDfromPDGID(PDGID);
			if(ID != -1){
				antiID = getAntiParticleID(ID);	
				if(antiID != ID){
					antiFlag = 1;
				}
				else{
					antiFlag = 0;
				}
				inFile >> n_decayChannels;
				particles[ID].n_decayChannels = n_decayChannels;
				particles[ID].dChannels = new SingleParticle::decayChannel[n_decayChannels];

				if(antiFlag){
					particles[antiID].n_decayChannels = n_decayChannels;
					particles[antiID].dChannels = new SingleParticle::decayChannel[n_decayChannels];
				}
				for (int j = 0; j < n_decayChannels; ++j)
				{	
					//Storing the daughterIDs for each channel of the decayed hadron
					inFile >> BR >> daughters;
					particles[ID].dChannels[j].branchingRatio = BR;
					particles[ID].dChannels[j].daughters = daughters;
					particles[ID].dChannels[j].daughterID = new int[daughters];

					if (antiFlag)
					{
						particles[antiID].dChannels[j].branchingRatio = BR;
						particles[antiID].dChannels[j].daughters = daughters;
						particles[antiID].dChannels[j].daughterID = new int[daughters];

					}
					for (int k = 0; k < daughters; ++k)
					{

						inFile >> PDGID;
						daughterID = getIDfromPDGID(PDGID);
						particles[ID].dChannels[j].daughterID[k] = daughterID;
						if (antiFlag)
						{
							antidaughterID = getAntiParticleID(daughterID);
							particles[antiID].dChannels[j].daughterID[k] = antidaughterID;
						}	
					}
				}
			}
			else{
				inFile >> n_decayChannels;
				for (int j = 0; j < n_decayChannels; ++j)
				{	
					inFile >> BR >> daughters;
					for (int k = 0; k < daughters; ++k)
					{
						inFile >> PDGID;
					}
				}
			}
		}
	}
	decaysLoaded = true;
	inFile.close();
}

void ParticleSystem::loadDefaultData()
{
	loadPDGTable();
	loadDecayChannels();
}

void ParticleSystem::printDataToFile(std::string fname)
{
	std::ofstream outFile;
	outFile.open(fname);
	if(outFile)
	{
		if (!PDGDataLoaded)
		{	
			outFile << "No PDG Data Loaded..." << '\n';
			return;
		}
		outFile << '\n';
		outFile << "                                    #######################################################                                            " << '\n';
		outFile << "                                    #                        PDG data                     #                                            " << '\n';
		outFile << "                                    #      Default File name: data/pdg_hadrons.dat        #                                            " << '\n';
		outFile << "                                    #######################################################                                            " << '\n';
		outFile << '\n';
		outFile << '\n';
		outFile << "---------------------------------------------------------------------------------------------------------------------------------------" << '\n';
		outFile << "ID     PDGID          Name                isStable  mass(in GeV)  Spin Degeneracy     Baryon number(B)  Charge(Q)  Strangeness(S)      " << '\n';
		outFile << "---------------------------------------------------------------------------------------------------------------------------------------" << '\n';
		for (int i = 0; i < species; ++i)
		{
			outFile << std::left << std::setw(7) << particles[i].ID;
			outFile << std::setw(15) << particles[i].PDGID; 
			outFile << std::setw(20) << particles[i].name;
			outFile << std::setw(10) << particles[i].isStable;
			outFile << std::setw(14) << particles[i].mass; 
			outFile << std::setw(20) << particles[i].spinDegeneracy;
			outFile << std::setw(18) << particles[i].baryonNumber;
			outFile << std::setw(11) << particles[i].charge;
			outFile << std::setw(14) << particles[i].strangeness << '\n';
		}
		outFile << "---------------------------------------------------------------------------------------------------------------------------------------" << '\n';
		if(!decaysLoaded)
		{
			outFile << "No Decays Loaded..." << '\n';
			return;
		}
		outFile << '\n';
		outFile << '\n';
		outFile << "                                    #######################################################                                            " << '\n';
		outFile << "                                    #                     Decays data                     #                                            " << '\n';
		outFile << "                                    #          Default File name: data/decays.dat         #                                            " << '\n';
		outFile << "                                    #-----------------------------------------------------#                                            " << '\n';
		outFile << "                                    #         In the table if -1 is given instead         #                                            " << '\n';
		outFile << "                                    #      of name, it means the PDG ID has not been      #                                            " << '\n';
		outFile << "                                    #                 found in the table.                 #                                            " << '\n';
		outFile << "                                    #######################################################                                            " << '\n';
		outFile << '\n';
		outFile << "---------------------------------------------------------------------------------------------------------------------------------------" << '\n';
		outFile << "ParentHadron        BranchingRatio(%)   Daughter1      Daughter2      Daughter3      Daughter4      Daughter5                          " << '\n';
		outFile << "---------------------------------------------------------------------------------------------------------------------------------------" << '\n';
		for (int i = 0; i < species; ++i)
		{
			if (particles[i].n_decayChannels > 0)
			{
				outFile << std::left << std::setw(20) << particles[i].name << '\n';
				for (int j = 0; j < particles[i].n_decayChannels; ++j)
				{
					outFile << std::left << std::setw(20) << " " <<std::setw(20) << particles[i].dChannels[j].branchingRatio*100;
					for (int k = 0; k < particles[i].dChannels[j].daughters; ++k)
					{
						if (particles[i].dChannels[j].daughterID[k] < 0)
						{
							outFile << std::setw(15) << particles[i].dChannels[j].daughterID[k];
						}
						else{
							outFile << std::setw(15) << particles[particles[i].dChannels[j].daughterID[k]].name;
						}
					}
					outFile << '\n';
				}
				outFile << '\n';
			}
		}
		outFile << "---------------------------------------------------------------------------------------------------------------------------------------" << '\n';
		outFile << "------------------------------------------------------------xxxxxxxxxxxxxxxx-----------------------------------------------------------" << '\n';
	}
	outFile.close();
}