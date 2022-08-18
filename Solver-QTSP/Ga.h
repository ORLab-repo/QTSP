#pragma once
#include "Solution.h"
#include <ctime>

class GA
{
public:
	int n;
	double bestCost;
	/*const int nClose = 2;
	const int nElite = 1;	*/
	const int nClose = 3;
	const int nElite = 10;	
	/*int nPop = 15;
	int delta = 25;*/
	int nPop = 20;
	int delta = 40;
	/*const int ItSCP = 1000;
	const int ItNI = 3000;
	const int nMut = 10;*/	
	const int ItNI = 9000;
	const int totalIT = 50000;
	//const int ItNI = 30000;
	int nMut = 10;
	//const int ItMut = 500;
	double pMinMut = 0.3;
	const double pMaxMut = 0.6;
	//const int delta = 83;
	vector<Solution*> pop;// [nPop + 200];
	Solution* valPop;
	vector<double> adapt;// [nPop + 200] ;
	double sumAdapt;
	vector<II> sortAdapt;// [nPop + 200] ;	
	/// 	
	int curNPop;
	int numNotCha;
	int ddID[1000];
	int omega = 0;
	int threshold;	
	double pM = 1.1;
	Param* pr;
	GA();
	~GA();
	void init(Param* _pr) {
		pr = _pr;
		n = pr->numLoc - 1;
		nPop = pr->nPop;
		delta = pr->delta;
		pM = pr->rateMut;
		nMut = pr->nMut;		
		//pMinMut = pr->rateMut;
		//init for SCP		
		for (int i = 0; i <= nPop + 2*delta; ++i) {
		//for (int i = 0; i <= 6; ++i) {			
			pop.push_back(new Solution(_pr));
			adapt.push_back(0);
			sortAdapt.push_back(II(0, 0));
		}
		Solution* val = new Solution(_pr);
		valPop = val;
	}	
	///
	// Evaluates the biased fitness of all individuals in the population
	void updateBiasedFitnesses();
	void removeWorstIndv();
	///

    bool checkIdSol(Solution* u);
	int broken_pairs(Solution* u, Solution* v);
    bool CheckEqual(Solution* u, Solution* v);
    void equalSol(Solution* u, Solution* v);    
    void FindAdapt();    
    int getChild();    
    void choose(int& u, int& v);    
    void uni(Solution* u, Solution* v, Solution* u1, Solution* v1, int numga);
    void insertNew(Solution* u);
    void DelPopu();
	void InitPopu(bool isEdu);	
    void DiversifyPopu(Solution* bestSol);        
	void findGasSol(int maxNumGas = 5000);

private:
};

