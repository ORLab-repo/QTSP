#pragma once
#include "RandomGenerator.h"
#include "Location.h"

class Location;
class RandomGenerator;

class Param {
public:
	std::string nameIns;
	ofstream fileOut;
	ofstream fl;
	int Q;// capacity
	int T;// time horizon	
	int numVeh = int(1e9);// number of vehicles
	int maxVeh;// max number of vehicles (used for Moccia 2012 instances).
	int numClient;//number of client (containing depot)	
	int numLoc;// number of node	
	vector<Location> listLoc;
	vector<vector<int> > corDis;// correlation measure	
	//vector<vector<vector<double> >> costs;// distance	
	double costs[204][204][204];
	double fSup[204][204];
	double fMain[204][204];
	int predSup[204][204];
	int predMain[204][204];
	vector<vector<int> > new_costs;// distance	
	vector<vector<int> > times;// distance		
	vector<vector<vector<int> > > correlatedNodes;
	bool isDebug = false;
	bool isTurnCkSol = false;
	int lambda = 2;
	clock_t start, end;
	double total = 0;
	int nPop = 25;
	int delta = 40;
	int nMut = 10;
	double totalIntra = 0;
	double rateMut = 1.1;
	double rate4Opt = 0.3;
	int initItSCP = 2000;	
	bool use4Opt = false;
	//int maxNeibor = 20;// max size of neigbor vertex set for granular search (can be used for dynamically sertification).
	int maxNeibor = 20;// max size of neigbor
	//int maxNeibor = 30;// max size of neigbor
	int worstDeg = 3;
	int ldTw = 5;// coff for TW vertex set for granular search (can be used for dynamically sertification).
	int nbR = 20;
	int nbIls = 100;
	int nbF = 20;
	int sizeSub = 10;// max size of subsequence for using concatenation
	int TL = 1800;// time limit (default for 120-cus ins)
	bool bi = true;//cheking best improvement
	int nI = 10;
	int nC = 10;
	int pMin = 1;
	int pMax = 2;
	//LNS-based params:
	int maxRmv = 40;
	int minRmv = 5;	
	double minRateSmallRmv = 0.01;
	double maxRateSmallRmv = 0.1;
	double minRateBigRmv = 0.1;
	double maxRateBigRmv = 0.3;
	int debugLS = false;
	bool m2Opt = false;
	bool crosver = false;
	bool swap2 = false;
	bool moveor = false;
	bool stringex = false;
	
	string pathOut;
	RandomGenerator Rng;
	~Param() {
		cout << "clear param\n";		
		listLoc.clear();
		/*for (int i = 0; i < costs.size(); ++i)costs[i].clear();
		costs.clear();		*/
	};
};
