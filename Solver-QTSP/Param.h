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
	vector<vector<vector<double> >> costs;// distance	
	vector<vector<int> > new_costs;// distance	
	vector<vector<int> > times;// distance		
	vector<vector<int>> correlatedNodes;
	bool isDebug = false;
	bool isTurnCkSol = false;
	int lambda = 2;
	clock_t start, end;
	double total = 0;
	int nPop = 40;
	int delta = 80;
	int nMut = 10;
	double totalIntra = 0;
	double rateMut = 0.2;
	int initItSCP = 2000;
	//int maxNeibor = 40;// max size of neigbor vertex set for granular search (can be used for dynamically sertification).
	int maxNeibor = 120;// max size of nei
	int ldTw = 5;// coff for TWgbor vertex set for granular search (can be used for dynamically sertification).
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
	int debugLS = false;
	bool m2Opt = false;
	bool crosver = false;
	bool swap2 = false;
	bool moveor = false;
	bool stringex = false;
	string pathOut;
	RandomGenerator Rng;
	Param() {
		cout << "clear param\n";		
		listLoc.clear();
		for (int i = 0; i < costs.size(); ++i)costs[i].clear();
		costs.clear();		
	};
};
