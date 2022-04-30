#pragma once
#include "Solution.h"
#include <ctime>

class ILS
{
public:
	int n;	
	Param* pr;
	double acceptRate = 0.1;
	double coolingRate = 0.9995;
	int step = 10000;
	int totalIt = 100000;
	Solution* bestSol;
	Solution* oriSol;
	Solution* curSol;		
	ILS();
	~ILS();
	void init(Param* _pr) {
		pr = _pr;
		n = pr->numLoc - 1;
		bestSol = new Solution(pr);
		oriSol = new Solution(pr);
		curSol = new Solution(pr);
		coolingRate = pow(coolingRate, (double)step / totalIt);
	}	
	void equalSol(Solution* u, Solution* v);
	void runAlgo();
private:
};
