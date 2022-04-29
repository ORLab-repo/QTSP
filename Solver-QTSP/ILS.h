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
	}	
	void equalSol(Solution* u, Solution* v);
	void runAlgo();
private:
};
