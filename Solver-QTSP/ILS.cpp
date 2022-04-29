#include "ILS.h"

ILS::ILS()
{
}

ILS::~ILS()
{
}

void ILS::equalSol(Solution* u, Solution* v)
{
    u->cost = v->cost;    
    for (int i = 1; i <= n; ++i) {
        u->giantT[i] = v->giantT[i];        
    }
}

void ILS::runAlgo()
{
	curSol->cheapestIns();
	curSol->updateObj();
	equalSol(oriSol, curSol);
	equalSol(bestSol, curSol);
	while (true)
	{
		equalSol(curSol, oriSol);
		curSol->calCost();
		curSol->pertubation();
		curSol->updateObj();
		if (curSol->cost < bestSol->cost) {
			equalSol(bestSol, curSol);
			equalSol(oriSol, curSol);
		}
		else if ((curSol->cost - oriSol->cost) / oriSol->cost < acceptRate) {
			equalSol(oriSol, curSol);
		}
		else {
			equalSol(oriSol, bestSol);
		}
		acceptRate *= coolingRate;
	}
}
