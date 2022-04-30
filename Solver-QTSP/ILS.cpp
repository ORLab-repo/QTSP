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
	//curSol->cheapestIns();
	curSol->genGiantT();	
	curSol->genGiantT();
	curSol->exportGiantT();
	curSol->calCost();
	cout << "ori cost: " << curSol->cost << "\n";
	curSol->updateObj();
	cout << "improved cost: " << curSol->cost << "\n";		
	cout << "check sol: " << boolalpha << curSol->checkSol() << "\n";	
	cout << curSol->calCostWtUpdate() << "\n";	
	pr->fileOut.close();
	equalSol(oriSol, curSol);
	equalSol(bestSol, curSol);
	cout << "initial cost: " << bestSol->cost << "\n";
	for (int i = 1; i <= totalIt; ++i)
	{
		//cout << i << "\n";
		equalSol(curSol, oriSol);
		curSol->calCost();
		curSol->pertubation();
		curSol->updateObj();
		if (bestSol->cost - curSol->cost > MY_EPSILON) {						
			equalSol(bestSol, curSol);
			cout << "new best: " << bestSol->cost << "\n";
			equalSol(oriSol, curSol);
		}
		else if ((curSol->cost - oriSol->cost)  -  oriSol->cost*acceptRate > EP) {
			equalSol(oriSol, curSol);
		}
		else {
			equalSol(oriSol, bestSol);
		}
		acceptRate *= coolingRate;
	}
}
