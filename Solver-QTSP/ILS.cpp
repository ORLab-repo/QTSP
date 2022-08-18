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
	/*curSol->genGiantTOri();			
	curSol->calCost();*/
	curSol->cheapestIns();
	//curSol->randomIns();
	//curSol->nearestIns();
	cout << "ori cost: " << curSol->cost << "\n";
	cout << curSol->calCostWtUpdate() << "\n";
	cout << boolalpha << curSol->checkSol() << "\n";
	curSol->updateObj();
	/*curSol->exportGiantT();	
	pr->fileOut.close();*/
	cout << "improved cost: " << curSol->cost << "\n";		
	cout << "check sol: " << boolalpha << curSol->checkSol() << "\n";	
	cout << curSol->calCostWtUpdate() << "\n";	
	//curSol->exportGiantT();	
	equalSol(oriSol, curSol);
	equalSol(bestSol, curSol);
	cout << "initial cost: " << bestSol->cost << "\n";
	int iter = 0;
	for (int i = 1; i <= totalIt; ++i)
	{
		//cout << i << "\n";
		equalSol(curSol, oriSol);
		curSol->calCost();
		curSol->pertubation(true);
	/*	if (iter < omega) {
			equalSol(curSol, oriSol);
			curSol->calCost();
			curSol->pertubation(false);
			iter++;
		}
		else {
			equalSol(curSol, bestSol);
			curSol->calCost();
			curSol->pertubation(true);
			iter = 0;
		}*/
		curSol->updateObj();
		//double oriCost = curSol->cost;
		/*if ((1 + lsRate) * bestSol->cost - curSol->cost > MY_EPSILON) {
			curSol->updateObj();
		}*/
		if (i % 100 == 0) {
			cout << i << " " << "new best: " << bestSol->cost << "\n";
		}		
		/*if (oriSol->cost - curSol->cost > MY_EPSILON) {
			equalSol(oriSol, curSol);
		}	
		if (bestSol->cost - curSol->cost > MY_EPSILON) {						
			equalSol(bestSol, curSol);			
			iter = 0;
		}*/
		if (bestSol->cost - curSol->cost > MY_EPSILON) {						
			equalSol(bestSol, curSol);
			//cout <<i<<" "<< "new best: " << bestSol->cost << "\n";
			equalSol(oriSol, curSol);
		}
		else if ((curSol->cost - oriSol->cost)  -  oriSol->cost*acceptRate > MY_EPSILON) {
			equalSol(oriSol, bestSol);
		}
		else {
			equalSol(oriSol, curSol);
		}
		acceptRate *= coolingRate;
	}

	
}

void ILS::RandR(Solution* baseSol)
{
	// TODO: insert return statement here
	equalSol(curSol, baseSol);
	equalSol(oriSol, curSol);
	equalSol(bestSol, curSol);
	for (int i = 1; i <= totalIt; ++i)
	{
		//cout << i << "\n";
		equalSol(curSol, oriSol);
		curSol->calCost();
		curSol->pertubation(true);

		if (bestSol->cost - curSol->cost > MY_EPSILON) {
			equalSol(bestSol, curSol);
			//cout <<i<<" "<< "new best: " << bestSol->cost << "\n";
			equalSol(oriSol, curSol);
		}
		else if ((curSol->cost - oriSol->cost) - oriSol->cost * acceptRate > MY_EPSILON) {
			equalSol(oriSol, bestSol);
		}
		else {
			equalSol(oriSol, curSol);
		}
		acceptRate *= coolingRate;
	}

}
