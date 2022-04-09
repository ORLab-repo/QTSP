#include "lib.h"
#include "ReadData.h"
#include "Solution.h"

using namespace std;

int main() {		
	Param* pr = read_Ins("PointSets\\PointSet_25_1.tsp", "ag");	
	Solution bestSol(pr);
	bestSol.genGiantT();
	bestSol.calCost();
	cout << bestSol.cost << "\n";
	bestSol.updateObj();
	cout << bestSol.cost << "\n";
	return 0;
}