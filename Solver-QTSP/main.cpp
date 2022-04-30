#include "lib.h"
#include "ReadData.h"
//#include "Solution.h"
#include "Ga.h"
#include "ILS.h"

using namespace std;

int arr[] = { 0, 50, 29, 65, 54, 66, 47, 76, 43, 40, 31, 69, 60, 19, 11, 30, 9, 41, 53, 28, 51, 62, 72, 67, 37, 42, 32, 35, 58, 10, 2, 48, 38, 8, 14, 1, 23, 22, 20, 26, 4, 77, 25, 55, 24, 70, 33, 7, 63, 78, 56, 21, 44, 13, 57, 18, 46, 6, 73, 5, 71, 61, 16, 52, 49, 79, 45, 74, 27, 36, 75, 68, 12, 3, 15, 64, 39, 17, 34, 59, };
int seed[] = {
	18319894,
	23390422,
	36197069,
	45945346,
	54500951,
	63196367,
	71110057,
	89578146,
	96527670,
	10415237,
};
string nameIns = "PointSet_80_3.tsp";
string type = "ag";
int main(int argc, char* argv[]) {
	for (int i = 1; i < argc; ++i) {		
		if (string(argv[i]) == "-nameIns") {
			nameIns = argv[i + 1];
		}
		if (string(argv[i]) == "-type") {
			type = argv[i + 1];

		}
	}
	string pathIn = "PointSets\\" + nameIns;
	string pathOut = "solution\\" + nameIns + "_" + type + ".sol";
	Param* pr = read_Ins(pathIn, type);		
	pr->fileOut.open(pathOut);
	cout << setprecision(5) << fixed;		
	//pr->isDebug = true;
	Solution bestSol(pr);
	for (int i = 0; i <= bestSol.n; ++i)bestSol.giantT[i] = arr[i];
	bestSol.calCost();
	//bestSol.cheapestIns();	
	cout << bestSol.cost << "\n";
	bestSol.updateObj();
	cout << bestSol.cost << "\n";
	//bestSol.worstRmv(5);	
	//bestSol.cheapestIns();
	//cout << bestSol.cost << "\n";
	//cout << bestSol.calCostWtUpdate() << "\n";
	////pertubation...
	exit(0);
	ILS ilsAlgo;
	ilsAlgo.init(pr);
	ilsAlgo.runAlgo();
	exit(0);
	//double bestCost = oo;		
	//for (int i = 1; i <= 1000; ++i) {		
	//	cout << i << "\n";
	//	try {
	//		bestSol.genGiantT();
	//		bestSol.calCost();
	//		bestSol.updateObj();
	//		bestCost = min(bestCost, bestSol.cost);
	//	}
	// 
	//	catch (const char* msg) {
	//		cerr << msg << endl;
	//		exit(0);
	//		//system("pause");
	//	}
	//	catch (...) {
	//		cout << "error\n";
	//		exit(0);
	//		//system("pause");
	//	}		
	//}
	//cout << "nbMove1: " << bestSol.nbMove1 << "\n";
	//cout << "nbMove2: " << bestSol.nbMove2 << "\n";
	//cout << "nbMove3: " << bestSol.nbMove3 << "\n";
	//cout << "nbMove4: " << bestSol.nbMove4 << "\n";
	//cout << "nbMove5: " << bestSol.nbMove5 << "\n";
	//cout << "nbMove6: " << bestSol.nbMove6 << "\n";
	//cout << "nbMove7: " << bestSol.nbMove7 << "\n";
	//cout << bestCost << "\n";
	GA Algo;
	Algo.init(pr);
	cout << "finish init\n";
	double minCost = oo;
	double sumCost = 0;
	for (int numRun = 0; numRun < 10; ++numRun) {
		pr->Rng.config(seed[numRun]);
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		try {
			Algo.findGasSol();
		}
		catch (const char* msg) {
			cerr << msg << endl;
			exit(0);
			//system("pause");
		}
		catch (...) {
			cout << "error\n";
			exit(0);
			//system("pause");
		}
		minCost = min(minCost, Algo.bestCost);
		sumCost += Algo.bestCost;
		cout << "name ins: " << nameIns << " " << numRun << " " << Algo.bestCost << "\n";
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::time_t end_time = std::chrono::system_clock::to_time_t(end);		
		pr->fileOut << "elapsed time: " << elapsed_seconds.count() << "s\n\n";
	}
	pr->fileOut << fixed << setprecision(2) << "best run: " << minCost << "\n";
	pr->fileOut << fixed << setprecision(2) << "avg run: " << (double)sumCost / 10 << "\n";
	pr->fileOut.close();
	//system("pause");
	return 0;	
}