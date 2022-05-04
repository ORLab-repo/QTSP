#include "lib.h"
#include "ReadData.h"
//#include "Solution.h"
#include "Ga.h"
#include "ILS.h"

using namespace std;

int arr[] = { 0, 44, 21, 71, 34, 17, 49, 70, 65, 60, 8, 72, 74, 51, 9, 41, 53, 47, 46, 67, 37, 18, 50, 42, 16, 56, 43, 30, 76, 73, 45, 33, 26, 15, 5, 78, 40, 32, 57, 14, 69, 38, 3, 64, 2, 79, 39, 6, 28, 11, 75, 61, 68, 12, 54, 48, 52, 27, 1, 55, 25, 23, 22, 20, 7, 58, 59, 31, 13, 77, 36, 35, 62, 63, 66, 10, 4, 24, 19, 29, };
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
	//Solution bestSol(pr);
	////for (int i = 0; i <= bestSol.n; ++i)bestSol.giantT[i] = arr[i];
	//bestSol.genGiantT();
	//bestSol.calCost();
	//bestSol.exportGiantT();
	////bestSol.cheapestIns();	
	//cout << bestSol.cost << "\n";
	//bestSol.updateObj();
	//cout << bestSol.cost << "\n";
	//cout << "check sol: " << bestSol.calCostWtUpdate() << "\n";		
	//bestSol.worstRmv(5);	
	//bestSol.cheapestIns();
	//cout << bestSol.cost << "\n";
	//cout << bestSol.calCostWtUpdate() << "\n";
	////pertubation...
	//exit(0);
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