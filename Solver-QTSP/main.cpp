#include "lib.h"
#include "ReadData.h"
//#include "Solution.h"
#include "Ga.h"
#include "ILS.h"

using namespace std;

int arr[] = { 0, 15, 49, 137, 6, 63, 71, 141, 61, 27, 84, 1, 9, 55, 80, 16, 32, 40, 142, 69, 105, 112, 140, 36, 62, 118, 138, 68, 3, 23, 94, 119, 132, 117, 37, 44, 116, 57, 99, 120, 109, 13, 147, 115, 89, 73, 111, 24, 101, 22, 45, 102, 54, 107, 74, 7, 64, 91, 149, 75, 90, 17, 86, 59, 65, 82, 29, 81, 128, 108, 10, 20, 26, 104, 110, 114, 70, 95, 144, 146, 52, 18, 58, 135, 97, 50, 51, 47, 25, 85, 2, 100, 134, 5, 35, 14, 125, 123, 72, 83, 103, 126, 124, 131, 43, 60, 96, 133, 145, 136, 11, 106, 12, 78, 130, 139, 31, 67, 93, 53, 77, 79, 4, 30, 56, 28, 38, 121, 48, 39, 92, 46, 98, 41, 129, 88, 21, 143, 113, 8, 127, 19, 33, 122, 148, 76, 66, 42, 34, 87, };
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
string nameIns = "PointSet_150_3.tsp";
//string nameIns = "PointSet_50_1.tsp";
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
	pr->Rng.config(seed[0]);
	cout << setprecision(5) << fixed;		
	//pr->isDebug = true;
	//Solution bestSol(pr);	
	//for (int i = 0; i <= bestSol.n; ++i)bestSol.giantT[i] = arr[i];
	//bestSol.genGiantT();
	//bestSol.calCost();
	//bestSol.oldCalCost();
	//////bestSol.exportGiantT();
	////////bestSol.cheapestIns();	
	//cout << boolalpha << bestSol.checkSol() << "\n";
	//cout << bestSol.cost << " " << bestSol.calCostWtUpdate() << "\n";
	//bestSol.updateObj();
	//bestSol.exportGiantT();
	//cout << bestSol.cost << "\n";
	//cout << "check sol: " << bestSol.calCostWtUpdate() << "\n";		
	//bestSol.exportGiantT();
	////bestSol.worstRmv(5);	
	////bestSol.cheapestIns();
	////cout << bestSol.cost << "\n";
	////cout << bestSol.calCostWtUpdate() << "\n";
	//////pertubation...
	//cout << "nbMove1: " << bestSol.nbMove1 << "\n";
	//cout << "nbMove2: " << bestSol.nbMove2 << "\n";
	//cout << "nbMove3: " << bestSol.nbMove3 << "\n";
	//cout << "nbMove4: " << bestSol.nbMove4 << "\n";
	//cout << "nbMove5: " << bestSol.nbMove5 << "\n";
	//cout << "nbMove6: " << bestSol.nbMove6 << "\n";
	//cout << "nbMove7: " << bestSol.nbMove7 << "\n";
	//exit(0);
	/*ILS ilsAlgo;
	ilsAlgo.init(pr);
	ilsAlgo.runAlgo();
	exit(0);*/
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
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::time_t end_time = std::chrono::system_clock::to_time_t(end);		
		pr->fileOut << "elapsed time: " << elapsed_seconds.count() << "s\n\n";
		cout << "name ins: " << nameIns << " " << numRun << " " << Algo.bestCost <<" "<< elapsed_seconds.count() << "\n";
	}
	pr->fileOut << fixed << setprecision(2) << "best run: " << minCost << "\n";
	pr->fileOut << fixed << setprecision(2) << "avg run: " << (double)sumCost / 10 << "\n";
	pr->fileOut.close();
	//system("pause");
	return 0;	
}