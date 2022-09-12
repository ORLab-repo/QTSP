#include "lib.h"
#include "ReadData.h"
//#include "Solution.h"
#include "Ga.h"
#include "ILS.h"

using namespace std;

int arr[] = { 0, 11, 15, 29, 10, 1, 5, 28, 32, 56, 68, 34, 18, 58, 53, 31, 21, 39, 70, 7, 44, 33, 65, 69, 42, 12, 62, 3, 25, 37, 49, 35, 9, 46, 40, 41, 23, 73, 8, 63, 2, 16, 30, 20, 27, 24, 61, 50, 22, 47, 60, 36, 54, 43, 45, 71, 72, 17, 66, 6, 13, 64, 38, 19, 74, 67, 55, 51, 14, 4, 26, 52, 48, 57, 59, };
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
//string nameIns = "PointSet_75_4.tsp";
string nameIns = "PointSet_150_3.tsp";
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
	////for (int i = 0; i <= bestSol.n; ++i)bestSol.giantT[i] = arr[i];
	//bestSol.genGiantT();	
	//bestSol.genGiantT();
	//bestSol.calCost();
	//cout << bestSol.cost << "\n";
	//int i1, i2, i3, i4;
	//cout << bestSol.fastDoubleBridge(i1, i2, i3, i4) << "\n";	
	//bestSol.apply4Opt(i1, i2, i3, i4);
	//cout << bestSol.cost << "\n";
	/*bestSol.doubleBridge(false, i1, i2, i3, i4);
	//cout << i1 << " " << i2 << " " << i3 << " " << i4 << "\n";
	//bestSol.doubleBridge(true, i1, i2, i3, i4);
	//cout << bestSol.cost << "\n";	*/
	//exit(0);
	//bestSol.updateObj();
	////bestSol.oldCalCost();
	////////bestSol.exportGiantT();
	//////////bestSol.cheapestIns();	
	//cout << boolalpha << bestSol.checkSol() << "\n";		
	//cout << bestSol.cost << " " << bestSol.calCostWtUpdate() << "\n";
	//bestSol.exportGiantT();
	//exit(0);
	//bestSol.updateObj();
	//cout << bestSol.cost << " " << bestSol.calCostWtUpdate() << "\n";	
	////bestSol.exportGiantT();
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