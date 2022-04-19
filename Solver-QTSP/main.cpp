#include "lib.h"
#include "ReadData.h"
//#include "Solution.h"
#include "Ga.h"

using namespace std;

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
int main() {			
	Param* pr = read_Ins("PointSets\\PointSet_50_6.tsp", "ag");	
	pr->Rng.config(seed[0]);
	cout << setprecision(5) << fixed;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	/*pr->isDebug = true;*/
	//Solution bestSol(pr);
	//bestSol.cheapestIns();	
	////pertubation...
	//exit(0);
	//double bestCost = oo;		
	//for (int i = 1; i <= 1000; ++i) {		
	//	cout << i << "\n";
	//	try {
	//		bestSol.genGiantT();
	//		bestSol.calCost();
	//		bestSol.updateObj();
	//		bestCost = min(bestCost, bestSol.cost);
	//	}
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
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	cout << "elapsed time: " << elapsed_seconds.count() << "s\n\n";
	system("pause");
	/*return 0;	*/
}