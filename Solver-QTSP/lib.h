#pragma once

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <stdio.h>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <queue>
#include <cstdlib>   
#include <chrono>
#include <ctime>
#include "Param.h"
#include "Util.h"
#include "Node.h"
//#include <ilcplex/ilocplex.h>

#define pb push_back
#define mp make_pair
#define ft first
#define sc second
#define EP 1e-6
#define oo 1e9
#define MY_EPSILON 0.00001 // Precision parameter, used to avoid numerical instabilities

//using namespace std;

typedef pair<int, int> II;
typedef pair<II, int> III;
typedef pair<double, double> DD;
typedef pair<double, int> DI;
typedef pair<double, II> DII;
//typedef IloArray<IloNumVarArray> IloNumVarArray2;//2D
//typedef IloArray<IloNumVarArray2> IloNumVarArray3;//3D
//typedef IloArray<IloNumVarArray3> IloNumVarArray4;//4D