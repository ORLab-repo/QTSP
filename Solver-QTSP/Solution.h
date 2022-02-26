#pragma once
#include "lib.h"
#include "Node.h"

class Solution {    
public:    
    vector<int> giantT;//giant tour
	vector<Node*> nodes;// for nodes
    Param* pr;    
    int cost;// objective
    int n;// number of customers    
    string typeIns;
    //GA parts:
    vector<int> predecessors; // store prev node of each client  
    vector<int> successors; // store next node of each client  
    multiset<pair<int, Solution*>> indivsPerProximity;
    double biasedFitness; // Biased fitness of the solution
    vector<int> ordNodeLs; // using for LS
    Node* depot;
    Solution(Param* _pr) {
        pr = _pr;
        n = pr->numClient - 1;//not cointain depot
        for (int i = 0; i <= n; ++i) {
            giantT.pb(0);            
        }
        for (int i = 1; i <= n + 3; ++i)nodes.pb(new Node());
        for (int i = 1; i <= n; ++i) {
            nodes[i]->idxClient = i;
            ordNodeLs.push_back(i);
        }
        nodes[n + 1]->idxClient = 0;        
        depot = nodes[n + 1];        
        depot->posInSol = 0;
    }
    void genGiantT() {
        for (int i = 1; i <= n; ++i)giantT[i] = i;
        shuffle(giantT.begin() + 1, giantT.end(), pr->Rng.generator);
        //random_shuffle(giantT.begin() + 1, giantT.end());
    }

    /*diversity contribution in GA*/
    void removeProximity(Solution* indiv)
    {
        auto it = indivsPerProximity.begin();
        while (it->second != indiv) ++it;
        indivsPerProximity.erase(it);
    }


    double brokenPairsDistance(Solution* valSol) {
        int differences = 0;
        for (int j = 1; j <= n; j++)
        {
            if (successors[j] != valSol->successors[j] && successors[j] != valSol->predecessors[j]) differences++;
            if (predecessors[j] == 0 && valSol->predecessors[j] != 0 && valSol->successors[j] != 0) differences++;
        }
        return (double)differences / (double)n;
    }

    double averageBrokenPairsDistanceClosest(int nbClosest) {
        double result = 0;
        int maxSize = min<int>(nbClosest, indivsPerProximity.size());
        auto it = indivsPerProximity.begin();
        for (int i = 0; i < maxSize; i++)
        {
            result += it->first;
            ++it;
        }
        return result / (double)maxSize;
    }
    void updateInfo() {
        Node* valNode = depot;        
        while (valNode != depot->pred)
        {    
            valNode = valNode->suc;
            if (valNode->pred == depot)valNode->culCost = 0;
            else valNode->culCost = valNode->pred->culCost + pr->costs[valNode->pred->pred->idxClient][valNode->pred->idxClient][valNode->idxClient];
            valNode->posInSol = valNode->pred->posInSol + 1;
        }
    }
    /**/
    void calCost() {        
        cost = pr->costs[giantT[n]][giantT[0]][giantT[1]] + pr->costs[giantT[n - 1]][giantT[n]][giantT[0]];
        for (int i = 0; i <= n-2; ++i) {
            cost += pr->costs[giantT[i]][giantT[i + 1]][giantT[i + 2]];
        }
        depot->suc = nodes[giantT[1]];
        depot->pred = nodes[giantT[n]];
        for (int i = 2; i <= n; ++i) {
            int idV = giantT[i], idU = giantT[i - 1];
            nodes[idU]->suc = nodes[idV];
            nodes[idV]->pred = nodes[idU];
        }
        updateInfo();
    }
};