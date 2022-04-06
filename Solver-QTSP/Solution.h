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
    vector<SeqData*> mySeq;
    SeqData* seqSet;
    SeqData* seqDep;
    SeqData* valSeq[7];    
    //GA parts:
    vector<int> predecessors; // store prev node of each client  
    vector<int> successors; // store next node of each client  
    multiset<pair<int, Solution*>> indivsPerProximity;
    double biasedFitness; // Biased fitness of the solution
    vector<int> ordNodeLs; // using for LS
    Node* depot;
    Node* nodeU;
    Node* uPred;
    Node* uSuc;
    Node* nodeV;
    Node* vPred;
    Node* vSuc;
    int nodeUIdx, nodeVIdx, predUIdx, predVIdx, sucUIdx, sucVIdx;
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
        int maxSizeSeq = 4 * n + 1;
        SeqData* myseqDatas = new SeqData[maxSizeSeq];
        for (int i = 0; i < maxSizeSeq; ++i)myseqDatas[i].pr = _pr;
        seqSet = myseqDatas;
        int posSeq = 0;
        for (int i = 1; i <= n; ++i) {
            nodes[i]->seq0_i = &myseqDatas[posSeq];
            nodes[i]->seqi_0 = &myseqDatas[posSeq + 1];
            nodes[i]->seqi_n = &myseqDatas[posSeq + 2];
            nodes[i]->seqn_i = &myseqDatas[posSeq + 3];
            posSeq += 4;
        }
        //only contain depot:
        seqDep = new SeqData(_pr);
        seqDep->init(0);
        for (int i = 1; i <= 7; ++i)valSeq[i] = new SeqData(_pr);
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
        valNode->seq0_i = seqDep;
        valNode->seqi_0 = seqDep;
        int u, v;
        //update seqdata (0->i, i->0):
        while (valNode != depot->pred)
        {    
            valNode = valNode->suc;
            v = valNode->idxClient;
            /*if (valNode->pred == depot)valNode->culCost = 0;
            else valNode->culCost = valNode->pred->culCost + pr->costs[valNode->pred->pred->idxClient][valNode->pred->idxClient][valNode->idxClient]; */
            valNode->seq0_i->concatOneAfter(valNode->pred->seq0_i, v);
            valNode->seqi_0->concatOneBefore(valNode->pred->seqi_0, v);
            valNode->posInSol = valNode->pred->posInSol + 1;
        }
        //update seqdata(i->n, n->i):
        valNode = depot;
        valNode->seqi_n = seqDep;
        valNode->seqn_i = seqDep;
        while (valNode != depot->suc)
        {
            valNode = valNode->pred;
            v = valNode->idxClient;
            valNode->seqn_i->concatOneAfter(valNode->suc->seqn_i, v);
            valNode->seqi_n->concatOneAfter(valNode->suc->seqi_n, v);
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
        reinitNegiborSet();
    }

    void reinitNegiborSet() {
        set<DI> sDis;
        int u, v, vSuc;
        for (int i = 1; i <= n; ++i) {
            sDis.clear();
            u = nodes[i]->idxClient;
            for (int j = 1; j <= n; j++)if (i != j) {
                v = nodes[j]->idxClient;
                vSuc = nodes[j]->suc->idxClient;
                if (vSuc == u)vSuc = nodes[j]->suc->suc->idxClient;
                sDis.insert(DI(pr->costs[u][v][vSuc], v));
            }
            nodes[i]->moves.clear();
            for (auto val : sDis) {
                nodes[i]->moves.push_back(val.sc);
                if (nodes[i]->moves.size() == pr->maxNeibor)break;                    
            }
        }        
    }

    // insert node u after node v
    void insertNode(Node* u, Node* v)
    {
        if (u->pred != v && u != v)
        {
            u->pred->suc = u->suc;
            u->suc->pred = u->pred;
            v->suc->pred = u;
            u->pred = v;
            u->suc = v->suc;
            v->suc = u;
            /*u->rou = v->rou;*/
        }
    }
    
    /*SeqData*/
    void constructSeqData(Node* startNode, Node* endNode, SeqData* resSeq) {
        if (startNode == endNode) {
            resSeq->init(startNode->idxClient);
            return;
        }
        resSeq->firstnode = startNode->idxClient;
        resSeq->lastnode = endNode->idxClient;
        resSeq->afterFiNode = startNode->suc->idxClient;
        resSeq->beforeLaNode = endNode->pred->idxClient;
        //cal cost based on seq0_i, seqi_n
        if (startNode->posInSol > endNode->posInSol) {
            resSeq->afterFiNode = startNode->pred->idxClient;
            resSeq->beforeLaNode = endNode->suc->idxClient;
            //symetric case:
            resSeq->cost = startNode->seq0_i - endNode->suc->seq0_i;
        }
        else {
            resSeq->cost = endNode->seq0_i - startNode->suc->seq0_i;
        }

    }
    //Or-opt
    bool move1() {
        //
        return true;
    }
    //String exchange
    //2-opt

    void setLocalValU() {        
        uSuc = nodeU->suc;
        uPred = nodeU->pred;
        nodeUIdx = nodeU->idxClient;
        predUIdx = uPred->idxClient;
        sucUIdx = uSuc->idxClient;
    }

    void setLocalValV() {
        vSuc = nodeV->suc;
        vPred = nodeV->pred;
        nodeVIdx = nodeV->idxClient;
        predVIdx = vPred->idxClient;
        sucVIdx = vSuc->idxClient;
    }

    void updateObj() {
        shuffle(ordNodeLs.begin(), ordNodeLs.end(), pr->Rng.generator);
        bool isFinished = false;
        while (!isFinished) {
            isFinished = true;
            for (int posU = 0; posU < ordNodeLs.size(); ++posU) {                
                nodeU = nodes[ordNodeLs[posU]];
                for (int posV = 0; posV < nodeU->moves.size(); ++posV) {
                    nodeV = nodes[nodeU->moves[posV]];
                    vSuc = nodeV->suc;
                    vPred = nodeV->pred;                    
                }
            }
        }
    }
    
};