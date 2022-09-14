#pragma once
#include "lib.h"
#include "Node.h"

class Solution {
public:
    vector<int> giantT;//giant tour
    vector<Node*> nodes;// for nodes
    Param* pr;
    double cost;// objective
    int n;// number of customers    
    bool isFinished;//used for LS
    string typeIns;
    vector<SeqData*> mySeq;
    SeqData* seqSet;
    SeqData* seqDep;
    SeqData* valSeq[12];
    //LNS-based pertubation part:
    vector<int> setIdRmv;
    set<DI> setRmv;
    //GA parts:
    vector<int> predecessors; // store prev node of each client  
    vector<int> successors; // store next node of each client  
    multiset<pair<double, Solution*>> indivsPerProximity;
    double biasedFitness; // Biased fitness of the solution
    vector<int> ordNodeLs; // using for LS
    Node* depot;
    Node* nodeU;
    Node* uPred;
    Node* uSuc;
    Node* nodeV;
    Node* vPred;
    Node* vSuc;
    Node* oriNodeU;
    Node* oriNodeV;
    Node* curNodeU;
    Node* curNodeV;
    int nodeUIdx, nodeVIdx, predUIdx, predVIdx, sucUIdx, sucVIdx;
    int nbMove1 = 0, nbMove2 = 0, nbMove3 = 0, nbMove4 = 0, nbMove5 = 0, nbMove6 = 0, nbMove7 = 0, nbMove8 = 0, nbMove9 = 0;
    int lastInsType;
    Solution(Param* _pr) {
        pr = _pr;
        n = pr->numLoc - 1;//not cointain depot
        for (int i = 0; i <= n + 1; ++i) {
            giantT.pb(0);
        }
        for (int i = 1; i <= n + 3; ++i)nodes.pb(new Node());
        setIdRmv.clear();
        for (int i = 1; i <= n; ++i) {
            nodes[i]->idxClient = i;
            ordNodeLs.push_back(i);
            setIdRmv.push_back(i);
        }
        ordNodeLs.push_back(n + 1);
        nodes[n + 1]->idxClient = 0;
        depot = nodes[n + 1];
        depot->isDepot = true;
        depot->posInSol = 0;
        nodes[0] = depot;
        int maxSizeSeq = 7 * (n + 1) + 1;
        SeqData* myseqDatas = new SeqData[maxSizeSeq];
        for (int i = 0; i < maxSizeSeq; ++i)myseqDatas[i].pr = _pr;
        seqSet = myseqDatas;
        int posSeq = 0;
        for (int i = 1; i <= n + 1; ++i) {
            nodes[i]->seq0_i = &myseqDatas[posSeq];
            nodes[i]->seqi_0 = &myseqDatas[posSeq + 1];
            nodes[i]->seqi_n = &myseqDatas[posSeq + 2];
            nodes[i]->seqn_i = &myseqDatas[posSeq + 3];
            nodes[i]->seqi_i = &myseqDatas[posSeq + 4];
            nodes[i]->seqi_j = &myseqDatas[posSeq + 5];
            nodes[i]->seqj_i = &myseqDatas[posSeq + 6];
            nodes[i]->seqi_i->init(i);
            posSeq += 7;
        }
        //only contain depot:
        seqDep = new SeqData(_pr);
        seqDep->init(0);
        depot->seqi_i = seqDep;
        for (int i = 1; i <= 7; ++i)valSeq[i] = new SeqData(_pr);
        predecessors.resize(n + 2);
        successors.resize(n + 2);
    }

    void genGiantT() {
        for (int i = 1; i <= n; ++i)giantT[i] = i;
        shuffle(giantT.begin() + 1, giantT.begin() + n + 1, pr->Rng.generator);
        setIdRmv.clear();
        //random_shuffle(giantT.begin() + 1, giantT.end());
    }

    void genGiantTOri() {
        for (int i = 1; i <= n; ++i)giantT[i] = i;
        setIdRmv.clear();
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
        /*for (int j = 1; j <= n+1; j++)
        {
            if (successors[j] != valSol->successors[j] && successors[j] != valSol->predecessors[j]) differences++;
            else if (successors[j] == valSol->successors[j] && predecessors[j] != valSol->predecessors[j]) differences++;
            else if (successors[j] == valSol->predecessors[j] && predecessors[j] != valSol->successors[j]) differences++;
        }
        return (double)differences / (double)(n + 1);*/
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
        //TODO: edit update after adding artificial depot.
        Node* valNode = depot;
        valNode->seq0_i = seqDep;
        valNode->seqi_0 = seqDep;
        valNode->seqi_j->concatOneAfter(valNode->seqi_i, valNode->suc->idxClient);
        valNode->seqj_i->concatOneAfter(valNode->suc->seqi_i, valNode->idxClient);
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
            valNode->seqi_j->concatOneAfter(valNode->seqi_i, valNode->suc->idxClient);
            valNode->seqj_i->concatOneAfter(valNode->suc->seqi_i, valNode->idxClient);
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
            valNode->seqi_n->concatOneBefore(valNode->suc->seqi_n, v);
        }
    }
    /**/
    // insert depot into its best position...
    //void calCost() {
    //    cost = pr->costs[giantT[n - 1]][giantT[n]][giantT[1]] + pr->costs[giantT[n]][giantT[1]][giantT[2]];
    //    for (int i = 1; i <= n - 2; ++i) {
    //        cost += pr->costs[giantT[i]][giantT[i + 1]][giantT[i + 2]];
    //    }
    //    int prvPos = n;
    //    int sucPos = 1;
    //    double minChangeCost = -(pr->costs[giantT[n - 1]][giantT[n]][giantT[1]] + pr->costs[giantT[n]][giantT[1]][giantT[2]])
    //        + (pr->costs[giantT[n]][giantT[0]][giantT[1]] + pr->costs[giantT[n - 1]][giantT[n]][giantT[0]] + pr->costs[giantT[0]][giantT[1]][giantT[2]]);        
    //    double valChangeCost;
    //    int prvI;
    //    int sucI1;
    //    for (int i = 1; i <= n - 1; ++i) {
    //        prvI = i - 1;
    //        if (prvI == 0)prvI = n;
    //        sucI1 = i + 2;
    //        if (sucI1 == n + 1)sucI1 = 1;
    //        valChangeCost = -(pr->costs[giantT[prvI]][giantT[i]][giantT[i + 1]] + pr->costs[giantT[i]][giantT[i + 1]][giantT[sucI1]])
    //            + (pr->costs[giantT[i]][giantT[0]][giantT[i + 1]] + pr->costs[giantT[prvI]][giantT[i]][giantT[0]] + pr->costs[giantT[0]][giantT[i + 1]][giantT[sucI1]]);
    //        if (valChangeCost < minChangeCost) {
    //            minChangeCost = valChangeCost;
    //            prvPos = i;
    //            sucPos = i + 1;
    //        }
    //    }
    //    cost += minChangeCost;
    //    //construct route
    //    nodes[giantT[1]]->pred = nodes[giantT[n]];
    //    nodes[giantT[n]]->suc = nodes[giantT[1]];
    //    for (int i = 2; i <= n; ++i) {
    //        int idV = giantT[i], idU = giantT[i - 1];
    //        nodes[idU]->suc = nodes[idV];
    //        nodes[idV]->pred = nodes[idU];
    //    }
    //    //insert depot into its best position:
    //    nodes[giantT[prvPos]]->suc = depot;
    //    depot->pred = nodes[giantT[prvPos]];
    //    nodes[giantT[sucPos]]->pred = depot;
    //    depot->suc = nodes[giantT[sucPos]];
    //    /*cost = pr->costs[giantT[n]][giantT[0]][giantT[1]] + pr->costs[giantT[n - 1]][giantT[n]][giantT[0]];
    //    for (int i = 0; i <= n - 2; ++i) {
    //        cost += pr->costs[giantT[i]][giantT[i + 1]][giantT[i + 2]];
    //    }
    //    depot->suc = nodes[giantT[1]];
    //    nodes[giantT[1]]->pred = depot;
    //    depot->pred = nodes[giantT[n]];
    //    nodes[giantT[n]]->suc = depot;
    //    for (int i = 2; i <= n; ++i) {
    //        int idV = giantT[i], idU = giantT[i - 1];
    //        nodes[idU]->suc = nodes[idV];
    //        nodes[idV]->pred = nodes[idU];
    //    }*/
    //    updateInfo();
    //    reinitNegiborSet();
    //    for (int i = 1; i <= n + 1; ++i) {
    //        predecessors[i] = nodes[i]->pred->idxClient;
    //        successors[i] = nodes[i]->suc->idxClient;
    //    }
    //}

    void calCost() {
        cost = pr->costs[giantT[n]][giantT[0]][giantT[1]] + pr->costs[giantT[n - 1]][giantT[n]][giantT[0]];
        //cost = pr->costs[giantT[n]][giantT[1]][giantT[2]] + pr->costs[giantT[n - 1]][giantT[n]][giantT[1]]; //artificial case
        for (int i = 0; i <= n - 2; ++i) {
            cost += pr->costs[giantT[i]][giantT[i + 1]][giantT[i + 2]];
        }
        depot->suc = nodes[giantT[1]];
        nodes[giantT[1]]->pred = depot;
        depot->pred = nodes[giantT[n]];
        nodes[giantT[n]]->suc = depot;
        for (int i = 2; i <= n; ++i) {
            int idV = giantT[i], idU = giantT[i - 1];
            nodes[idU]->suc = nodes[idV];
            nodes[idV]->pred = nodes[idU];
        }
        updateInfo();
        reinitNegiborSet();
        for (int i = 1; i <= n + 1; ++i) {
            predecessors[i] = nodes[i]->pred->idxClient;
            successors[i] = nodes[i]->suc->idxClient;
        }
    }

    bool checkSol() {
        int* dd = new int[n + 1];
        for (int i = 1; i <= n; ++i)dd[i] = 0;
        Node* val = depot;
        int pos = 0;
        do {            
            if (val->idxClient != 0) {
                dd[val->idxClient] = 1;
                pos++;
            }            
            val = val->suc;
        } while (val != depot);
        for (int i = 1; i <= n; ++i)if (dd[i] == 0) {
            return false;
        }
        if (pos != n) {
            return false;
        }
        delete[] dd;
        return true;
    }

    bool checkGiantT() {
        int* dd = new int[n + 1];
        for (int i = 1; i <= n; ++i)dd[i] = 0;
        for (int i = 1; i <= n; ++i)dd[giantT[i]] = 1;
        for (int i = 1; i <= n; ++i)if (dd[i] == 0) {
            delete[] dd;
            return false;
        }        
        delete[] dd;
        return true;        
    }

    void cvGiantT() {
        Node* val = depot;
        int pos = 0;
        do {            
            if (val->idxClient != 0)giantT[++pos] = val->idxClient;
            val = val->suc;
        } while (val != depot);        
        /*if (pr->isDebug) {
            int* dd = new int[n + 1];
            for (int i = 1; i <= n; ++i)dd[i] = 0;
            for (int i = 1; i <= n; ++i)dd[giantT[i]] = 1;
            for (int i = 1; i <= n; ++i)if (dd[i] == 0) {
                throw "Wrong here 1";
            }
            if (pos != n) {
                throw "Wrong here 2";
            }
            delete[] dd;
        }*/
    }

    double calCostWtUpdate() {
        cvGiantT();
        double res = pr->costs[giantT[n]][giantT[0]][giantT[1]] + pr->costs[giantT[n - 1]][giantT[n]][giantT[0]];
        for (int i = 0; i <= n - 2; ++i) {
            res += pr->costs[giantT[i]][giantT[i + 1]][giantT[i + 2]];
        }
        return res;
    }

    void reinitNegiborSet() {
        /*set<DI> sDis;
        int u, v, vSuc;
        for (int i = 1; i <= n; ++i) {
            sDis.clear();
            u = nodes[i]->idxClient;
            for (int j = 1; j <= n + 1; j++)if (i != j) {
                v = nodes[j]->idxClient;
                vSuc = nodes[j]->suc->idxClient;
                if (vSuc == u)vSuc = nodes[j]->suc->suc->idxClient;
                sDis.insert(DI(pr->costs[v][u][vSuc], j));
            }
            if(nodes[i]->moves.size() != 0)nodes[i]->moves.clear();
            for (auto val : sDis) {
                nodes[i]->moves.push_back(val.sc);
                if (nodes[i]->moves.size() == pr->maxNeibor)break;
            }
        }        */
    }

    //swap two nodes u and v
    void swapNode(Node* u, Node* v) {
        Node* oldVPred = v->pred;
        Node* oldVSuc = v->suc;
        Node* oldUPred = u->pred;
        Node* oldUSuc = u->suc;

        oldUPred->suc = v;
        oldUSuc->pred = v;
        oldVPred->suc = u;
        oldVSuc->pred = u;

        u->pred = oldVPred;
        u->suc = oldVSuc;
        v->pred = oldUPred;
        v->suc = oldUSuc;
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

    // insert node u (not belong to tour) after node v
    void insertNodeNotInTour(Node* u, Node* v)
    {
        v->suc->pred = u;
        u->pred = v;
        u->suc = v->suc;
        v->suc = u;
    }

    void showR() {
        Node* val = this->depot;
        do {
            cout << val->idxClient << " ";
            val = val->suc;
        } while (val != this->depot);
        cout << "\n";
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
            resSeq->cost = startNode->seq0_i->cost - endNode->suc->seq0_i->cost;
        }
        else {
            resSeq->cost = endNode->seq0_i->cost - startNode->suc->seq0_i->cost;
        }
    }
    //relocate v after u
    bool move1() {
        if (pr->isDebug)cout << "move1\n";
        //nodeU in this case can be depot (arrival depot)
        if (nodeV == uSuc)return false;
        if (nodeV->posInSol < nodeU->posInSol) {
            valSeq[1]->copy(vPred->seq0_i);
            constructSeqData(vSuc, nodeU, valSeq[2]);
            valSeq[3]->copy(nodeV->seqi_i);
            valSeq[4]->copy(uSuc->seqi_n);
        }
        else {
            valSeq[1]->copy(nodeU->seq0_i);
            valSeq[2]->copy(nodeV->seqi_i);
            constructSeqData(uSuc, vPred, valSeq[3]);
            valSeq[4]->copy(vSuc->seqi_n);
        }
        mySeq.clear();
        for (int i = 1; i <= 4; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        insertNode(nodeV, nodeU);
        updateInfo();
        nbMove1++;
        isFinished = false;
        cost = newCost;
        if (pr->isDebug && abs(cost - calCostWtUpdate()) > MY_EPSILON) {
            throw "bug in move1";
        }
        return true;
    }

    //relocate v vSuc after u
    bool move2() {
        if (pr->isDebug)cout << "move2\n";
        //nodeU in this case can be depot (arrival depot)
        if (nodeV == uSuc || nodeU == vSuc || vSuc->isDepot)return false;
        if (nodeV->posInSol < nodeU->posInSol) {
            valSeq[1]->copy(vPred->seq0_i);
            constructSeqData(vSuc->suc, nodeU, valSeq[2]);
            valSeq[3]->copy(nodeV->seqi_j);
            valSeq[4]->copy(uSuc->seqi_n);
        }
        else {
            valSeq[1]->copy(nodeU->seq0_i);
            valSeq[2]->copy(nodeV->seqi_j);
            constructSeqData(uSuc, vPred, valSeq[3]);
            valSeq[4]->copy(vSuc->suc->seqi_n);
        }
        mySeq.clear();
        for (int i = 1; i <= 4; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        insertNode(nodeV, nodeU);
        insertNode(vSuc, nodeV);
        updateInfo();
        nbMove2++;
        isFinished = false;
        cost = newCost;
        if (pr->isDebug && abs(cost - calCostWtUpdate()) > MY_EPSILON) {
            throw "bug in move2";
        }
        return true;
    }
    //relocate v vSuc vSucSuc after u
    bool move9() {
        if (pr->isDebug)cout << "move9\n";
        //nodeU in this case can be depot (arrival depot)
        if (nodeV == uSuc || nodeU == vSuc ||  vSuc->isDepot
            || nodeU == vSuc->suc || vSuc->suc->isDepot)return false;
        if (nodeV->posInSol < nodeU->posInSol) {
            valSeq[1]->copy(vPred->seq0_i);
            constructSeqData(vSuc->suc->suc, nodeU, valSeq[2]);
            constructSeqData(nodeV, vSuc->suc, valSeq[3]);            
            valSeq[4]->copy(uSuc->seqi_n);
        }
        else {
            valSeq[1]->copy(nodeU->seq0_i);            
            constructSeqData(nodeV, vSuc->suc, valSeq[2]);
            constructSeqData(uSuc, vPred, valSeq[3]);
            valSeq[4]->copy(vSuc->suc->suc->seqi_n);
        }
        mySeq.clear();
        for (int i = 1; i <= 4; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        Node* vSucSuc = vSuc->suc;
        insertNode(nodeV, nodeU);
        insertNode(vSuc, nodeV);        
        insertNode(vSucSuc, vSuc);
        updateInfo();
        nbMove9++;
        isFinished = false;
        cost = newCost;
        if (pr->isDebug && abs(cost - calCostWtUpdate()) > MY_EPSILON) {            
            throw "bug in move9";
        }
        return true;
    }
    //relocate vSuc v after U
    bool move3() {
        if (pr->isDebug)cout << "move3\n";
        //nodeU in this case can be depot (arrival depot)
        if (nodeV == uSuc || nodeU == vSuc || vSuc->isDepot)return false;
        if (nodeV->posInSol < nodeU->posInSol) {
            valSeq[1]->copy(vPred->seq0_i);
            constructSeqData(vSuc->suc, nodeU, valSeq[2]);
            valSeq[3]->copy(nodeV->seqj_i);
            valSeq[4]->copy(uSuc->seqi_n);
        }
        else {
            valSeq[1]->copy(nodeU->seq0_i);
            valSeq[2]->copy(nodeV->seqj_i);
            constructSeqData(uSuc, vPred, valSeq[3]);
            valSeq[4]->copy(vSuc->suc->seqi_n);
        }
        mySeq.clear();
        for (int i = 1; i <= 4; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        insertNode(vSuc, nodeU);
        insertNode(nodeV, vSuc);
        updateInfo();
        nbMove3++;
        isFinished = false;
        cost = newCost;
        if (pr->isDebug && abs(cost - calCostWtUpdate()) > MY_EPSILON) {
            throw "bug in move3";
        }
        return true;
    }
    //Swaps moves:
    //swap u and v (u->pos < v->pos)    
    bool move4() {
        if (pr->isDebug)cout << "move4\n";
        if (nodeU == vPred || nodeU == vSuc)return false;
        valSeq[1]->copy(uPred->seq0_i);
        valSeq[2]->copy(nodeV->seqi_i);
        constructSeqData(uSuc, vPred, valSeq[3]);
        valSeq[4]->copy(nodeU->seqi_i);
        valSeq[5]->copy(vSuc->seqi_n);
        mySeq.clear();
        for (int i = 1; i <= 5; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        swapNode(nodeU, nodeV);
        updateInfo();
        nbMove4++;
        isFinished = false;
        cost = newCost;
        if (pr->isDebug && abs(cost - calCostWtUpdate()) > MY_EPSILON) {
            throw "bug in move4";
        }
        return true;
    }

    //swap (u, uSuc) and v    
    bool move5() {
        if (pr->isDebug)cout << "move5\n";
        if (nodeU == vPred || uSuc == vPred || nodeU == vSuc || uSuc->isDepot)return false;
        if (nodeU->posInSol < nodeV->posInSol) {
            valSeq[1]->copy(uPred->seq0_i);
            valSeq[2]->copy(nodeV->seqi_i);
            constructSeqData(uSuc->suc, vPred, valSeq[3]);
            valSeq[4]->copy(nodeU->seqi_j);
            valSeq[5]->copy(vSuc->seqi_n);
        }
        else {
            valSeq[1]->copy(vPred->seq0_i);
            valSeq[2]->copy(nodeU->seqi_j);
            constructSeqData(vSuc, uPred, valSeq[3]);
            valSeq[4]->copy(nodeV->seqi_i);
            valSeq[5]->copy(uSuc->suc->seqi_n);
        }
        mySeq.clear();
        for (int i = 1; i <= 5; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        swapNode(nodeU, nodeV);
        insertNode(uSuc, nodeU);
        updateInfo();
        nbMove5++;
        isFinished = false;
        cost = newCost;
        if (pr->isDebug && abs(cost - calCostWtUpdate()) > MY_EPSILON) {
            throw "bug in move5";
        }
        return true;
    }

    //swap u and (v, vSuc)
    bool move8() {
        if (pr->isDebug)cout << "move8\n";
        if (nodeV == uPred || vSuc == uPred || nodeV == uSuc || vSuc->isDepot)return false;
        if (nodeU->posInSol < nodeV->posInSol) {
            valSeq[1]->copy(uPred->seq0_i);
            valSeq[2]->copy(nodeV->seqi_j);
            constructSeqData(uSuc, vPred, valSeq[3]);
            valSeq[4]->copy(nodeU->seqi_i);
            valSeq[5]->copy(vSuc->suc->seqi_n);
        }
        else {
            valSeq[1]->copy(vPred->seq0_i);
            valSeq[2]->copy(nodeU->seqi_i);
            constructSeqData(vSuc->suc, uPred, valSeq[3]);
            valSeq[4]->copy(nodeV->seqi_j);
            valSeq[5]->copy(uSuc->seqi_n);
        }
        mySeq.clear();
        for (int i = 1; i <= 5; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        swapNode(nodeU, nodeV);
        insertNode(vSuc, nodeV);
        updateInfo();
        nbMove8++;
        isFinished = false;
        cost = newCost;
        if (pr->isDebug && abs(cost - calCostWtUpdate()) > MY_EPSILON) {
            throw "bug in move8";
        }
        return true;
    }

    //swap (u, Suc) and (v, vSuc) (u->pos < v->pos)
    bool move6() {
        if (pr->isDebug)cout << "move6\n";
        if (uSuc->isDepot || vSuc->isDepot || vSuc == uPred || nodeU == vSuc || uSuc == nodeV || nodeV == uSuc->suc)return false;
        valSeq[1]->copy(uPred->seq0_i);
        valSeq[2]->copy(nodeV->seqi_j);
        constructSeqData(uSuc->suc, vPred, valSeq[3]);
        valSeq[4]->copy(nodeU->seqi_j);
        valSeq[5]->copy(vSuc->suc->seqi_n);
        mySeq.clear();
        for (int i = 1; i <= 5; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        swapNode(nodeU, nodeV);
        swapNode(uSuc, vSuc);
        updateInfo();
        nbMove6++;
        isFinished = false;
        cost = newCost;
        if (pr->isDebug && abs(cost - calCostWtUpdate()) > MY_EPSILON) {
            throw "bug in move6";
        }
        return true;
    }

    //2-opt
    //reverse u...v to v...u
    //...oriUuv... and ...oriUuxv... are trivial cases
    bool move7(bool isTurn) {
        if (pr->isDebug)cout << "move7\n";
        if (uSuc == nodeV || uSuc == vPred)return false;
        double newCost;
        if (!isTurn) {//oriUu...v....
            /*if (nodeU->posInSol > nodeV->posInSol)return false;            */
            valSeq[1]->copy(uPred->seq0_i);
            constructSeqData(nodeV, nodeU, valSeq[2]);
            valSeq[3]->copy(vSuc->seqi_n);
            mySeq.clear();
            for (int i = 1; i <= 3; ++i)mySeq.push_back(valSeq[i]);
            newCost = seqDep->evaluation(mySeq);
        }
        else {//v...oriUu
            valSeq[1]->copy(nodeU->seqn_i);
            constructSeqData(vSuc, uPred, valSeq[2]);
            valSeq[3]->copy(nodeV->seqi_0);
            mySeq.clear();
            for (int i = 1; i <= 3; ++i)mySeq.push_back(valSeq[i]);
            newCost = seqDep->evaluation(mySeq);
        }
        if (newCost - cost > -MY_EPSILON)return false;
        //update route structure:
        Node* startNode = uPred;
        Node* endNode = nodeV;
        while (startNode != nodeU)
        {
            Node* tempNode = endNode->pred;
            startNode->suc = endNode;
            endNode->pred = startNode;
            startNode = endNode;
            endNode = tempNode;
        }
        nodeU->suc = vSuc;
        vSuc->pred = nodeU;
        updateInfo();
        nbMove7++;
        isFinished = false;
        cost = newCost;
        if (pr->isDebug && abs(cost - calCostWtUpdate()) > MY_EPSILON) {
            throw "bug in move7";
        }
        return true;
    }

    ////2-opt
    //bool move7() {
    //    if (pr->isDebug)cout << "move7\n";
    //    if (nodeU->posInSol > nodeV->posInSol)return false;
    //    if (uSuc == nodeV)return false;
    //    valSeq[1]->copy(nodeU->seq0_i);
    //    constructSeqData(nodeV, uSuc, valSeq[2]);
    //    valSeq[3]->copy(vSuc->seqi_n);
    //    mySeq.clear();
    //    for (int i = 1; i <= 3; ++i)mySeq.push_back(valSeq[i]);
    //    double newCost = seqDep->evaluation(mySeq);
    //    if (newCost - cost > -MY_EPSILON)return false;
    //    //update route structure:
    //    Node* startNode = nodeU;
    //    Node* endNode = nodeV;
    //    while (startNode != uSuc)
    //    {
    //        Node* tempNode = endNode->pred;
    //        startNode->suc = endNode;
    //        endNode->pred = startNode;
    //        startNode = endNode;
    //        endNode = tempNode;
    //    }
    //    uSuc->suc = vSuc;
    //    vSuc->pred = uSuc;
    //    updateInfo();
    //    nbMove7++;
    //    isFinished = false;
    //    cost = newCost;
    //    if (pr->isDebug && abs(cost - calCostWtUpdate()) > MY_EPSILON) {
    //        throw "bug in move7";
    //    }
    //    return true;
    //}

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
        int oriCost = cost;
        shuffle(ordNodeLs.begin(), ordNodeLs.end(), pr->Rng.generator);
        isFinished = false;
        bool isTurn = false;
        int bestI1, bestI2, bestJ1, bestJ2;                
        while (!isFinished) {
            isFinished = true;
            //cout << "updated cost: " << cost << "\n";
            for (int posU = 0; posU < ordNodeLs.size(); ++posU) {
                nodeU = nodes[ordNodeLs[posU]];
                oriNodeU = nodeU;
                int idNodeU = ordNodeLs[posU];
                if (idNodeU == n + 1)idNodeU = 0;
                int idPrvU = nodeU->pred->idxClient;
                if (pr->Rng.getIntRand() % pr->maxNeibor == 0) {
                    shuffle(pr->correlatedNodes[idPrvU][idNodeU].begin(), pr->correlatedNodes[idPrvU][idNodeU].end(), pr->Rng.generator);
                }
                for (int posV = 0; posV < pr->correlatedNodes[idPrvU][idNodeU].size(); ++posV) {
                    isTurn = false;
                    nodeV = nodes[pr->correlatedNodes[idPrvU][idNodeU][posV]];
                    oriNodeV = nodeV;
                    nodeU = oriNodeU;
                    setLocalValU();
                    setLocalValV();
                    if (move1())continue;//relocate 1
                    if (move2())continue;//relocate 2
                    //if (move9())continue;//relocate 3
                    //if (move3())continue;//relocate 2 (reversed)                    
                    if (nodeU->suc == nodeV)continue;
                    nodeU = nodeU->suc;
                    setLocalValU();
                    if (!nodeU->isDepot) {
                        if (move5())continue;// swap 2, 1                    
                        if (move8())continue;// swap 1, 2                    
                    }
                    if (nodeU->posInSol > nodeV->posInSol) {
                        nodeV = nodeU;
                        nodeU = oriNodeV;
                        setLocalValU();
                        setLocalValV();
                        isTurn = true;
                    }
                    if (!nodeU->isDepot) {
                        if (move4())continue;//swap 1, 1                    
                        if (move6())continue;//swap 2, 2          
                    }
                    if (isTurn) {
                        swap(nodeU, nodeV);
                        setLocalValU();
                        setLocalValV();
                    }
                    if (move7(isTurn)) continue;
                }
            }                     
        }
        double oldCost = cost;        
        cvGiantT();
        calCost();       
        if (pr->use4Opt) {
            //double test4opt = doubleBridge(false, bestI1, bestI2, bestJ1, bestJ2);
            if (n <= 5)return;// not use 4-opt for 5 customers
            double test4opt = fastDoubleBridge(bestI1, bestI2, bestJ1, bestJ2);
            //cout << oldCost << " " << cost << " " << test4opt << "\n";
            //if (oldCost - cost > MY_EPSILON)updateObj();
            if (test4opt > -MY_EPSILON) return;
            //doubleBridge(true, bestI1, bestI2, bestJ1, bestJ2);
            apply4Opt(bestI1, bestI2, bestJ1, bestJ2);
            updateObj();
        }
    }

    ///ELSALGO:
    void exchange() {
        int posU = pr->Rng.getNumInRan(1, n);
        int posV = pr->Rng.getNumInRan(1, n);
        while (posU == posV)
        {
            posV = pr->Rng.getNumInRan(1, n);
        }
        swap(giantT[posU], giantT[posV]);
    }

    void interchange() {
        int posU = pr->Rng.getNumInRan(1, n);
        int posV = pr->Rng.getNumInRan(1, n);
        while (posU == posV)
        {
            posV = pr->Rng.getNumInRan(1, n);
        }
        if (posU < posV) {
            for (int i = posU + 1; i <= posV; ++i) {
                swap(giantT[i], giantT[i - 1]);
            }
        }
        else {
            for (int i = posU - 1; i >= posV; --i) {
                swap(giantT[i], giantT[i + 1]);
            }
        }
    }

    void mutate(int numP) {
        for (int i = 1; i <= numP; ++i) {
            exchange();
        }
    }

    void exportGiantT() {
        pr->fileOut << "giant tour:\n";
        for (int i = 0; i <= n; ++i)pr->fileOut << giantT[i] << ", ";
    }

    DI findBestPos(Node* insNod) {
        //check the difference of cost between before and after
        Node* valNode = depot->pred;
        double newMinCost = oo;
        int bestInsPos = -1;
        double valNewCost;
        while (true)
        {
            valNode = valNode->suc;
            valSeq[1]->copy(valNode->seq0_i);
            valSeq[2]->copy(insNod->seqi_i);
            valSeq[3]->copy(valNode->suc->seqi_n);
            mySeq.clear();
            for (int j = 1; j <= 3; ++j)mySeq.push_back(valSeq[j]);
            valNewCost = seqDep->evaluation(mySeq);
            if (newMinCost > valNewCost) {
                newMinCost = valNewCost;
                bestInsPos = valNode->idxClient;
                if (valNode->isDepot)bestInsPos = n + 1;
            }
            if (valNode == depot->pred)break;
        }
        return DI(newMinCost, bestInsPos);
    }   

    //construction heuristics:
    void cheapestIns() {        
        int* check = new int[n + 1];
        for (int i = 1; i <= n; ++i)check[i] = 0;
        int nbIns = setIdRmv.size();        
        if (nbIns == n) { //first construction case:
            depot->suc = depot;
            depot->pred = depot;
            int u = pr->Rng.getNumInRan(1, n);
            /*int v = pr->Rng.getNumInRan(1, n);
            while (u == v)
            {
                v = pr->Rng.getNumInRan(1, n);
            }*/
            check[u] = 1;
            //check[v] = 1;
            depot->suc = nodes[u];
            depot->pred = nodes[u];
            nodes[u]->suc = depot;
            nodes[u]->pred = depot;
            //insertNodeNotInTour(nodes[v], nodes[u]);
            updateInfo();
            /*int idDepot = depot->idxClient, idNodeU = nodes[u]->idxClient, idNodeV = nodes[v]->idxClient;
            cost = pr->costs[idDepot][idNodeU][idNodeV] + pr->costs[idNodeU][idNodeV][idDepot] + pr->costs[idNodeV][idDepot][idNodeU];
            nbIns -= 2;*/
            nbIns--;
            cost = 0;
        }
        double minCost = oo;
        int bestId = -1;
        int bestInsPos = -1;        
        for (int i = 1; i <= nbIns; ++i) {             
            minCost = oo;
            for (auto idIns: setIdRmv)if (check[idIns] == 0) {
                DI valInfo = findBestPos(nodes[idIns]);
                if (minCost > valInfo.first) {
                    minCost = valInfo.first;
                    bestId = idIns;
                    bestInsPos = valInfo.second;
                }
            }            
            check[bestId] = 1;
            insertNodeNotInTour(nodes[bestId], nodes[bestInsPos]);
            updateInfo();            
            cost = minCost;                       
        }        
        cvGiantT();
        calCost();
        setIdRmv.clear();
        delete[] check;
    }
    
    void recoverSetRmv() {
        setIdRmv.clear();
        for (int i = 1; i <= n; ++i) {            
            setIdRmv.push_back(i);
        }
    }
    void randomIns() {
        int nbIns = setIdRmv.size();
        int idFirst = -1;
        if (nbIns == n) { //first construction case:
            depot->suc = depot;
            depot->pred = depot;
            int u = pr->Rng.getNumInRan(1, n);
            idFirst = u;
            /*int v = pr->Rng.getNumInRan(1, n);
            while (u == v)
            {
                v = pr->Rng.getNumInRan(1, n);
            }*/            
            //check[v] = 1;
            depot->suc = nodes[u];
            depot->pred = nodes[u];
            nodes[u]->suc = depot;
            nodes[u]->pred = depot;
            //insertNodeNotInTour(nodes[v], nodes[u]);
            updateInfo();
            /*int idDepot = depot->idxClient, idNodeU = nodes[u]->idxClient, idNodeV = nodes[v]->idxClient;
            cost = pr->costs[idDepot][idNodeU][idNodeV] + pr->costs[idNodeU][idNodeV][idDepot] + pr->costs[idNodeV][idDepot][idNodeU];
            nbIns -= 2;*/
            nbIns--;
            cost = 0;
        }
        shuffle(setIdRmv.begin(), setIdRmv.end(), pr->Rng.generator);
        for (auto idIns : setIdRmv) {            
            if (idIns == idFirst)continue;
            DI valInfo = findBestPos(nodes[idIns]);
            insertNodeNotInTour(nodes[idIns], nodes[valInfo.sc]);
            updateInfo();
            cost = valInfo.ft;
        }        
        cvGiantT();
        calCost();
        setIdRmv.clear();
    }

    double calEdgeCost(Node* u, Node* v, Node* k) {
        return pr->costs[u->idxClient][v->idxClient][k->idxClient];
    }

    void nearestIns() {        
        int* check = new int[n + 1];
        for (int i = 1; i <= n; ++i)check[i] = 0;
        int nbIns = setIdRmv.size();
        int idFirst = -1;
        if (nbIns == n) { //first construction case:
            depot->suc = depot;
            depot->pred = depot;
            int u = pr->Rng.getNumInRan(1, n);
            idFirst = u;           
            depot->suc = nodes[u];
            depot->pred = nodes[u];
            nodes[u]->suc = depot;
            nodes[u]->pred = depot;            
            check[u] = 1;            
            updateInfo();            
            nbIns--;
            cost = 0;
        }
        Node* firstT = depot;
        Node* lastT = depot->pred;
        double minCost = oo;
        int bestId = -1;
        int bestInsPos = -1;
        for (int i = 1; i <= nbIns; ++i) {
            minCost = oo;
            for (auto idIns : setIdRmv)if (check[idIns] == 0) {                
                if (calEdgeCost(nodes[idIns], firstT, firstT->suc) < minCost) {
                    minCost = calEdgeCost(nodes[idIns], firstT, firstT->suc);
                    bestId = idIns;
                    bestInsPos = firstT->idxClient;
                }
                if (calEdgeCost(lastT->pred, lastT, nodes[idIns]) < minCost) {
                    minCost = calEdgeCost(lastT->pred, lastT, nodes[idIns]);
                    bestId = idIns;
                    bestInsPos = lastT->idxClient;
                }                
            }
            if (bestInsPos == firstT->idxClient){                        
                addEdge(nodes[bestId], firstT);
                firstT = nodes[bestId];
            }
            else {                
                addEdge(lastT, nodes[bestId]);
                lastT = nodes[bestId];
            }
            check[bestId] = 1;
            /*insertNodeNotInTour(nodes[bestId], nodes[bestInsPos]);
            updateInfo();*/
            /*cost = minCost;*/
        }
        addEdge(lastT, firstT);
        cvGiantT();
        calCost();
        setIdRmv.clear();
        delete[] check;
    }


    void rmvNode(Node* u) {
        u->pred->suc = u->suc;
        u->suc->pred = u->pred;
    }

    void calNewCost() {
        if (depot->pred->posInSol < 2) {
            cost = 0;
            return;
        }
        cost = depot->pred->seq0_i->cost + pr->costs[depot->pred->pred->idxClient][depot->pred->idxClient][depot->idxClient]
            + pr->costs[depot->pred->idxClient][depot->idxClient][depot->suc->idxClient];
    }

    void randomRmv(int nbRmv) {
        if (setIdRmv.size() == n)setIdRmv.clear();
        int* rmvIdx = new int[n + 2];
        for (int i = 1; i <= n; ++i)rmvIdx[i] = i;
        shuffle(rmvIdx + 1, rmvIdx + n + 1, pr->Rng.generator);
        for (int i = 1; i <= nbRmv; ++i) {
            rmvNode(nodes[rmvIdx[i]]);
            setIdRmv.push_back(rmvIdx[i]);
        }
        updateInfo();
        calNewCost();
        delete[] rmvIdx;
    }

    void blockRmv(int nbRmv) {
        if (setIdRmv.size() == n)setIdRmv.clear();
        int idRmv = pr->Rng.getNumInRan(1, n);        
        for (int i = 1; i <= nbRmv; ++i) {
            rmvNode(nodes[giantT[idRmv]]);
            setIdRmv.push_back(giantT[idRmv]);
            idRmv++;
            if (idRmv == n + 1)idRmv = 1;
        }
        updateInfo();
        calNewCost();        
    }

    void changeRmvOnSet(Node* changedNode) {
        if (changedNode->isDepot)return;
        setRmv.erase(setRmv.find(DI(changedNode->rmvCost, changedNode->idxClient)));
        changedNode->calCurRmvCost(pr);
        setRmv.insert(DI(changedNode->rmvCost, changedNode->idxClient));
    }
    void worstRmv(int nbRmv, int worstDeg) {                        
        if (setIdRmv.size() == n)setIdRmv.clear();
        setRmv.clear();
        for (int i = 1; i <= n; ++i) {                        
            nodes[i]->calCurRmvCost(pr);            
            setRmv.insert(DI(nodes[i]->rmvCost, i));            
        }
        double randRate;
        int idInSet, idRmv, countId;      
        DI valRmv;
        for (int i = 1; i <= nbRmv; ++i) {            
            randRate = pow(pr->Rng.genRealInRang01(), worstDeg) * setRmv.size();
            idInSet = min((int)floor(randRate), (int)setRmv.size() - 1);
            countId = 0;
            //DI valRmv = *setRmv.begin();            
            for (auto val : setRmv) {
                if (countId == idInSet) {
                    valRmv = val;
                    break;
                }
                countId++;
            }
            idRmv = valRmv.second;
            rmvNode(nodes[idRmv]);            
            setIdRmv.push_back(idRmv);            
            changeRmvOnSet(nodes[idRmv]->pred);
            changeRmvOnSet(nodes[idRmv]->pred->pred);            
            changeRmvOnSet(nodes[idRmv]->suc);            
            changeRmvOnSet(nodes[idRmv]->suc->suc);
            setRmv.erase(valRmv);
        }
        updateInfo();
        calNewCost();
        setRmv.clear();
    }
    
    void addEdge(Node* U, Node* V) {        
        U->suc = V;
        V->pred = U;
    }

    double cal4optCost(int i1, int i2, int j1, int j2) {
        mySeq.clear();
        valSeq[1]->copy(nodes[giantT[i1]]->seq0_i);
        constructSeqData(nodes[giantT[j1 + 1]], nodes[giantT[j2]], valSeq[2]);
        constructSeqData(nodes[giantT[i2 + 1]], nodes[giantT[j1]], valSeq[3]);
        constructSeqData(nodes[giantT[i1 + 1]], nodes[giantT[i2]], valSeq[4]);
        valSeq[5]->copy(nodes[giantT[j2 + 1]]->seqi_n);
        for (int i = 1; i <= 5; ++i)mySeq.push_back(valSeq[i]);
        return seqDep->evaluation(mySeq) - cost;
    }

    double doubleBridge(bool isChanged,int &valI1, int& valI2, int& valJ1, int& valJ2) {        
        int bestI1 = valI1, bestI2 = valI2, bestJ1 = valJ1, bestJ2 = valJ2;
        double changedCost = oo;
        if (!isChanged) {            
            for (int i1 = 0; i1 <= n; ++i1)
                for (int i2 = i1 + 1; i2 <= n; ++i2)
                    for (int j1 = i2 + 1; j1 <= n; ++j1)
                        for (int j2 = j1 + 1; j2 <= n - 1; ++j2) {                                                        
                            double valChangedCost = cal4optCost(i1, i2, j1, j2);
                            if (changedCost - valChangedCost > -MY_EPSILON) {
                                changedCost = valChangedCost;
                                bestI1 = i1;
                                bestI2 = i2;
                                bestJ1 = j1;
                                bestJ2 = j2;
                            }
                        }
            valI1 = bestI1; valI2 = bestI2; valJ1 = bestJ1; valJ2 = bestJ2;
        }
        if (isChanged) {            
            apply4Opt(bestI1, bestI2, bestJ1, bestJ2);
        }
        return changedCost;
        /*int step = (n + 1) / 4 - 1;
        int pos0 = 0;
        int pos0_p = n;

        int pos1 = 1 + pr->Rng.getNumInRan(0, step);
        int pos1_p = pos1 - 1;
        
        int pos2 = 1 + pos1 + pr->Rng.getNumInRan(0, step);
        int pos2_p = pos2 - 1;

        int pos3 = 1 + pos2 + pr->Rng.getNumInRan(0, step);
        int pos3_p = pos3 - 1;

        addEdge(nodes[giantT[pos2_p]], nodes[giantT[pos0]]);
        addEdge(nodes[giantT[pos0_p]], nodes[giantT[pos2]]);
        addEdge(nodes[giantT[pos1_p]], nodes[giantT[pos3]]);
        addEdge(nodes[giantT[pos3_p]], nodes[giantT[pos1]]);
        
        cvGiantT();
        calCost();        */
    }
    void apply4Opt(const int& bestI1, const int& bestI2, const int& bestJ1, const int& bestJ2) {
        addEdge(nodes[giantT[bestI1]], nodes[giantT[bestJ1 + 1]]);
        addEdge(nodes[giantT[bestJ2]], nodes[giantT[bestI2 + 1]]);
        addEdge(nodes[giantT[bestJ1]], nodes[giantT[bestI1 + 1]]);
        addEdge(nodes[giantT[bestI2]], nodes[giantT[bestJ2 + 1]]);
        cvGiantT();
        /*if (!checkGiantT()) {
            cout << "bug here\n";
        }*/
        calCost();
    }
    
    double costEGiantT(const int &i, const int &j, const int &k) {       
        return pr->costs[giantT[i]][giantT[j]][giantT[k]];
    }    
    double cal2DOpt(const int &valI, const int &valJ) {        
        //valI-1 = -1 -> valI = n
        //valJ+2 = n+1 -> valJ = 0
        int valIsub1 = valI - 1;
        if (valIsub1 == -1) valIsub1 = n;
        return costEGiantT(valIsub1, valI, valJ + 1) + costEGiantT(valI, valJ + 1, valJ + 2)
            + costEGiantT(valJ, valI + 1, valI + 2) + costEGiantT(valJ - 1, valJ, valI + 1)
            - costEGiantT(valIsub1, valI, valI + 1) - costEGiantT(valI, valI + 1, valI + 2)
            - costEGiantT(valJ - 1, valJ, valJ + 1) - costEGiantT(valJ, valJ + 1, valJ + 2);
    }    

    double fastDoubleBridge(int& valI1, int& valI2, int& valJ1, int& valJ2) {
    //double fastDoubleBridge() {
        int bestI1 = 0, bestI2 = 0, bestJ1 = 0, bestJ2 = 0;
        double changedCost = oo;
        for (int j = 4; j <= n - 3; ++j) {
            pr->predSup[2][j] = 0;//save i1 for best i2, j1
            pr->fSup[2][j] = cal2DOpt(0, j);
            for (int i = 3; i <= j - 2; ++i) {                
                pr->predSup[i][j] = i - 2;
                pr->fSup[i][j] = cal2DOpt(i - 2, j);
                if (pr->fSup[i][j] - pr->fSup[i - 1][j] > -MY_EPSILON) {
                    pr->fSup[i][j] = pr->fSup[i - 1][j];
                    pr->predSup[i][j] = pr->predSup[i - 1][j];
                }
            }
        }
        double valCost;
        for (int i = 2; i <= n - 5; ++i) {
            pr->predMain[i][i + 4] = i + 2;//save j1 for best i2, j2
            pr->fMain[i][i + 4] = pr->fSup[i][i + 2];
            for (int j = i + 5; j <= n - 1; ++j) {
                pr->predMain[i][j] = j - 2;
                pr->fMain[i][j] = pr->fSup[i][j - 2];
                if (pr->fMain[i][j] - pr->fMain[i][j - 1] > -MY_EPSILON) {
                    pr->fMain[i][j] = pr->fMain[i][j - 1];
                    pr->predMain[i][j] = pr->predMain[i][j - 1];
                }
                valCost = pr->fMain[i][j] + cal2DOpt(i, j);
                if (changedCost - valCost > -MY_EPSILON) {
                    changedCost = valCost;
                    bestI2 = i;
                    bestJ2 = j;
                    bestJ1 = pr->predMain[i][j];
                    bestI1 = pr->predSup[i][bestJ1];
                }
            }
        }
        /*for (int i1 = 0; i1 <= n; ++i1)
            for (int i2 = i1 + 1; i2 <= n; ++i2)
                for (int j1 = i2 + 1; j1 <= n; ++j1)
                    for (int j2 = j1 + 1; j2 <= n - 1; ++j2) */
        int i1, i2, j1, j2;
        
        //O(n^3)
        //for (i1 = 0; i1 < n; ++i1) {
        //    //i1+1 = i2
        //    i2 = i1 + 1;
        //    for (j1 = i2 + 1; j1 < n; ++j1)
        //        for (int j2 = j1 + 1; j2 < n; ++j2) {
        //            valCost = cal4optCost(i1, i2, j1, j2);
        //            if (changedCost - valCost > -MY_EPSILON) {
        //                changedCost = valCost;
        //                bestI1 = i1;
        //                bestI2 = i2;
        //                bestJ1 = j1;
        //                bestJ2 = j2;
        //            }
        //        }
        //    
        //    //i2+1 = j1
        //    for (i2 = i1 + 1; i2 < n; ++i2) {
        //        j1 = i2 + 1;
        //        for (j2 = j1 + 1; j2 < n; ++j2) {
        //            valCost = cal4optCost(i1, i2, j1, j2);
        //            if (changedCost - valCost > -MY_EPSILON) {
        //                changedCost = valCost;
        //                bestI1 = i1;
        //                bestI2 = i2;
        //                bestJ1 = j1;
        //                bestJ2 = j2;
        //            }
        //        }
        //    }

        //    //j1+1 = j2
        //    for (i2 = i1 + 1; i2 < n; ++i2)
        //        for (j1 = i2 + 1; j1 < n; ++j1) {
        //            j2 = j1 + 1;
        //            if (j2 == n)continue;
        //            valCost = cal4optCost(i1, i2, j1, j2);
        //            if (changedCost - valCost > -MY_EPSILON) {
        //                changedCost = valCost;
        //                bestI1 = i1;
        //                bestI2 = i2;
        //                bestJ1 = j1;
        //                bestJ2 = j2;
        //            }
        //        }
        //}

        //O(n^2)
        for (i1 = 0; i1 < n-3; ++i1) {
            //i1+1=i2 i2+1=j1
            i2 = i1 + 1;
            j1 = i2 + 1;
            for (int j2 = j1 + 2; j2 < n; ++j2) {
                valCost = cal4optCost(i1, i2, j1, j2);
                if (changedCost > valCost) {
                    changedCost = valCost;
                    bestI1 = i1;
                    bestI2 = i2;
                    bestJ1 = j1;
                    bestJ2 = j2;
                }
            }

            //i1+1=i2 j1+1=j2
            for (j1 = i2 + 2; j1 < n - 1; ++j1) {
                j2 = j1 + 1;
                valCost = cal4optCost(i1, i2, j1, j2);
                if (changedCost > valCost) {
                    changedCost = valCost;
                    bestI1 = i1;
                    bestI2 = i2;
                    bestJ1 = j1;
                    bestJ2 = j2;
                }
            }        

            //i2+1=j1 j1+1=j2
            for (i2 = i1 + 2; i2 < n - 2; ++i2) {
                j1 = i2 + 1;
                j2 = j1 + 1;
                valCost = cal4optCost(i1, i2, j1, j2);
                if (changedCost > valCost) {
                    changedCost = valCost;
                    bestI1 = i1;
                    bestI2 = i2;
                    bestJ1 = j1;
                    bestJ2 = j2;
                }
            }
            //i1+1=i2 i2+1=j1 j1+1=j2
            i2 = i1 + 1; j1 = i2 + 1; j2 = j1 + 1;
            valCost = cal4optCost(i1, i2, j1, j2);            
            if (changedCost > valCost) {                
                changedCost = valCost;
                bestI1 = i1;
                bestI2 = i2;
                bestJ1 = j1;
                bestJ2 = j2;
            }
        }
        valI1 = bestI1; valI2 = bestI2; valJ1 = bestJ1; valJ2 = bestJ2;
        //cout << bestI1 << " " << bestI2 << " " << bestJ1 << " " << bestJ2 << "\n";
        //return cost + changedCost;
        return changedCost;
    }

    //LNS-based pertubation
    void pertubation(bool isLarge) {        
        if (n <= 5)return;
        //remove
        int typeRmv = pr->Rng.getNumInRan(0, 1);
       /* pr->minRmv = max(1, (int)(n * pr->minRateSmallRmv));
        pr->maxRmv = max(1, (int)(n * pr->maxRateSmallRmv));
        if (isLarge) {
            pr->minRmv = max(1, (int)(n * pr->minRateBigRmv));
            pr->maxRmv = max(1, (int)(n * pr->maxRateBigRmv));
        }*/
        int nbRmv = pr->Rng.getNumInRan(pr->minRmv, pr->maxRmv);                
        switch (typeRmv)
        {
        case 0:
            randomRmv(nbRmv);            
            break;
        case 1:
            worstRmv(nbRmv, pr->worstDeg);            
            break;
        case 2:
            blockRmv(nbRmv);
            break;
        default:
            break;
        }
        //insertion        
       int typeIns = pr->Rng.getNumInRan(0, 1);              
       lastInsType = typeIns;       
        switch (typeIns)
        {
        case 0:            
            cheapestIns();
            break;
        case 1:
            randomIns();            
            break;
        case 2: 
            nearestIns();
            break;
        default:
            break;
        }
    }      
    //deconstructor:
    ~Solution() {
        /*for (int k = 0; k <= m; ++k) {
            delete[] F[k];
            delete[] pred[k];
        }*/        
        delete[] seqSet;        
        giantT.clear();        
        ordNodeLs.clear();
        //for (int i = 0; i < n + 2 * pr->numVeh + 1; ++i)delete nodes[i];
        //for (int i = 1; i <= m; ++i)delete setR[i];
        //delete pr;
    }
    
};