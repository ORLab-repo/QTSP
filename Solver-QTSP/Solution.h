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
    bool isFinished;//used for LS
    string typeIns;
    vector<SeqData*> mySeq;
    SeqData* seqSet;
    SeqData* seqDep;
    SeqData* valSeq[12];    
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
        n = pr->numLoc - 1;//not cointain depot
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
        int maxSizeSeq = 7 * n + 1;
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
    }

    void reinitNegiborSet() {
        set<DI> sDis;
        int u, v, vSuc;
        for (int i = 1; i <= n; ++i) {
            sDis.clear();
            u = nodes[i]->idxClient;
            for (int j = 1; j <= n + 1; j++)if (i != j) {
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

    //swap two nodes u and v
    void swapNode(Node *u, Node* v){
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
    //relocate u after v
    bool move1() {
        cout << "move1\n";
        //nodeV in this case can be depot (arrival depot)
        if (uSuc == nodeV)return false;        
        int posUInSol = nodeU->posInSol;
        int posVInSol = nodeV->posInSol;
        if (nodeV->isDepot)posVInSol = nodeV->pred->posInSol + 1;
        if (nodeU->posInSol < nodeV->posInSol) {
            valSeq[1] = uPred->seq0_i;
            constructSeqData(uSuc, vPred, valSeq[2]);            
            valSeq[3] = nodeU->seqi_i;
            valSeq[4] = nodeV->seqi_n;
        }
        else {
            valSeq[1] = vPred->seq0_i;
            valSeq[2] = nodeU->seqi_i;
            constructSeqData(nodeV, uPred, valSeq[3]);
            valSeq[4] = uSuc->seqi_n;
        }
        mySeq.clear();
        for (int i = 1; i <= 4; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        insertNode(nodeU, nodeV);
        updateInfo();
        isFinished = false;
        cost = newCost;        
        return true;
    }

    //relocate u Suc after V
    bool move2() {
        cout << "move2\n";
        if (nodeU == vSuc || nodeV == uSuc || uSuc->suc == nodeV || uSuc->isDepot)return false;
        int posUInSol = nodeU->posInSol;
        int posVInSol = nodeV->posInSol;
        if (nodeV->isDepot)posVInSol = nodeV->pred->posInSol + 1;
        if (nodeU->posInSol < nodeV->posInSol) {
            valSeq[1] = uPred->seq0_i;
            constructSeqData(uSuc->suc, vPred, valSeq[2]);
            valSeq[3] = nodeU->seqi_j;
            valSeq[4] = nodeV->seqi_n;
        }
        else {
            valSeq[1] = vPred->seq0_i;
            valSeq[2] = nodeU->seqi_j;
            constructSeqData(nodeV, uPred, valSeq[3]);
            valSeq[4] = uSuc->suc->seqi_n;
        }
        mySeq.clear();
        for (int i = 1; i <= 4; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        insertNode(uSuc, nodeV);
        insertNode(nodeU, uSuc);
        updateInfo();
        isFinished = false;
        cost = newCost;
        return true;
    }
    //relocate sucU u after V
    bool move3() {
        cout << "move3\n";
        if (nodeU == vSuc || nodeV == uSuc || uSuc->suc == nodeV || uSuc->isDepot)return false;
        int posUInSol = nodeU->posInSol;
        int posVInSol = nodeV->posInSol;
        if (nodeV->isDepot)posVInSol = nodeV->pred->posInSol + 1;
        if (nodeU->posInSol < nodeV->posInSol) {
            valSeq[1] = uPred->seq0_i;
            constructSeqData(uSuc->suc, vPred, valSeq[2]);
            valSeq[3] = nodeU->seqj_i;
            valSeq[4] = nodeV->seqi_n;
        }
        else {
            valSeq[1] = vPred->seq0_i;
            valSeq[2] = nodeU->seqj_i;
            constructSeqData(nodeV, uPred, valSeq[3]);
            valSeq[4] = uSuc->suc->seqi_n;
        }
        mySeq.clear();
        for (int i = 1; i <= 4; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        insertNode(nodeU, nodeV);
        insertNode(uSuc, nodeU);
        updateInfo();
        isFinished = false;
        cost = newCost;
        return true;
    }
    //Swaps moves:
    //swap u and v (u->pos < v->pos)
    bool move4() {
        cout << "move4\n";
        if (nodeU == vPred || nodeU == vSuc)return false;        
        valSeq[1] = uPred->seq0_i;
        valSeq[2] = nodeV->seqi_i;
        constructSeqData(uSuc, vPred, valSeq[3]);
        valSeq[4] = nodeU->seqi_i;
        valSeq[5] = vSuc->seqi_n;
        mySeq.clear();
        for (int i = 1; i <= 5; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        swapNode(nodeU, nodeV);
        updateInfo();
        isFinished = false;
        cost = newCost;
        return true;
    }

    //swap (u, uSuc) and v
    bool move5() {
        cout << "move5\n";
        if (nodeU == vPred || uSuc == vPred || nodeU == vSuc || uSuc->isDepot)return false;
        int posUInSol = nodeU->posInSol;
        int posVInSol = nodeV->posInSol;
        if (nodeV->isDepot)posVInSol = nodeV->pred->posInSol + 1;
        if (nodeU->posInSol < nodeV->posInSol) {
            valSeq[1] = uPred->seq0_i;
            valSeq[2] = nodeV->seqi_i;
            constructSeqData(uSuc->suc, vPred, valSeq[3]);
            valSeq[4] = nodeU->seqi_j;
            valSeq[5] = vSuc->seqi_n;
        }
        else {
            valSeq[1] = vPred->seq0_i;
            valSeq[2] = nodeU->seqi_j;
            constructSeqData(vSuc, uPred, valSeq[3]);
            valSeq[4] = nodeV->seqi_i;
            valSeq[5] = uSuc->suc->seqi_n;
        }
        mySeq.clear();
        for (int i = 1; i <= 5; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        swapNode(nodeU, nodeV);
        insertNode(uSuc, nodeU);
        updateInfo();
        isFinished = false;
        cost = newCost;
        return true;
    }

    //swap (u, Suc) and (v, vSuc) (u->pos < v->pos)
    bool move6() {
        cout << "move6\n";
        if (uSuc->isDepot || vSuc->isDepot || vSuc == uPred || nodeU == vSuc || uSuc == nodeV || nodeV == uSuc->suc)return false;
        valSeq[1] = uPred->seq0_i;
        valSeq[2] = nodeV->seqi_j;
        constructSeqData(uSuc->suc, vPred, valSeq[3]);
        valSeq[4] = nodeU->seqi_j;
        valSeq[5] = vSuc->suc->seqi_n;
        mySeq.clear();
        for (int i = 1; i <= 5; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        swapNode(nodeU, nodeV);
        swapNode(uSuc, vSuc);
        updateInfo();
        isFinished = false;
        cost = newCost;
        return true;
    }
    
    //2-opt
    bool move7() {
        cout << "move7\n";
        if (nodeU->posInSol > nodeV->posInSol)return false;
        if (uSuc == nodeV)return false;
        valSeq[1] = nodeU->seq0_i;
        constructSeqData(nodeV, uSuc, valSeq[2]);
        valSeq[3] = vSuc->seqi_n;
        mySeq.clear();
        for (int i = 1; i <= 3; ++i)mySeq.push_back(valSeq[i]);
        double newCost = seqDep->evaluation(mySeq);
        if (newCost - cost > -MY_EPSILON)return false;
        //update route structure:
        Node* startNode = nodeU;
        Node* endNode = nodeV;
        while (startNode != uSuc)
        {
            Node* tempNode = endNode->pred;
            startNode->suc = endNode;
            endNode->pred = startNode;
            startNode = endNode;
            endNode = tempNode;
        }
        uSuc->suc = vSuc;
        vSuc->pred = uSuc;
        updateInfo();     
        isFinished = false;
        cost = newCost;        
        return true;
    }

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
        isFinished = false;
        while (!isFinished) {
            isFinished = true;
            for (int posU = 0; posU < ordNodeLs.size(); ++posU) {                
                nodeU = nodes[ordNodeLs[posU]];
                for (int posV = 0; posV < nodeU->moves.size(); ++posV) {
                    cout << posU << " " << posV << "\n";
                    nodeV = nodes[nodeU->moves[posV]];
                    setLocalValU();
                    setLocalValV();
                    if (move1())continue;
                    if (move2())continue;
                    if (move3())continue;
                    if (nodeV->isDepot) continue;
                    if (nodeUIdx < nodeVIdx &&move4())continue;
                    if (move5())continue;
                    if (nodeUIdx < nodeVIdx && move6())continue;    
                    if (move7())continue;
                }
            }
        }
    }
    
};