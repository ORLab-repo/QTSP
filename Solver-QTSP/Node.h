#pragma once
#include "SeqData.h"
#include <vector>

using namespace std;

class Node {
public:
    bool isDepot;    
    int posInSol; 
    Node* pred;
    Node* suc; 
    SeqData* seq0_i;
    SeqData* seqi_0;
    SeqData* seqi_n;
    SeqData* seqn_i;
    SeqData* seqi_i;
    SeqData* seqi_j;
    SeqData* seqj_i;
    //vector <SeqData*> seqi_j; // data for (i,j) with j > i
    //vector <SeqData*> seqj_i; // data for (j,i) (for the same subsequence as i_j, but reversed)    
    int demand;    
    /*double culCost;    */    
    int idxClient = -1;
    /*int idxLoc = -1;*/

    Node() {
        isDepot = false;        
    }

    bool ckNearDepot() {
        return (pred->isDepot|| suc->isDepot);
    }

    ~Node() {                
    };
};
