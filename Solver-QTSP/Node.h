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
    //vector <SeqData*> seqi_j; // data for (i,j) with j > i
    //vector <SeqData*> seqj_i; // data for (j,i) (for the same subsequence as i_j, but reversed)    
    int demand;    
    double culCost;    
    vector<int> moves; //init based on correlation measure for locations (still contain the index of customer)    
    int idxClient = -1;
    /*int idxLoc = -1;*/

    Node() {
        isDepot = false;        
    }

    bool ckNearDepot() {
        return (pred->isDepot|| suc->isDepot);
    }

    ~Node() {        
        moves.clear();        
    };
};
