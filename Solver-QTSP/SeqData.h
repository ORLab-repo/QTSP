#pragma once
#include "Param.h"

/*
* Based on paper "A unified solution framework for multi-attribute vehicle routing problems"
*/

class SeqData
{
public:
	//attributes:
	int firstnode;//index location
	int afterFiNode;// index after first node
	int lastnode;//index location
	int beforeLaNode;// index before last node
	double cost;// will set to oo if infeasible sub-sequence		
	Param* pr;

	//constructor
	SeqData() {};	
	SeqData(Param* _pr) {
		pr = _pr;
	};
	~SeqData() {
	};

	void showSeq() {
		cout << pr->listLoc[firstnode].idxClient << " " << pr->listLoc[lastnode].idxClient << "\n";
	}
	void showSeqLoc() {
		cout << firstnode << " " << lastnode << "\n";
	}
	//method:
	/*
	* Construct sequence containing only one node
	*/
	void init(int idxLoc) {
		firstnode = idxLoc;
		afterFiNode = -1;
		lastnode = idxLoc;
		beforeLaNode = -1;		
		cost = 0;				
	}

	void copy(SeqData* val) {
		firstnode = val->firstnode;
		lastnode = val->lastnode;
		beforeLaNode = val->beforeLaNode;
		afterFiNode = val->afterFiNode;
		cost = val->cost;
	}

	/**/
	void concatOneAfter(SeqData* seq, int idxLoc) {
		if (seq == NULL) init(idxLoc);		
		
		if (seq->beforeLaNode != -1) {
			cost = seq->cost + pr->costs[seq->beforeLaNode][seq->lastnode][idxLoc];
		}
		else cost = 0;
		lastnode = idxLoc;
		firstnode = seq->firstnode;
		if (seq->afterFiNode == -1) {
			afterFiNode = idxLoc;
		}
		else afterFiNode = seq->afterFiNode;
		beforeLaNode = seq->lastnode;
	}
	/**/
	void concatOneBefore(SeqData* seq, int idxLoc) {
		if (seq == NULL) init(idxLoc);		

		if (seq->afterFiNode != -1) {
			cost = seq->cost + pr->costs[idxLoc][seq->firstnode][seq->afterFiNode];
		}
		else cost = 0;
		firstnode = idxLoc;
		lastnode = seq->lastnode;
		if (seq->beforeLaNode == -1) {
			beforeLaNode = idxLoc;
		}
		else beforeLaNode = seq->beforeLaNode;
		afterFiNode = seq->firstnode;
	}
	/**/
	double evaluation(vector<SeqData*> seqs) {
		if (seqs.front() == NULL)seqs.erase(seqs.begin());//remove null sequence
		if (seqs.back() == NULL)seqs.pop_back();//remove null sequence
		int costR = seqs[0]->cost;		
		int u = -1, v = -1;
		int uPred = -1, vSuc = -1;
		for (int i = 0; i < seqs.size() - 1; ++i) {
			u = seqs[i]->lastnode;		
			uPred = seqs[i]->beforeLaNode;
			if (uPred == -1) {
				if (i != 0)uPred = seqs[i - 1]->lastnode;
				else uPred = -1;
			}
			v = seqs[i + 1]->firstnode;
			vSuc = seqs[i + 1]->afterFiNode;
			if (vSuc == -1) {
				if (i + 2 != seqs.size())vSuc = seqs[i + 2]->afterFiNode;
				else vSuc = -1;
			}
			costR += seqs[i + 1]->cost;			
			//totalF = (totalF & seqs[i + 1]->F) & (totalE + pr->times[u][v] <= seqs[i + 1]->L) & (loadR <= pr->Q) & (costR < oo);
			//if (!totalF)return oo;// only use for checking feasible solution			
			return costR;
		}
	}

private:
};