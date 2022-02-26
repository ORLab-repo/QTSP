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
	int cost;// will set to oo if infeasible sub-sequence	
	bool F;//check feasibility;
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
		F = true;
	}

	/**/
	void concatOneAfter(SeqData* seq, int idxLoc) {
		if (seq == NULL) init(idxLoc);
		F = seq->F;

		int t_12 = pr->times[seq->lastnode][idxLoc];		
		//cost = min(seq->cost + pr->costs[seq->lastnode][idxLoc], int(oo));				
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
		F = seq->F;

		int t_12 = pr->times[idxLoc][seq->firstnode];		
		//cost = min(seq->cost + pr->costs[idxLoc][seq->firstnode], int(oo));		
		firstnode = idxLoc;
		lastnode = seq->lastnode;
		if (seq->beforeLaNode == -1) {
			beforeLaNode = idxLoc;
		}
		else beforeLaNode = seq->beforeLaNode;
		afterFiNode = seq->firstnode;
	}
	/**/
	int evaluation(vector<SeqData*> seqs) {
		if (seqs.front() == NULL)seqs.erase(seqs.begin());//remove null sequence
		if (seqs.back() == NULL)seqs.pop_back();//remove null sequence
		int costR = seqs[0]->cost;
		bool totalF = seqs[0]->F;
		int u = -1, v = -1;
		for (int i = 0; i < seqs.size() - 1; ++i) {
			u = seqs[i]->lastnode;
			v = seqs[i + 1]->firstnode;
			//costR += (pr->costs[u][v] + seqs[i + 1]->cost);
			//totalF = (totalF & seqs[i + 1]->F) & (totalE + pr->times[u][v] <= seqs[i + 1]->L) & (loadR <= pr->Q) & (costR < oo);
			//if (!totalF)return oo;// only use for checking feasible solution			
			return costR;
		}
	}

private:
};