#include "Ga.h"

GA::GA()
{
}

GA::~GA()
{
}


void GA::updateBiasedFitnesses()
{
   /* vector<DI> ranking;
    for (int i = 1; i <= curNPop; i++)
    {
        ranking.push_back({ -pop[i]->averageBrokenPairsDistanceClosest(nClose), i});
    }
    sort(ranking.begin(), ranking.end());    
    
    for (int i = 0; i < curNPop; ++i) {
        double divRank = (double)(i + 1) / (curNPop);
        double fitRank = (double)ranking[i].second / (curNPop);
        if (curNPop <= nElite)pop[ranking[i].second]->biasedFitness = fitRank;
        else pop[ranking[i].second]->biasedFitness = fitRank + (1.0 - (double)nElite / curNPop) * divRank;
    }*/
    for (int i = 1; i <= curNPop; ++i) {
        //double fitRank = (double)(i) / (curNPop);        
        double fitRank = i;
        pop[i]->biasedFitness = fitRank;
    }
    FindAdapt();
}

void GA::removeWorstIndv()
{
    updateBiasedFitnesses();
    Solution* worstIndividual = NULL;
    int worstIndividualPosition = -1;
    bool isWorstIndividualClone = false;
    double worstIndividualBiasedFitness = -1.e30;
    for (int i = 1; i <= curNPop; i++)
    {
        bool isClone = (pop[i]->averageBrokenPairsDistanceClosest(1) < MY_EPSILON); // A distance equal to 0 indicates that a clone exists
        if ((isClone && !isWorstIndividualClone) || (isClone == isWorstIndividualClone && pop[i]->biasedFitness > worstIndividualBiasedFitness))
        {
            worstIndividualBiasedFitness = pop[i]->biasedFitness;
            isWorstIndividualClone = isClone;
            worstIndividualPosition = i;
            worstIndividual = pop[i];
        }
    }
    for (int i = worstIndividualPosition; i < curNPop; ++i) {
        swap(pop[i], pop[i + 1]);
    }
    curNPop--;
    for (int i = 1; i <= curNPop; ++i)pop[i]->removeProximity(worstIndividual);
}

bool GA::checkIdSol(Solution* u)
{
    for (int i = 1; i <= n; ++i)ddID[i] = 0;
    for (int i = 1; i <= n; ++i)ddID[u->giantT[i]] = 1;
    for (int i = 1; i <= n; ++i)if (ddID[i] == 0)return false;
    return true;
}

int GA::broken_pairs(Solution* u, Solution* v)
{
    //computing broken pair dis of u with v
    int* nxt = new int[n + 1];
    for (int i = 0; i < n - 1; ++i)nxt[u->giantT[i]] = u->giantT[i + 1];
    nxt[u->giantT[n - 1]] = -1; 
    int numEqualEdge = 0;
    for (int i = 0; i < n - 1; ++i)if (v->giantT[i + 1] == nxt[v->giantT[i]])numEqualEdge++;
    delete[] nxt;
    return n - 1 - numEqualEdge;
}

bool GA::CheckEqual(Solution* u, Solution* v)
{    
    //if (abs(u->cost - v->cost) <= omega) return true;
    /*for (int i = 1; i <= n; ++i)
        if(u->giantT[i] != v->giantT[i])return false;*/
    return false;
    /*return(abs(u->cost - v->cost) <= omega)
        || (broken_pairs(u, v) < threshold);*/
}

void GA::equalSol(Solution* u, Solution* v)
{
    u->cost = v->cost;
    u->biasedFitness = v->biasedFitness;
    for (int i = 1; i <= n; ++i) {
        u->giantT[i] = v->giantT[i];
        u->predecessors[i] = v->predecessors[i];
        u->successors[i] = v->successors[i];
    }
}

//For Roulette Wheel Selection
void GA::FindAdapt()
{    
    for (int i = 1; i <= curNPop; ++i)
        adapt[i] = curNPop + 1 - pop[i]->biasedFitness;
        //adapt[i] = (double)(curNPop + 1 - pop[i]->biasedFitness) / curNPop;        
        //adapt[i] = 1.0 / pop[i]->biasedFitness;
    sumAdapt = 0;
    for (int i = 1; i <= curNPop; ++i)
        sumAdapt += adapt[i];
}

int GA::getChild()
{
    /// binary tournament
    //int u = pr->Rng.getNumInRan(1, curNPop);
    //int v = pr->Rng.getNumInRan(1, curNPop);
    ///*return (pop[u]->cost < pop[v]->cost) ? u : v;*/
    //return (pop[u]->biasedFitness < pop[v]->biasedFitness) ? u : v;
    ///routle wheel selection
    double percent = pr->Rng.genRealInRang01();
    double sum1 = 0;
    for (int i = 1; i <= curNPop; ++i) {
        sum1 += adapt[i];
        if (sumAdapt * percent <= sum1)return i;
    }
    return 1;
}

void GA::choose(int& u, int& v)
{
    updateBiasedFitnesses();    
    /// new selection
    u = getChild();
    v = getChild();
    /*do {
        v = getChild();
    } while (u == v);*/
    /*double avg=sumAdapt/2;
    percent=randomdN(0,1);
    double sum1=percent*avg;
    u=1;sumAdapt=0;*/   
}

void GA::uni(Solution* u, Solution* v, Solution* u1, Solution* v1)
{
    deque<int> q;
    q.clear();
    int* idu = new int[n + 1];
    int* idv = new int[n + 1];
    for (int i = 1; i <= n; ++i) {
        idu[i] = u->giantT[i];
        idv[i] = v->giantT[i];
    }
    /*for(int i=0;i<n;++i)out<<u.id[i]<<" ";
    out<<endl;
    for(int i=0;i<n;++i)out<<v.id[i]<<" ";
    out<<endl;*/
    //int rn=randomR(2*n/5,3*n/5);
    //int be=randomR(0,n-rn);
    /// test:
    //int en=be+rn-1;
    int be = pr->Rng.getNumInRan(1, n);
    //int en = pr->Rng.getNumInRan(be, n);
    int en = pr->Rng.getNumInRan(1, n);    
    while (be == en)
    {
        en = pr->Rng.getNumInRan(1, n);
    }
    int* dd = new int[n + 1];
    int vt;
    //born u:
    for (int i = 0; i <= n; ++i)
        dd[i] = 0;   
    int idI = be;
    /*for (int i = be; i <= en; ++i) {        
        u1->giantT[i] = idu[i];
        dd[idu[i]] = 1;
    }*/
    while (true)
    {        
        if ((en + 1) % n == idI % n)break;
        u1->giantT[idI] = idu[idI];
        dd[idu[idI]] = 1;
        idI++;        
        if (idI > n) idI = 1;
    }    
    /*for (int i = 1; i <= n; ++i)
        if (dd[idv[i]] == 0)
            q.push_back(idv[i]);
    vt = en + 1;
    while (!q.empty())
    {
        if (vt > n)
            vt = 1;
        u1->giantT[vt] = q.front();
        q.pop_front();
        vt++;
    }*/    
    int curId = en;
    for (int i = 1; i <= n; ++i) {                
        curId++;
        if (curId > n)curId = 1;
        if (dd[idv[curId]] == 0) {
            u1->giantT[idI] = idv[curId];
            idI = (idI + 1) % n;
            if (idI == 0)idI = n;
        }
    }
    //born v:
    /*for (int i = 0; i <= n; ++i)
        dd[i] = 0;
    for (int i = be; i <= en; ++i) {
        v1->giantT[i] = idv[i];
        dd[idv[i]] = 1;
    }
    for (int i = 1; i <= n; ++i)
        if (dd[idu[i]] == 0)
            q.push_back(idu[i]);
    vt = en + 1;
    while (!q.empty())
    {
        if (vt > n)
            vt = 1;
        v1->giantT[vt] = q.front();
        q.pop_front();
        vt++;
    }*/  

    if (pr->Rng.genRealInRang01_muta() > 1-pM) {
        for (int i = 1; i <= nMut; ++i)u1->exchange();
    }
    u1->calCost();
    u1->updateObj();
    //v1->Split(); v1->updateTotal();
    /*if(rand()%2==0){
        u1.updateObj();v1.updateObj();
    }else{
        u1.newSplit();v1.newSplit();
    }*/
    delete[] dd;
    delete[] idu;
    delete[] idv;
}

/*diversity constribution*/


void GA::insertNew(Solution* u)
{
    int maxObj = oo;
    int posIns = curNPop + 1;    

    if (curNPop == 0) {
        equalSol(pop[++curNPop], u);     
        pop[curNPop]->indivsPerProximity.clear();
        return;
    }    
    
    if (u->cost < pop[1]->cost)posIns = 1;
    else {
        /*for (int i = 1; i <= curNPop; ++i)
        {
            if (CheckEqual(u, pop[i]))
            {                                
                return;
            }
        }*/

        for (int i = 1; i <= curNPop; ++i)
        {
            if (pop[i]->cost > u->cost) {
                posIns = i;
                break;
            }
        }
    }        
    //pop[curNPop + 1] = u;
    equalSol(pop[curNPop + 1], u);
    pop[curNPop + 1]->indivsPerProximity.clear();
    for (int i = 1; i <= curNPop; ++i) {
        double disBrokPair = pop[curNPop + 1]->brokenPairsDistance(pop[i]);
        pop[i]->indivsPerProximity.insert({ disBrokPair, pop[curNPop + 1] });
        pop[curNPop + 1]->indivsPerProximity.insert({disBrokPair, pop[i]});
    }

    for (int i = curNPop; i >= posIns; --i) {
        swap(pop[i], pop[i + 1]);
    }
    curNPop++;
    if (curNPop > nPop + delta) {
        DelPopu();
    }
}

void GA::DelPopu()
{
    while (curNPop > nPop) {
        removeWorstIndv();
    }
}

void GA::InitPopu(bool isEdu = true)
{        
    cout << "start init: " << curNPop << "\n";
    while (curNPop != nPop) {       
        cout << curNPop << "\n";
        valPop->genGiantT();
        valPop->calCost();
        if(isEdu)valPop->updateObj();
        //cout<<valPop->cost<<endl;
        insertNew(valPop);
    }
}

void GA::DiversifyPopu(Solution* bestSol)
{
    int newPop = nPop / 3;    
    while (curNPop > newPop) {
        removeWorstIndv();
    }
    InitPopu(false);
    if (bestSol->cost > pop[1]->cost) {
        //pop[i].printSol();
        equalSol(bestSol, pop[1]);
    }
}

void GA::findGasSol(int maxNumGas)
{    
    int idFa, idMo;
    Solution* bestSol = new Solution(pr);
    clock_t be = clock();
    // generate population with 50 sols
    curNPop = 0;
    InitPopu(false);
    /*for(int i=1;i<=25;++i){
        for(int j=0;j<n;++j)out<<pop[i].id[j]<<" ";
        out<<endl;
    }*/
    // hybrid 50 times:
    bestSol->cost = oo;
    equalSol(bestSol, pop[1]);
    //maxNumGas=500;
    numNotCha = 0;
    threshold = (n - 1) / 2;
    Solution* child1 = new Solution(pr);
    Solution* child2 = new Solution(pr);        
    for (int numga = 1;; ++numga)
    {
        /*cout << "time: " << (double)(clock() - be) / CLOCKS_PER_SEC << "\n";*/
        /*if (numNotCha == 0)pM = pMinMut;
        else if (numNotCha % ItMut == 0)pM = min(pM + 0.1, pMaxMut);*/
        numNotCha++;
        cout<<numga<<":"<<endl;
        //cout << "name ins: " << pr->nameIns << "\n";
        /*cout << numNotCha << " " << numga << "{" << endl;*/
        //out<<numga<<"{\n";

        //exit(0);
        // calc adaptation
        //FindAdapt();
        // selection
        choose(idFa, idMo);
        /*cout<<"parent"<<endl;
        cout<<pop[idFa].obj<<" "<<pop[idMo].obj<<endl;*/
        /*cout<<idFa<<" "<<idMo<<endl;
        for(int i=0;i<n;++i)cout<<pop[idFa].id[i]<<" ";
        cout<<endl;
        for(int i=0;i<n;++i)cout<<pop[idMo].id[i]<<" ";
        cout<<endl;*/
        // hybrid
        uni(pop[idFa], pop[idMo], child1, child2);
        /*cout<<"children"<<endl;
        cout<<child1.obj<<" "<<child2.obj<<endl;*/
        /*for(int i=0;i<n;++i)cout<<child1.id[i]<<" ";
        cout<<endl;
        for(int i=0;i<n;++i)cout<<child2.id[i]<<" ";
        cout<<endl;*/
        //insert child:
        insertNew(child1);        
        /*insertNew(child2);
        if (nPop1 == nPop + delta)
            DelPopu();*/
        if (bestSol->cost > pop[1]->cost) {
            numNotCha = 0;
            equalSol(bestSol, pop[1]);
            cout << "iteration: " << numga << "\n";
            cout << "new best: " << bestSol->cost << "\n";
            pr->fileOut << "itreation: " << numga << "\n";
            pr->fileOut << "new best: " << bestSol->cost << "\n";            
            //pr->fileOut << (double)(clock() - be) / CLOCKS_PER_SEC << "\n";
        }
        // if bestSol don't change 100 times change 25 worst sols by 25 new sols
        //if(numNotCha>=200){numNotCha=0;addPopu();}
        //if (curNumRou >= maxNumRou) {      
        /*
        if (numNotCha == (int)(0.4*ItNI))
        {
            //threshold = min(threshold - 1, 1);
            int oldBestObj = bestSol->cost;
            DiversifyPopu(bestSol);
            if (bestSol->cost != oldBestObj) {
                numNotCha = 0;
                pr->fileOut << (double)(clock() - be) / CLOCKS_PER_SEC << "\n";
            }
        }*/       
        //cout<<"best obj:\n";cout<<bestSol.obj<<endl<<endl;
        if (numNotCha == ItNI || (double)(clock() - be) / CLOCKS_PER_SEC > pr->TL) { //original termination        
            bestSol->calCost();
            /*bestSol->ckSol();*/
            /*bestSol.printSol();
            cout << id_test << " " << bestSol.obj << " " << (double)(clock() - be) / CLOCKS_PER_SEC << endl;
            fl << bestSol.obj << " " << (double)(clock() - be) / CLOCKS_PER_SEC << "\n";*/
            bestCost = bestSol->cost;
            //cout << bestSol->cost;
            cout << "Cost: " << bestSol->cost << "\n";
            pr->fileOut << "Cost: " << bestSol->cost << "\n";            
            pr->fileOut << "giantTour: ";
            for (int i = 1; i <= n; ++i)pr->fileOut << bestSol->giantT[i] << ", ";
            pr->fileOut << "\n";
            //pr->fileOut <<(double)(clock() - be) / CLOCKS_PER_SEC << "\n";
            /*pr->fileOut << "num_iterations: " << numga << "\n";*/
            break;
        }
        //if((double)(clock()-be)/CLOCKS_PER_SEC>=600)break;
        /*cout << curNumRou << "\n";
        cout << "SCP sol: " << scpSol << "\n";
        cout << bestSol->cost << "}" << endl;*/
    }
}
