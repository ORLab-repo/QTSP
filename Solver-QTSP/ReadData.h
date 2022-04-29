#pragma once
#include "lib.h"
#include "Param.h"

using namespace std;

Param* read_Ins(string path, string type) {
    Param* pr = new Param();
    ifstream filein(path);
    int countLine = 0;
    vector<string> lineCont;// content of Line 
    pr->listLoc.clear();
    while (filein) {
        string line;
        getline(filein, line);        
        lineCont = Util::splitString(line, " ");
        if (lineCont[0] == "NAME:") {
            pr->nameIns = lineCont[1];
            continue;
        }
        if (lineCont[0] == "TYPE:") {
            continue;
        }
        if (lineCont[0] == "DIMENSION:") {
            pr->numLoc = stoi(lineCont[1]);
            continue;
        }
        if (lineCont[0] == "EDGE_WEIGHT_TYPE:") {
            continue;
        }
        //1 NODE
        if (lineCont[0] == "NODE_COORD_SECTION") {
            countLine++;
            continue;
        }

        if (lineCont[0] == "EOF")break;

        //get coordinate:
        if (countLine == 1) {            
            pr->listLoc.push_back(Location(stoi(lineCont[1]),
                stoi(lineCont[2])
            )
            );            
        }        
    }
    pr->maxRmv = min(pr->maxRmv, pr->numLoc / 3);
    cout << "finish reading\n";
    //cal distance:    
    double angle;
    pr->correlatedNodes = vector<vector<int> >(pr->numLoc + 1);
    for (int i = 1; i < pr->numLoc; ++i) {
        for (int j = 1; j <= pr->numLoc; ++j) if (j != i)pr->correlatedNodes[i].push_back(j);
        shuffle(pr->correlatedNodes[i].begin(), pr->correlatedNodes[i].end(), pr->Rng.generator);
    }
    for (int i = 0; i < pr->numLoc; ++i) {                
        pr->costs.push_back(vector<vector<double>>(pr->numLoc, vector<double> (pr->numLoc)));                       
        for (int j = 0; j < pr->numLoc; ++j) {            
            for (int k = 0; k < pr->numLoc; ++k) {        
                pr->costs[i][j][k] = 0.;               
                if (i == j || j == k || i == k) {
                    continue;
                }
                else {                    
                    angle = (pr->listLoc[j].x - pr->listLoc[i].x) * (pr->listLoc[k].x - pr->listLoc[j].x)
                        + (pr->listLoc[j].y - pr->listLoc[i].y) * (pr->listLoc[k].y - pr->listLoc[j].y);
                    angle /= sqrt(pow((pr->listLoc[j].x - pr->listLoc[i].x), 2) + pow((pr->listLoc[j].y - pr->listLoc[i].y), 2));
                    angle /= sqrt(pow((pr->listLoc[k].x - pr->listLoc[j].x), 2) + pow((pr->listLoc[k].y - pr->listLoc[j].y), 2));
                }
                angle = max(-1.0, min(1.0, angle));                
                angle = acos(angle);
                pr->costs[i][j][k] = angle;
                if (type == "ag")pr->costs[i][j][k] *= 1000;
                else if (type == "ag-dis")pr->costs[i][j][k] = 100 * (40 * pr->costs[i][j][k]
                    + 0.5 * (pr->listLoc[i].calDis(pr->listLoc[j]) + pr->listLoc[j].calDis(pr->listLoc[k])));
            }
        }
    }
    return pr;
}

