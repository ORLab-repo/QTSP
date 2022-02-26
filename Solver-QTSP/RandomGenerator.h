#pragma once
#include <cstdlib>
#include<random>

class RandomGenerator {
public:
    unsigned seed;
    std::mt19937 generator;// using for gen num in range [0,1]    
    std::uniform_real_distribution<double> distribution;
    std::uniform_real_distribution<double> distribution1;

    void config(unsigned int sd) {
        seed = sd;
        generator.seed(seed);
        std::uniform_real_distribution<double> valDis(0.0, 1.0);
        distribution = valDis;
        std::uniform_real_distribution<double> valDis1(0.0, 1.0);
        distribution1 = valDis1;
    }

    //Gen real number in range [0,1]    
    double genRealInRang01() {
        return distribution(generator);
    }

    //using for mutation_reversal    
    double genRealInRang01_muta() {
        return distribution1(generator);
    }

    int getIntRand() {
        return generator();
    }

    //gen number in range [u,v]
    int getNumInRan(int u, int v) {
        return (generator() % (v - u + 1)) + u;
    }
    /*
    inline double randomdN(int n)
    {
        return (double)rand() / ((double)RAND_MAX / n);
    }
    inline int randomN(int n)
    {
        return (rand() % n) + 1;
    }
    inline int randomR(int u, int v)
    {
        return (rand() % (v - u + 1)) + u;
    }*/
};
