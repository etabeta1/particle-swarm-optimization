#include <iostream>
#include <cmath>
#include <vector>
#include <memory>

#include "chaosmap.hpp"
#include "point.hpp"
#include "swarm.hpp"
#include "function.hpp"
#include "utils.hpp"

int main(){
    using namespace swarm;
    int nN, nC;
    int max_iterations=1000;

    std::cout<<"Enter the number of normal particles."<<std::endl;
    std::cin>> nN;
    std::cout<<"\n"<<"Enter the number of chaotic particles."<<std::endl;
    std::cin>>nC;

    using T = float;
    const int dim = 2;

    Function<T, dim> fitness;

    ChaosMap<T, float, dim> chaosMap;

    Swarm<T, dim> swarm(fitness);

    for (int i = 0; i < nN; ++i) {
        swarm.addParticle(std::make_unique<NormalParticle<T,dim>>());
    }

    for (int i = 0; i < nC; ++i) {
        swarm.addParticle(std::make_unique<ChaoticParticleParticle<T,dim>>());
    }

    for(int i = 0; i < max_iterations; ++i){
        swarm.findGlobalBest();
        swarm.updateEveryone();
    }

    std::cout << "Best value = " << swarm.getGlobalBestValue() << std::endl;
    std::cout << "Best position = " << swarm.getGlobalBest() << std::endl;

    
}
