#include <iostream>
#include <cmath>
#include <vector>
#include <memory>

#include "chaosmap.hpp"
#include "point.hpp"
#include "swarm.hpp"
#include "function.hpp"
#include "utils.hpp"

int main()
{
    int nN, nC;
    int max_iterations = 7;

    std::cout << "Enter the number of normal particles." << std::endl;
    std::cin >> nN;
    std::cout << "\n"
              << "Enter the number of chaotic particles." << std::endl;
    std::cin >> nC;

    using T = float;
    constexpr int dim = 2;

    std::unique_ptr<Swarm::Function<T, dim>> fitness = std::make_unique<Swarm::DropwaveFunction<T, dim>>();

    Swarm::Point<T, dim> map_a(-1.f);
    Swarm::Point<T, dim> map_b(1.f);
    Swarm::MapFunction<T, dim> f = [](const Swarm::Point<T, dim> &p, int k)
    { return (p.arccos() * k).cos(); };

    Swarm::Point<T, dim> a(-1000.f);
    Swarm::Point<T, dim> b(+1000.f);

    Swarm::ChaosMap<T, dim> chaosMap(f, map_a, map_b);

    Swarm::Swarm<T, dim> swarm(fitness, a, b, max_iterations);

    for (int i = 0; i < nN; ++i)
    {
        swarm.addParticle(std::make_unique<Swarm::NormalParticle<T, dim>>());
    }

    for (int i = 0; i < nC; ++i)
    {
        swarm.addParticle(std::make_unique<Swarm::ChaoticParticle<T, dim>>(chaosMap));
    }

    for (int i = 0; i < max_iterations; ++i)
    {
        std::cout << "Iteration: " << i << std::endl;
        swarm.updateEveryone();
        std::cout << "Best value = " << swarm.getGlobalBestValue() << std::endl;
        std::cout << "Best position = " << swarm.getGlobalBest() << std::endl;
    }

    std::cout << "Best value = " << swarm.getGlobalBestValue() << std::endl;
    std::cout << "Best position = " << swarm.getGlobalBest() << std::endl;
}
