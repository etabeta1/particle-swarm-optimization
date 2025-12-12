#include <iostream>
#include <cmath>
#include <vector>
#include <memory>
#include <chrono>

#include "chaosmap.hpp"
#include "point.hpp"
#include "swarm.hpp"
#include "function.hpp"
#include "utils.hpp"

int main()
{
    int nN = 10000, nC = 10000;
    int max_iterations = 10000;

    // std::cout << "Enter the number of iterations" << std::endl;
    // std::cin >> max_iterations;
    // std::cout << "Enter the number of normal particles." << std::endl;
    // std::cin >> nN;
    // std::cout << "\n"
    //           << "Enter the number of chaotic particles." << std::endl;
    // std::cin >> nC;

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

    auto start = std::chrono::high_resolution_clock::now();
    swarm.run();

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "\n---------FINAL_RESULTS----------\n"
              << std::endl;
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "target_function CPU time: "
              << elapsed.count() << " s\n";
    auto best = swarm.getGlobalBest();
    std::cout << "Best value = " << best.value << std::endl;
    std::cout << "Best position = " << best.point << std::endl;
}
