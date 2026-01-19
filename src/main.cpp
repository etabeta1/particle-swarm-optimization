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

std::chrono::duration<double> measure(int threads, int particles)
{
    using T = float;
    constexpr int dim = 2;

    int nN = particles, nC = particles;
    int max_iterations = 100;

    omp_set_num_threads(threads);
    std::unique_ptr<Swarm::ObjectiveFunction<T, dim>> fitness = std::make_unique<Swarm::DropwaveFunction<T, dim>>();

    Swarm::Point<T, dim> map_a(-1.f);
    Swarm::Point<T, dim> map_b(1.f);
    Swarm::MapFunction<T, dim> f = [](const Swarm::Point<T, dim> &p, int k)
    { return (p.arccos() * k).cos(); };

    Swarm::Point<T, dim> a(-1.f);
    Swarm::Point<T, dim> b(+1.f);
    Swarm::Point<T, dim> initial_best(.5f);

    Swarm::ChaosMap<T, dim> chaosMap(f, map_a, map_b);

    Swarm::Optimizers::CHOPSOOptimizer<T, dim> swarm(fitness, initial_best, a, b, nN, nC, chaosMap, max_iterations, true);

    auto start = std::chrono::high_resolution_clock::now();
    swarm.run();
    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "\n---------FINAL_RESULTS----------\n"
              << std::endl;
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Threads: " << threads << std::endl;
    std::cout << "target_function CPU time: "
              << elapsed.count() << " s\n";
    auto best = swarm.getGlobalBest();
    std::cout << "Best value = " << best.value << std::endl;
    std::cout << "Best position = " << best.point << std::endl;

    return elapsed;
}

int main()
{

    // std::cout << "Enter the number of iterations" << std::endl;
    // std::cin >> max_iterations;
    // std::cout << "Enter the number of normal particles." << std::endl;
    // std::cin >> nN;
    // std::cout << "\n"
    //           << "Enter the number of chaotic particles." << std::endl;
    // std::cin >> nC;

    measure(1, 100);

    // std::ofstream measurements;
    // measurements.open("time.txt");

    // measurements << "Threads,T10k,T40k,T160k" << std::endl;

    // for (int i = 1; i <= 28; i++)
    // {
    //     measurements << i << ",";

    //     std::chrono::duration<double> t10k = measure(i, 10000);
    //     std::chrono::duration<double> t20k = measure(i, 40000);
    //     std::chrono::duration<double> t40k = measure(i, 160000);

    //     measurements << t10k << "," << t20k << "," << t40k << std::endl;
    // }

    // measurements.close();
}
