#include <chaosmap.hpp>
#include <point.hpp>
#include <swarm.hpp>
#include <function.hpp>
#include <utils.hpp>

#include <nlohmann/json.hpp>

#include <chrono>
#include <iostream>

std::chrono::duration<double> measure(int threads, int particles, std::shared_ptr<Swarm::ObjectiveFunction<float, 2>> fitness, const Swarm::ChaosMap<float, 2> &chaosMap)
{
    using T = float;
    constexpr int dim = 2;

    int nN = particles, nC = particles;
    int max_iterations = 1000;

    omp_set_num_threads(threads);

    Swarm::Point<T, dim> a(-5.f);
    Swarm::Point<T, dim> b(+5.f);

    Swarm::Point<T, dim> initial_best(0.2f);

    Swarm::Optimizers::GENETICOptimizer<T, dim> swarm(fitness, initial_best, a, b, nN, nC, chaosMap, max_iterations);

    Swarm::Point<T, dim> center(0.f);
    T radius = 0.5f;
    auto dist_constraint = Swarm::makeMaxDistanceConstraint<T, dim>(center, radius);
    swarm.addConstraint(dist_constraint);

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
    using T = float;
    constexpr int dim = 2;

    constexpr int cpus = 2;

    int particles = 200;

    std::shared_ptr<Swarm::ObjectiveFunction<T, dim>> fitness = std::make_shared<Swarm::ObjectiveFunctions::DropwaveFunction<T, dim>>();

    std::cout << "\n============================\n"
              << "Benchmarking with " << (cpus) << " threads and " << particles << " particles per type.\n";

    auto elapsed = measure(cpus, particles, std::move(fitness), Swarm::ChaosFactory::Chebyshev<T, dim>());

    std::cout << "Time elapsed: " << elapsed.count() << " seconds.\n"
              << "============================\n";

    return 0;
}