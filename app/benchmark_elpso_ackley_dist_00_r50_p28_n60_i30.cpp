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
    int max_iterations = 30;

    omp_set_num_threads(threads);

    Swarm::Point<T, dim> a(-30.f);
    Swarm::Point<T, dim> b(+30.f);

    Swarm::Point<T, dim> initial_best(0.5f);

    Swarm::Optimizers::GENETICOptimizer<T, dim> swarm(fitness, initial_best, a, b, nN, nC, chaosMap, max_iterations, "/work/u10768804/benchmark/positions_benchmark_elpso_ackley_dist_00_r50_p28_n60_i30.txt");

    Swarm::Point<T, dim> center(0.f);
    T radius = 5.f;
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

    nlohmann::json benchmark;

    constexpr int cpus = 28;

    size_t j = 0;

    int particles = 60;

    std::shared_ptr<Swarm::ObjectiveFunction<T, dim>> fitness = std::make_shared<Swarm::ObjectiveFunctions::AckleyFunction<T, dim>>();

    std::cout << "\n============================\n"
              << "Benchmarking with " << (cpus) << " threads and " << particles << " particles per type.\n";

    auto elapsed = measure(cpus, particles, std::move(fitness), Swarm::ChaosFactory::Chebyshev<T, dim>());

    std::cout << "Time elapsed: " << elapsed.count() << " seconds.\n"
              << "============================\n";

    benchmark[std::to_string((cpus - 1) * 4 + j + 1)] = {
        {"threads", cpus},
        {"particles_per_type", particles},
        {"time_seconds", elapsed.count()}};

    std::ofstream benchmark_file("/work/u10768804/benchmark/benchmark_elpso_ackley_dist_00_r50_p28_n60_i30.json");
    benchmark_file << benchmark.dump(4);
    benchmark_file.close();

    return 0;
}