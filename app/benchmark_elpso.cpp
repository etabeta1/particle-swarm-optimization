#include <chaosmap.hpp>
#include <point.hpp>
#include <swarm.hpp>
#include <function.hpp>
#include <utils.hpp>

#include <nlohmann/json.hpp>

#include <chrono>
#include <iostream>

std::chrono::duration<double> measure(int threads, int particles, std::unique_ptr<Swarm::ObjectiveFunction<float, 2>> fitness, const Swarm::ChaosMap<float, 2> &chaosMap)
{
    using T = float;
    constexpr int dim = 2;

    int nN = particles, nC = particles;
    int max_iterations = 8192;

    omp_set_num_threads(threads);

    Swarm::Point<T, dim> a(-1.f);
    Swarm::Point<T, dim> b(+1.f);

    Swarm::Point<T, dim> initial_best(.5f);

    Swarm::Optimizers::GENETICOptimizer<T, dim> swarm(fitness, initial_best, a, b, nN, nC, chaosMap, max_iterations);

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

    for (size_t i = 0; i < cpus; i++)
    {
        for (size_t j = 0; j < 4; j++)
        {
            int particles = 16384 << j;

            std::unique_ptr<Swarm::ObjectiveFunction<T, dim>> fitness = std::make_unique<Swarm::ObjectiveFunctions::DropwaveFunction<T, dim>>();

            std::cout << "\n============================\n"
                      << "Benchmarking with " << (i + 1) << " threads and " << particles << " particles per type.\n";

            auto elapsed = measure(i + 1, particles, std::move(fitness), Swarm::ChaosFactory::Chebyshev<T, dim>());

            std::cout << "Time elapsed: " << elapsed.count() << " seconds.\n"
                      << "============================\n";

            benchmark[std::to_string(i * cpus + j + 1)] = {
                {"threads", i + 1},
                {"particles_per_type", particles},
                {"time_seconds", elapsed.count()}};
        }
    }

    std::ofstream benchmark_file("/work/u10822715/benchmark/benchmark_elpso.json");
    benchmark_file << benchmark.dump(4);
    benchmark_file.close();

    return 0;
}