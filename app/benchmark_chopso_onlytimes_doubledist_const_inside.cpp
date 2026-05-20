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
    int max_iterations = 8192;

    omp_set_num_threads(threads);

    Swarm::Point<T, dim> a(-5.f);
    Swarm::Point<T, dim> b(+5.f);

    Swarm::Point<T, dim> initial_best(.5f);

    Swarm::Optimizers::CHOPSOOptimizer<T, dim> swarm(fitness, initial_best, a, b, nN, nC, chaosMap, max_iterations);

    Swarm::Point<T, dim> center1(3.f);
    T radius1 = 5.f;
    auto dist_constraint1 = Swarm::makeMaxDistanceConstraint<T, dim>(center1, radius1);
    swarm.addConstraint(dist_constraint1);

    Swarm::Point<T, dim> center2(-3.f);
    T radius2 = 5.f;
    auto dist_constraint2 = Swarm::makeMaxDistanceConstraint<T, dim>(center2, radius2);
    swarm.addConstraint(dist_constraint2);

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
        for (size_t j = 0; j < 3; j++)
        {
            int particles = 10000 * (std::pow(2, 2 * j));

            std::shared_ptr<Swarm::ObjectiveFunction<T, dim>> fitness = std::make_shared<Swarm::ObjectiveFunctions::DropwaveFunction<T, dim>>();

            std::cout << "\n============================\n"
                      << "Benchmarking with " << (i + 1) << " threads and " << particles << " particles per type.\n";

            auto elapsed = measure(i + 1, particles, std::move(fitness), Swarm::ChaosFactory::Chebyshev<T, dim>());

            std::cout << "Time elapsed: " << elapsed.count() << " seconds.\n"
                      << "============================\n";

            benchmark[std::to_string(i * 3 + j + 1)] = {
                {"threads", i + 1},
                {"particles_per_type", particles},
                {"time_seconds", elapsed.count()}};
        }
    }

    std::ofstream benchmark_file("/work/u10768804/benchmark/benchmark_chopso_onlytimes_doubledist_const_inside.json");
    benchmark_file << benchmark.dump(4);
    benchmark_file.close();

    return 0;
}