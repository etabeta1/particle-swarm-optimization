#include <chaosmap.hpp>
#include <point.hpp>
#include <swarm.hpp>
#include <function.hpp>
#include <utils.hpp>

#include <nlohmann/json.hpp>

#include <chrono>
#include <iostream>

int main()
{
    using T = float;
    constexpr int dim = 2;

    int nN = 32, nC = 32;
    int max_iterations = 8192;

    Swarm::Point<T, dim> a(-5.f);
    Swarm::Point<T, dim> b(+5.f);

    Swarm::Point<T, dim> initial_best(5.f);

    std::unique_ptr<Swarm::ObjectiveFunction<T, dim>> fitness2 = std::make_unique<Swarm::ObjectiveFunctions::AckleyFunction<T, dim>>();
    Swarm::Optimizers::CHOPSOOptimizer<T, dim> swarm_alpine1(fitness2, initial_best, a, b, nN, nC, Swarm::ChaosFactory::Iterative<T, dim>(), max_iterations, "/work/u10822715/benchmark/chopso_convergence_ackley_iterative.txt");
    swarm_alpine1.run();

    std::unique_ptr<Swarm::ObjectiveFunction<T, dim>> fitness3 = std::make_unique<Swarm::ObjectiveFunctions::AckleyFunction<T, dim>>();
    Swarm::Optimizers::CHOPSOOptimizer<T, dim> swarm_dropwave(fitness3, initial_best, a, b, nN, nC, Swarm::ChaosFactory::LogisticMap<T, dim>(), max_iterations, "/work/u10822715/benchmark/chopso_convergence_ackley_logistic.txt");
    swarm_dropwave.run();

    std::unique_ptr<Swarm::ObjectiveFunction<T, dim>> fitness4 = std::make_unique<Swarm::ObjectiveFunctions::AckleyFunction<T, dim>>();
    Swarm::Optimizers::CHOPSOOptimizer<T, dim> swarm_ellipsoid(fitness4, initial_best, a, b, nN, nC, Swarm::ChaosFactory::Sine<T, dim>(), max_iterations, "/work/u10822715/benchmark/chopso_convergence_ackley_sine.txt");
    swarm_ellipsoid.run();

    std::unique_ptr<Swarm::ObjectiveFunction<T, dim>> fitness5 = std::make_unique<Swarm::ObjectiveFunctions::AckleyFunction<T, dim>>();
    Swarm::Optimizers::CHOPSOOptimizer<T, dim> swarm_quintic(fitness5, initial_best, a, b, nN, nC, Swarm::ChaosFactory::Singer<T, dim>(), max_iterations, "/work/u10822715/benchmark/chopso_convergence_ackley_singer.txt");
    swarm_quintic.run();

    std::unique_ptr<Swarm::ObjectiveFunction<T, dim>> fitness6 = std::make_unique<Swarm::ObjectiveFunctions::AckleyFunction<T, dim>>();
    Swarm::Optimizers::CHOPSOOptimizer<T, dim> swarm_sphere(fitness6, initial_best, a, b, nN, nC, Swarm::ChaosFactory::Sinusoidal<T, dim>(), max_iterations, "/work/u10822715/benchmark/chopso_convergence_ackley_sinusoidal.txt");
    swarm_sphere.run();

    return 0;
}