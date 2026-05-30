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

    std::shared_ptr<Swarm::ObjectiveFunction<T, dim>> fitness2 = std::make_shared<Swarm::ObjectiveFunctions::DropwaveFunction<T, dim>>();
    Swarm::Optimizers::SAOptimizer<T, dim> swarm_alpine1(fitness2, initial_best, a, b, nN, nC, Swarm::ChaosFactory::Chebyshev<T, dim>(), max_iterations, 395.0f, 0.95f, "/work/u10768804/benchmark/SAopt_convergence_dropwave_Chebyshev.txt");
    swarm_alpine1.run();

    std::shared_ptr<Swarm::ObjectiveFunction<T, dim>> fitness3 = std::make_shared<Swarm::ObjectiveFunctions::AckleyFunction<T, dim>>();
    Swarm::Optimizers::SAOptimizer<T, dim> swarm_dropwave(fitness3, initial_best, a, b, nN, nC, Swarm::ChaosFactory::Chebyshev<T, dim>(), max_iterations, 395.0f, 0.95f, "/work/u10768804/benchmark/SAopt_convergence_ackley_Chebyshev.txt");
    swarm_dropwave.run();

    std::shared_ptr<Swarm::ObjectiveFunction<T, dim>> fitness4 = std::make_shared<Swarm::ObjectiveFunctions::Alpine1Function<T, dim>>();
    Swarm::Optimizers::SAOptimizer<T, dim> swarm_ellipsoid(fitness4, initial_best, a, b, nN, nC, Swarm::ChaosFactory::Chebyshev<T, dim>(), max_iterations, 395.0f, 0.95f, "/work/u10768804/benchmark/SAopt_convergence_ellipsoid_Chebyshev.txt");
    swarm_ellipsoid.run();

    std::shared_ptr<Swarm::ObjectiveFunction<T, dim>> fitness5 = std::make_shared<Swarm::ObjectiveFunctions::EllipsoidFunction<T, dim>>();
    Swarm::Optimizers::SAOptimizer<T, dim> swarm_quintic(fitness5, initial_best, a, b, nN, nC, Swarm::ChaosFactory::Chebyshev<T, dim>(), max_iterations, 395.0f, 0.95f, "/work/u10768804/benchmark/SAopt_convergence_ellipsoid_Chebyshev.txt");
    swarm_quintic.run();

    std::shared_ptr<Swarm::ObjectiveFunction<T, dim>> fitness6 = std::make_shared<Swarm::ObjectiveFunctions::QuinticFunction<T, dim>>();
    Swarm::Optimizers::SAOptimizer<T, dim> swarm_sphere(fitness6, initial_best, a, b, nN, nC, Swarm::ChaosFactory::Chebyshev<T, dim>(), max_iterations, 395.0f, 0.95f, "/work/u10768804/benchmark/SAopt_convergence_sphere_Chebyshev.txt");
    swarm_sphere.run();

    std::shared_ptr<Swarm::ObjectiveFunction<T, dim>> fitness7 = std::make_shared<Swarm::ObjectiveFunctions::SphereFunction<T, dim>>();
    Swarm::Optimizers::SAOptimizer<T, dim> swarm_sphere2(fitness7, initial_best, a, b, nN, nC, Swarm::ChaosFactory::Chebyshev<T, dim>(), max_iterations, 395.0f, 0.95f, "/work/u10768804/benchmark/SAopt_convergence_sphere_Chebyshev.txt");
    swarm_sphere2.run();
    return 0;
}