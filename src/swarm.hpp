#ifndef _SWARM_HPP_
#define _SWARM_HPP_

#include <omp.h>

#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>
#include <random>

#include "chaosmap.hpp"
#include "constraint.hpp"
#include "funcs.hpp"
#include "function.hpp"
#include "particle.hpp"
#include "point.hpp"
#include "utils.hpp"

/**
 * \brief Macro to convert a value to a string.
 */
#define _TO_STRING(x) #x

/**
 * \brief Macro to convert a value to a string.
 */
#define TO_STRING(x) _TO_STRING(x)

/**
 * \brief Macro to enable OpenMP reduction for finding the best evaluated point.
 * \param T The type used to store the coordinates (defaults to `float`).
 * \param dim The number of dimensions for the vector (defaults to 2).
 *
 * This macro defines an OpenMP reduction operation to find the best evaluated point among particles
 * in the swarm.
 */
#define ENABLE_REDUCTION(T, dim)                                                                           \
    _Pragma(TO_STRING(                                                                                     \
        omp declare reduction(                                                                             \
            findBestPoint : Swarm::EvaluatedPoint<T, dim> : betterPointReduction<T, dim>(omp_out, omp_in)) \
            initializer(omp_priv = {                                                                       \
                            Swarm::EvaluatedPoint<T, dim>({Swarm::Point<T, dim>(T(0)),                     \
                                                           std::numeric_limits<T>::infinity()})})))

/**
 * \brief Function to perform reduction for finding the better evaluated point.
 * \tparam T The type used to store the coordinates (defaults to `float`).
 * \tparam dim The number of dimensions for the vector (defaults to 2).
 * \param inout The current best evaluated point.
 * \param in The evaluated point to compare against.
 */
template <typename T = float, int dim = 2>
void betterPointReduction(Swarm::EvaluatedPoint<T, dim> &inout, Swarm::EvaluatedPoint<T, dim> &in)
{
    if (in.value < inout.value)
    {
        inout = in;
    }
}

namespace Swarm
{

    template <typename T = float, int dim = 2>
    class Particle;

    namespace Optimizers
    {

        /**
         * \brief Definition of the swarm optimization algorithm.
         * \tparam T The type used to store the coordinates (defaults to `float`).
         * \tparam dim The number of dimensions for the vector (defaults to 2).
         *
         * This class manages a swarm of particles and performs optimization using the CHOPSO algorithm.
         */
        template <typename T = float, int dim = 2>
        class CHOPSOOptimizer
        {
        private:
            /**
             * \brief The particles in the swarm.
             */
            std::vector<std::unique_ptr<Particle<T, dim>>> particles;

            /**
             * \brief The fitness function used to evaluate particle positions.
             */
            std::unique_ptr<ObjectiveFunction<T, dim>> fitness_function;

            /**
             * \brief The best evaluated point found by the swarm.
             */
            EvaluatedPoint<T, dim> global_best;

            /**
             * \brief The lower bounds of the search space.
             */
            Point<T, dim> a;

            /**
             * \brief The upper bounds of the search space.
             */
            Point<T, dim> b;

            /**
             * \brief The current iteration of the swarm.
             */
            IterationType current_iteration;

            /**
             * \brief The maximum number of iterations for the swarm.
             */
            IterationType max_iterations;

            /**
             * \brief The constraints applied to the swarm.
             */
            std::vector<Constraint<T, dim>> constraints;

            /**
             * \brief File stream to log particle positions.
             */
            std::ofstream positions_file;

        public:
            /**
             * \brief Constructs a `CHOPSOOptimizer` with the given parameters.
             * \param p A unique pointer to the fitness function.
             * \param initial_best The initial position for the normal particles.
             * \param _a The lower bounds of the search space.
             * \param _b The upper bounds of the search space.
             * \param num_normal_particles The number of normal particles in the swarm.
             * \param num_chaotic_particles The number of chaotic particles in the swarm.
             * \param chaos_map The chaos map used for chaotic particles.
             * \param _max_iterations The maximum number of iterations for the swarm.
             * \param save_on_file Flag indicating whether to save particle positions to a file.
             * \throws std::invalid_argument if the initial position is not inside the search space.
             *
             * This constructor initializes the swarm with the provided fitness function, search space
             * bounds, and maximum iterations.
             */
            CHOPSOOptimizer(std::unique_ptr<ObjectiveFunction<T, dim>> &p,
                            const Point<T, dim> &initial_best,
                            const Point<T, dim> &_a,
                            const Point<T, dim> &_b,
                            size_t num_normal_particles,
                            size_t num_chaotic_particles,
                            const ChaosMap<T, dim> &chaos_map,
                            IterationType _max_iterations,
                            bool save_on_file)
                : particles(0),
                  fitness_function(std::move(p)),
                  global_best(initial_best, fitness_function->evaluate(initial_best)),
                  a(_a),
                  b(_b),
                  current_iteration(1),
                  max_iterations(_max_iterations)
            {
                if (save_on_file)
                {
                    positions_file.open("points_xy.txt");
                    if (!positions_file.is_open())
                    {
                        std::cerr << "Can't open the points file" << std::endl;
                    }
                }

                if (!initial_best.isInsideBox(a, b))
                {
                    throw std::invalid_argument("Initial position is not inside global search space");
                }

                for (size_t i = 0; i < num_normal_particles; i++)
                {
                    std::unique_ptr<Particle<T, dim>> particle =
                        std::make_unique<NormalParticle<T, dim>>();
                    Point<T, dim> random_pos([this](size_t j)
                                             { return generate_random(a[j], b[j]); });
                    particle->reinit(random_pos, initial_best, *fitness_function);
                    particles.emplace_back(std::move(particle));
                }

                for (size_t i = 0; i < num_chaotic_particles; i++)
                {
                    std::unique_ptr<Particle<T, dim>> particle =
                        std::make_unique<ChaoticParticle<T, dim>>(chaos_map);
                    Point<T, dim> random_pos([this](size_t j)
                                             { return generate_random(a[j], b[j]); });
                    particle->reinit(random_pos, initial_best, *fitness_function);
                    particles.emplace_back(std::move(particle));
                }

                constraints.emplace_back(makeBoxConstraint(a, b));
            }

            /**
             * \brief Destructor for the `CHOPSOOptimizer` class.
             *
             * Closes the positions file if it is open.
             */
            ~CHOPSOOptimizer()
            {
                if (positions_file.is_open())
                {
                    positions_file.close();
                }
            }

            /**
             * \brief Adds a new constraint to the swarm.
             * \param constraint The constraint to be added.
             * \throws std::invalid_argument if the constraint invalidates the current global best position.
             */
            void addConstraint(const Constraint<T, dim> &constraint)
            {
                if (!(constraint(global_best.point)))
                {
                    throw std::invalid_argument("Constraint would invalidate current global best position.");
                }

                constraints.emplace_back(constraint);
            }

            /**
             * \brief Runs the swarm optimization algorithm.
             *
             * This function iteratively updates the positions of all particles in the swarm until the
             * maximum number of iterations is reached.
             */
            void run()
            {
                while (current_iteration < max_iterations)
                {
                    updateEveryone();
                }
            }

            /**
             * \brief Updates the positions of all particles in the swarm and, if enabled, logs their positions to a file.
             */
            void updateEveryone()
            {
                EvaluatedPoint<T, dim> old_global_best = global_best;
                EvaluatedPoint<T, dim> best_in_iteration = old_global_best;

#pragma omp parallel for default(shared) schedule(static) reduction(findBestPoint : best_in_iteration)
                for (size_t i = 0; i < particles.size(); ++i)
                {
                    auto &particle = particles[i];
                    particle->updatePosition(old_global_best.point, a, b, current_iteration, max_iterations, constraints);
                    // We cannot check here whether a particle is inside or outside the constraints
                    // because each type of particles has a different behavior, so we make each particle
                    // update itself.
                    particle->updatePersonalBest(*fitness_function, constraints);

                    if (positions_file.is_open())
                    {
#pragma omp critical
                        {
                            for (size_t j = 0; j < dim; j++)
                            {
                                positions_file << particle->getPosition()[j] << " ";
                            }
                            positions_file << fitness_function->evaluate(particle->getPosition()) << " " // Z
                                           << current_iteration << " "                                   // Iterazione (k)
                                           << particle->getType() << "\n";                               // Type (0=Normal, 1=Chaotic)
                        }
                    }

                    float fitness = particles[i]->getPersonalBestValue();
                    best_in_iteration = {particles[i]->getPersonalBestPosition(), fitness};
                }

                // Update global best only if we found something better (monotonic improvement)
                if (best_in_iteration.value < global_best.value)
                {
                    global_best = best_in_iteration;
                }

                ++current_iteration;
            }

            /**
             * \brief Returns the best evaluated point found by the swarm.
             * \return The best evaluated point.
             */
            EvaluatedPoint<T, dim> getGlobalBest() const { return global_best; }
        };

        /**
         * \brief Definition of the swarm optimization algorithm.
         * \tparam T The type used to store the coordinates (defaults to `float`).
         * \tparam dim The number of dimensions for the vector (defaults to 2).
         *
         * This class manages a swarm of particles and performs optimization using the CHOPSO algorithm.
         */
        template <typename T = float, int dim = 2>
        class GENETICOptimizer
        {
        private:
            /**
             * \brief The particles in the swarm.
             */
            std::vector<std::unique_ptr<Particle<T, dim>>> particles;

            /**
             * \brief The fitness function used to evaluate particle positions.
             */
            std::unique_ptr<ObjectiveFunction<T, dim>> fitness_function;

            /**
             * \brief The best evaluated point found by the swarm.
             */
            EvaluatedPoint<T, dim> global_best;

            /**
             * \brief The lower bounds of the search space.
             */
            Point<T, dim> a;

            /**
             * \brief The upper bounds of the search space.
             */
            Point<T, dim> b;

            /**
             * \brief The current iteration of the swarm.
             */
            IterationType current_iteration;

            /**
             * \brief The maximum number of iterations for the swarm.
             */
            IterationType max_iterations;

            /**
             * \brief The constraints applied to the swarm.
             */
            std::vector<Constraint<T, dim>> constraints;

            /**
             * \brief File stream to log particle positions.
             */
            std::ofstream positions_file;

        public:
            /**
             * \brief Constructs a `GENETICOptimizer` with the given parameters.
             * \param p A unique pointer to the fitness function.
             * \param initial_best The initial position for the normal particles.
             * \param _a The lower bounds of the search space.
             * \param _b The upper bounds of the search space.
             * \param num_normal_particles The number of normal particles in the swarm.
             * \param num_chaotic_particles The number of chaotic particles in the swarm.
             * \param chaos_map The chaos map used for chaotic particles.
             * \param _max_iterations The maximum number of iterations for the swarm.
             * \param rep_size Repository size of Global Best and Personal Best.
             * \param save_on_file Flag indicating whether to save particle positions to a file.
             * \throws std::invalid_argument if the initial position is not inside the search space.
             *
             * This constructor initializes the swarm with the provided fitness function, search space
             * bounds, and maximum iterations.
             */
            GENETICOptimizer(std::unique_ptr<ObjectiveFunction<T, dim>> &p,
                             const Point<T, dim> &initial_best,
                             const Point<T, dim> &_a,
                             const Point<T, dim> &_b,
                             size_t num_normal_particles,
                             size_t num_chaotic_particles,
                             const ChaosMap<T, dim> &chaos_map,
                             IterationType _max_iterations,
                             bool save_on_file)
                : particles(0),
                  fitness_function(std::move(p)),
                  global_best(initial_best, fitness_function->evaluate(initial_best)),
                  a(_a),
                  b(_b),
                  current_iteration(1),
                  max_iterations(_max_iterations)
            {

                if (save_on_file)
                {
                    positions_file.open("points_xy.txt");
                    if (!positions_file.is_open())
                    {
                        std::cerr << "Can't open the points file" << std::endl;
                    }
                }

                if (!initial_best.isInsideBox(a, b))
                {
                    throw std::invalid_argument("Initial position is not inside global search space");
                }

                for (size_t i = 0; i < num_normal_particles; i++)
                {
                    std::unique_ptr<Particle<T, dim>> particle =
                        std::make_unique<NormalParticle<T, dim>>();
                    Point<T, dim> random_pos([this](size_t j)
                                             { return generate_random(a[j], b[j]); });
                    particle->reinit(random_pos, initial_best, *fitness_function);
                    particles.emplace_back(std::move(particle));
                }

                for (size_t i = 0; i < num_chaotic_particles; i++)
                {
                    std::unique_ptr<Particle<T, dim>> particle =
                        std::make_unique<ChaoticParticle<T, dim>>(chaos_map);
                    Point<T, dim> random_pos([this](size_t j)
                                             { return generate_random(a[j], b[j]); });
                    particle->reinit(random_pos, initial_best, *fitness_function);
                    particles.emplace_back(std::move(particle));
                }

                constraints.emplace_back(makeBoxConstraint(a, b));
            }

            /**
             * \brief Destructor for the `GENETICOptimizer` class.
             *
             * Closes the positions file if it is open.
             */
            ~GENETICOptimizer()
            {
                if (positions_file.is_open())
                {
                    positions_file.close();
                }
            }

            /**
             * \brief Adds a new constraint to the swarm.
             * \param constraint The constraint to be added.
             * \throws std::invalid_argument if the constraint invalidates the current global best position.
             */
            void addConstraint(const Constraint<T, dim> &constraint)
            {
                if (!(constraint(global_best.point)))
                {
                    throw std::invalid_argument("Constraint would invalidate current global best position.");
                }

                constraints.emplace_back(constraint);
            }

            /**
             * \brief Runs the swarm optimization algorithm.
             *
             * This function iteratively updates the positions of all particles in the swarm until the
             * maximum number of iterations is reached.
             */
            void run()
            {
                while (current_iteration < max_iterations)
                {
                    updateEveryone();
                }
            }

            /**
             * \brief Updates the positions of all particles in the swarm and, if enabled, logs their positions to a file.
             */
            void updateEveryone()
            {
                EvaluatedPoint<T, dim> old_global_best = global_best;
                EvaluatedPoint<T, dim> best_in_iteration = old_global_best;

#pragma omp parallel for default(shared) schedule(static) reduction(findBestPoint : best_in_iteration)
                for (size_t i = 0; i < particles.size(); ++i)
                {
                    auto &particle = particles[i];
                    particle->updatePosition(old_global_best.point, a, b, current_iteration, max_iterations, constraints);
                    // We cannot check here whether a particle is inside or outside the constraints
                    // because each type of particles has a different behavior, so we make each particle
                    // update itself.
                    particle->updatePersonalBest(*fitness_function, constraints);

                    if (positions_file.is_open())
                    {
#pragma omp critical
                        {
                            for (size_t j = 0; j < dim; j++)
                            {
                                positions_file << particle->getPosition()[j] << " ";
                            }
                            positions_file << fitness_function->evaluate(particle->getPosition()) << " " // Z
                                           << current_iteration << " "                                   // Iterazione (k)
                                           << particle->getType() << "\n";                               // Type (0=Normal, 1=Chaotic)
                        }
                    }

                    float fitness = particles[i]->getPersonalBestValue();
                    best_in_iteration = {particles[i]->getPersonalBestPosition(), fitness};
                }

                // Update global best only if we found something better (monotonic improvement)
                if (best_in_iteration.value < global_best.value)
                {
                    global_best = best_in_iteration;
                }

                mutateGlobalBest();

                ++current_iteration;
            }

            void mutateGlobalBest()
            {
                // Initialize with dummy values
                EvaluatedPoint<T, dim> uniform_mutation = {Point<T, dim>(T(0)), std::numeric_limits<float>::infinity()};
                EvaluatedPoint<T, dim> gaussian_mutation = {Point<T, dim>(T(0)), std::numeric_limits<float>::infinity()};
                EvaluatedPoint<T, dim> non_uniform_mutation = {Point<T, dim>(T(0)), std::numeric_limits<float>::infinity()};
                EvaluatedPoint<T, dim> dimension_wise_mutation = {Point<T, dim>(T(0)), std::numeric_limits<float>::infinity()};
                EvaluatedPoint<T, dim> opposition_based_mutation = {Point<T, dim>(T(0)), std::numeric_limits<float>::infinity()};

                // Parallelize the 5 independent mutations
                #pragma omp parallel sections
                {
                    #pragma omp section
                    uniform_mutation = uniformMutation();
                    
                    #pragma omp section
                    gaussian_mutation = GaussianMutation();
                    
                    #pragma omp section
                    non_uniform_mutation = nonUniformMutation();
                    
                    #pragma omp section
                    dimension_wise_mutation = dimensionWiseMutation();
                    
                    #pragma omp section
                    opposition_based_mutation = oppositionBasedMutation();
                }

                // Funzione lambda per controllare se un punto soddisfa tutti i constraints
                auto satisfiesConstraints = [this](const Point<T, dim>& point) -> bool {
                    for (const auto& constraint : constraints) {
                        if (!constraint(point)) {
                            return false;  // Violato un constraint
                        }
                    }
                    return true;  // Soddisfa tutti i constraints
                };

                // Verifica constraints e raccogli candidati validi
                std::vector<EvaluatedPoint<T, dim>> valid_mutations;

                if (satisfiesConstraints(uniform_mutation.point))
                    valid_mutations.push_back(uniform_mutation);
                
                if (satisfiesConstraints(gaussian_mutation.point))
                    valid_mutations.push_back(gaussian_mutation);
                
                if (satisfiesConstraints(non_uniform_mutation.point))
                    valid_mutations.push_back(non_uniform_mutation);
                
                if (satisfiesConstraints(dimension_wise_mutation.point))
                    valid_mutations.push_back(dimension_wise_mutation);
                
                if (satisfiesConstraints(opposition_based_mutation.point))
                    valid_mutations.push_back(opposition_based_mutation);

                // Se almeno una mutazione valida, scegli la migliore
                if (!valid_mutations.empty()) {
                    auto best_mutation = std::min_element(valid_mutations.begin(), valid_mutations.end(),
                        [](const EvaluatedPoint<T, dim>& a, const EvaluatedPoint<T, dim>& b) {
                            return a.value < b.value;
                        });
                    
                    // Aggiorna global_best se la mutazione è migliore
                    if (best_mutation->value < global_best.value) {
                        global_best = *best_mutation;
                    }
                }
            }
            /**
             * \brief Performs random 1d coordinate mutation
             * \return The uniform mutation point.
             */
            EvaluatedPoint<T, dim> uniformMutation()
            {
                Point<T, dim> gb_point = global_best.point;
                std::mt19937 rng(std::random_device{}());

                std::uniform_int_distribution<int> dist(0, dim - 1);
                int x = dist(rng);
                std::uniform_real_distribution<T> dist_t(a[x], b[x]);
                T value = dist_t(rng);
                gb_point[x] = value;

                return {gb_point, fitness_function->evaluate(gb_point)};
            }
            EvaluatedPoint<T, dim> GaussianMutation()
            {
                Point<T, dim> gb_point = global_best.point;
                std::mt19937 rng(std::random_device{}());

                std::uniform_int_distribution<int> dist(0, dim - 1);
                int x = dist(rng);

                T sigma = static_cast<T>(0.05) * (b[x] - a[x]);

                std::normal_distribution<T> gauss(gb_point[x], sigma);
                T value = gauss(rng);

                // Clamp the value to stay within bounds
                value = std::max(a[x], std::min(b[x], value));
                gb_point[x] = value;

                return {gb_point, fitness_function->evaluate(gb_point)};
            }
            EvaluatedPoint<T, dim> nonUniformMutation()
            {
                Point<T, dim> gb_point = global_best.point;
                std::mt19937 rng(std::random_device{}());

                // 1. Random dimension
                std::uniform_int_distribution<int> dist(0, dim - 1);
                int x = dist(rng);

                // 2. Random direction (+ or -)
                std::uniform_real_distribution<T> uni01(0.0, 1.0);
                bool increase = uni01(rng) < 0.5;

                // 3. Time-adaptive factor (decreases with iterations)
                T r = uni01(rng);
                T b_param = 5.0; // shape parameter

                T delta = (b[x] - a[x]) * (1.0 - std::pow(r, std::pow(1.0 - static_cast<T>(current_iteration) / static_cast<T>(max_iterations), b_param)));

                // 4. Apply mutation with clamping
                if (increase)
                    gb_point[x] = std::max(a[x], std::min(b[x], gb_point[x] + delta));
                else
                    gb_point[x] = std::max(a[x], std::min(b[x], gb_point[x] - delta));

                return {gb_point, fitness_function->evaluate(gb_point)};
            }
            EvaluatedPoint<T, dim> dimensionWiseMutation()
            {
                Point<T, dim> gb_point = global_best.point;
                std::mt19937 rng(std::random_device{}());

                // Probability of mutating each dimension
                T pm = 1.0 / static_cast<T>(dim);

                std::uniform_real_distribution<T> uni01(0.0, 1.0);
                std::uniform_real_distribution<T> val_dist(0.0, 1.0);

                // For each dimension, decide whether to mutate it
                for (int i = 0; i < dim; ++i)
                {
                    if (uni01(rng) < pm) // Mutation probability
                    {
                        // Gaussian perturbation centered at current value
                        T sigma = 0.05 * (b[i] - a[i]);
                        std::normal_distribution<T> gauss(gb_point[i], sigma);
                        T value = gauss(rng);

                        // Clamp to bounds
                        gb_point[i] = std::max(a[i], std::min(b[i], value));
                    }
                }

                return {gb_point, fitness_function->evaluate(gb_point)};
            }
            EvaluatedPoint<T, dim> oppositionBasedMutation()
            {
                Point<T, dim> gb_point = global_best.point;

                // Opposition-based mutation: reflect each coordinate around the center
                // opposite[d] = lower[d] + upper[d] - current[d]
                for (int i = 0; i < dim; ++i)
                {
                    gb_point[i] = a[i] + b[i] - gb_point[i];
                }

                return {gb_point, fitness_function->evaluate(gb_point)};
            }
            /**
             * \brief Returns the best evaluated point found by the swarm.
             * \return The best evaluated point.
             */
            EvaluatedPoint<T, dim> getGlobalBest() const { return global_best; }
        };
    };
}

ENABLE_REDUCTION(float, 2);

#endif
