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
#include "fuzzylogic.hpp"
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
        protected:
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
                            IterationType _max_iterations)
                : particles(0),
                  fitness_function(std::move(p)),
                  global_best(initial_best, fitness_function->evaluate(initial_best)),
                  a(_a),
                  b(_b),
                  current_iteration(1),
                  max_iterations(_max_iterations)
            {
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
             * \throws std::runtime_error if the output file cannot be opened.
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
                            std::string output_file)
                : CHOPSOOptimizer<T, dim>(p,
                                          initial_best,
                                          _a,
                                          _b,
                                          num_normal_particles,
                                          num_chaotic_particles,
                                          chaos_map,
                                          _max_iterations)
            {
                positions_file.open(output_file);

                if (!positions_file.is_open())
                {
                    throw std::runtime_error("Could not open output file for writing.");
                }
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
            virtual void run()
            {
                while (current_iteration < max_iterations)
                {
                    updateEveryone();
                }
            }

            /**
             * \brief Updates the positions of all particles in the swarm and, if enabled, logs their positions to a file.
             */
            virtual void updateEveryone()
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
        class GENETICOptimizer : protected CHOPSOOptimizer<T, dim>
        {

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
                             IterationType _max_iterations)
                : CHOPSOOptimizer<T, dim>(p,
                                          initial_best,
                                          _a,
                                          _b,
                                          num_normal_particles,
                                          num_chaotic_particles,
                                          chaos_map,
                                          _max_iterations) {};

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
             * \param output_filename The name of the output file to save particle positions.
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
                             std::string output_filename)
                : CHOPSOOptimizer<T, dim>(p,
                                          initial_best,
                                          _a,
                                          _b,
                                          num_normal_particles,
                                          num_chaotic_particles,
                                          chaos_map,
                                          _max_iterations,
                                          output_filename) {};

            /**
             * \brief Updates the positions of all particles in the swarm and, if enabled, logs their positions to a file.
             */
            using super = CHOPSOOptimizer<T, dim>;
            void updateEveryone() override
            {
                super::updateEveryone();
                mutateGlobalBest();
            }
            void run() override
            {
                super::run();
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
                auto satisfiesConstraints = [this](const Point<T, dim> &point) -> bool
                {
                    for (const auto &constraint : super::constraints)
                    {
                        if (!constraint(point))
                        {
                            return false; // Violato un constraint
                        }
                    }
                    return true; // Soddisfa tutti i constraints
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
                if (!valid_mutations.empty())
                {
                    auto best_mutation = std::min_element(valid_mutations.begin(), valid_mutations.end(),
                                                          [](const EvaluatedPoint<T, dim> &a, const EvaluatedPoint<T, dim> &b)
                                                          {
                                                              return a.value < b.value;
                                                          });

                    // Aggiorna global_best se la mutazione è migliore
                    if (best_mutation->value < super::global_best.value)
                    {
                        super::global_best = *best_mutation;
                    }
                }
            }
            /**
             * \brief Performs random 1d coordinate mutation
             * \return The uniform mutation point.
             */
            EvaluatedPoint<T, dim> uniformMutation()
            {
                Point<T, dim> gb_point = super::global_best.point;
                std::mt19937 rng(std::random_device{}());

                std::uniform_int_distribution<int> dist(0, dim - 1);
                int x = dist(rng);
                std::uniform_real_distribution<T> dist_t(super::a[x], super::b[x]);
                T value = dist_t(rng);
                gb_point[x] = value;

                return {gb_point, super::fitness_function->evaluate(gb_point)};
            }
            EvaluatedPoint<T, dim> GaussianMutation()
            {
                Point<T, dim> gb_point = super::global_best.point;
                std::mt19937 rng(std::random_device{}());

                std::uniform_int_distribution<int> dist(0, dim - 1);
                int x = dist(rng);

                T sigma = static_cast<T>(0.05) * (super::b[x] - super::a[x]);

                std::normal_distribution<T> gauss(gb_point[x], sigma);
                T value = gauss(rng);

                // Clamp the value to stay within bounds
                value = std::max(super::a[x], std::min(super::b[x], value));
                gb_point[x] = value;

                return {gb_point, super::fitness_function->evaluate(gb_point)};
            }
            EvaluatedPoint<T, dim> nonUniformMutation()
            {
                Point<T, dim> gb_point = super::global_best.point;
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

                T delta = (super::b[x] - super::a[x]) * (1.0 - std::pow(r, std::pow(1.0 - static_cast<T>(super::current_iteration) / static_cast<T>(super::max_iterations), b_param)));

                // 4. Apply mutation with clamping
                if (increase)
                    gb_point[x] = std::max(super::a[x], std::min(super::b[x], gb_point[x] + delta));
                else
                    gb_point[x] = std::max(super::a[x], std::min(super::b[x], gb_point[x] - delta));

                return {gb_point, super::fitness_function->evaluate(gb_point)};
            }
            EvaluatedPoint<T, dim> dimensionWiseMutation()
            {
                Point<T, dim> gb_point = super::global_best.point;
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
                        T sigma = 0.05 * (super::b[i] - super::a[i]);
                        std::normal_distribution<T> gauss(gb_point[i], sigma);
                        T value = gauss(rng);

                        // Clamp to bounds
                        gb_point[i] = std::max(super::a[i], std::min(super::b[i], value));
                    }
                }

                return {gb_point, super::fitness_function->evaluate(gb_point)};
            }
            EvaluatedPoint<T, dim> oppositionBasedMutation()
            {
                Point<T, dim> gb_point = super::global_best.point;

                // Opposition-based mutation: reflect each coordinate around the center
                // opposite[d] = lower[d] + upper[d] - current[d]
                for (int i = 0; i < dim; ++i)
                {
                    gb_point[i] = super::a[i] + super::b[i] - gb_point[i];
                }

                return {gb_point, super::fitness_function->evaluate(gb_point)};
            }
            /**
             * \brief Returns the best evaluated point found by the swarm.
             * \return The best evaluated point.
             */
            EvaluatedPoint<T, dim> getGlobalBest() const { return super::global_best; }
        };

        /**
         * \brief Definition of the Fuzzy Self Tuning PSO algorithm
         *
         * \tparam T The type used to store the coordinates (defaults to `float`).
         * \tparam dim The number of dimensions for the vector (defaults to 2).
         *
         * This class manages a swarm of particles and performs optimization using the FSTPSO algorithm.
         */
        template <typename T = float, int dim = 2>
        class FSTPSOOptimizer : protected CHOPSOOptimizer<T, dim>
        {
        protected:
            using super = CHOPSOOptimizer<T, dim>;

            /**
             * \brief Maximum distance in the search space
             */
            T delta_max;

            /**
             * \brief Worst fitness value
             */
            T worst_fitness;

            /**
             * \brief Flag in order to identify the first iteration
             *
             */
            bool first_iteration;

        public:
            /**
             * \brief Construct a new FSTPSOOptimizer object
             *
             * \param p pointer to the fitness function
             * \param initial_best initial position for the particles
             * \param _a  lower bound for the search space
             * \param _b  upper bound for the search space
             * \param num_particles number of fuzzy particles
             * \param _max_iterations the number of maximum iterations
             */
            FSTPSOOptimizer(std::unique_ptr<ObjectiveFunction<T, dim>> &p,
                            const Point<T, dim> &initial_best,
                            const Point<T, dim> &_a,
                            const Point<T, dim> &_b,
                            size_t num_particles,
                            IterationType _max_iterations)
                : CHOPSOOptimizer<T, dim>(p, initial_best, _a, _b, ChaosMap<T, dim>(), _max_iterations),
                  first_iteration(true)
            {
                // the standard implementation follows what defined by Hansen
                size_t N = (num_particles == 0) ? static_cast<size_t>(10 + 2 * std::sqrt(dim)) : num_particles;

                super::particles.clear();

                for (size_t i = 0; i < N; i++)
                {
                    auto particle = std::make_unique<FuzzyParticle<T, dim>>();
                    Point<T, dim> random_pos([this](size_t j)
                                             { return generate_random(super::a[j], super::b[j]); });
                    particle->reinit(random_pos, initial_best, *super::fitness_function);
                    super::particles.emplace_back(std::move(particle));
                }

                initializeFSTPSO();
            }

            /**
             * \brief Construct a new FSTPSOOptimizer object + file output
             *
             * \param p pointer to the fitness function
             * \param initial_best initial position for the particles
             * \param _a  lower bound for the search space
             * \param _b  upper bound for the search space
             * \param num_particles number of fuzzy particles
             * \param _max_iterations the number of maximum iterations
             */
            FSTPSOOptimizer(std::unique_ptr<ObjectiveFunction<T, dim>> &p,
                            const Point<T, dim> &initial_best,
                            const Point<T, dim> &_a,
                            const Point<T, dim> &_b,
                            size_t num_particles,
                            IterationType _max_iterations,
                            std::string output_file)
                : CHOPSOOptimizer<T, dim>(p, initial_best, _a, _b, ChaosMap<T, dim>(), _max_iterations, output_file),
                  first_iteration(true)
            {
                // the standard implementation follows what defined by Hansen
                size_t N = (num_particles == 0) ? static_cast<size_t>(10 + 2 * std::sqrt(dim)) : num_particles;

                super::particles.clear();

                for (size_t i = 0; i < N; i++)
                {
                    auto particle = std::make_unique<FuzzyParticle<T, dim>>();
                    Point<T, dim> random_pos([this](size_t j)
                                             { return generate_random(super::a[j], super::b[j]); });
                    particle->reinit(random_pos, initial_best, *super::fitness_function);
                    super::particles.emplace_back(std::move(particle));
                }

                initializeFSTPSO();
            }

            /**
             * \brief Updates all fuzzy particles using FST-PSO fuzzy logic.
             */
            void updateEveryone() override
            {
                if (first_iteration)
                {
                    worst_fitness = -std::numeric_limits<T>::infinity();
                    for (auto &particle : super::particles)
                    {
                        T fit = super::fitness_function->evaluate(particle->getPosition());
                        if (fit > worst_fitness)
                            worst_fitness = fit;

                        // Initialize previous fitness for FuzzyParticles
                        auto *fuzzy_particle = static_cast<FuzzyParticle<T, dim> *>(particle.get());
                        fuzzy_particle->updatePreviousFitness(fit);
                    }
                    first_iteration = false;
                }

                // crucial for phi calculation
                EvaluatedPoint<T, dim> old_global_best = super::global_best;
                EvaluatedPoint<T, dim> best_in_iteration = old_global_best;

#pragma omp parallel for default(shared) schedule(static) reduction(findBestPoint : best_in_iteration)
                for (size_t i = 0; i < super::particles.size(); ++i)
                {
                    auto &particle = super::particles[i];

                    // Cast to FuzzyParticle to access FST-PSO methods
                    auto *fuzzy_particle = static_cast<FuzzyParticle<T, dim> *>(particle.get());

                    // Update fuzzy parameters before position update
                    updateParticleFuzzyParameters(fuzzy_particle, old_global_best.point);

                    // Update position (particle uses its individual fuzzy parameters)
                    particle->updatePosition(old_global_best.point,
                                             super::a, super::b,
                                             super::current_iteration,
                                             super::max_iterations,
                                             super::constraints);

                    particle->updatePersonalBest(*super::fitness_function, super::constraints);

                    // Update previous fitness for next iteration's phi calculation
                    T current_fit = super::fitness_function->evaluate(particle->getPosition());
                    fuzzy_particle->updatePreviousFitness(current_fit);

                    if (super::positions_file.is_open())
                    {
#pragma omp critical
                        {
                            for (size_t j = 0; j < dim; j++)
                            {
                                super::positions_file << particle->getPosition()[j] << " ";
                            }
                            super::positions_file << current_fit << " "
                                                  << super::current_iteration << " "
                                                  << particle->getType() << "\n";
                        }
                    }

                    float fitness = super::particles[i]->getPersonalBestValue();
                    best_in_iteration = {super::particles[i]->getPersonalBestPosition(), fitness};
                }

                // Update global best only if we found something better
                if (best_in_iteration.value < super::global_best.value)
                {
                    super::global_best = best_in_iteration;
                }

                ++super::current_iteration;
            }

            /**
             * \brief Returns the best evaluated point found by the swarm.
             * \return The best evaluated point.
             */
            EvaluatedPoint<T, dim> getGlobalBest() const { return super::global_best; }

        private:
            /**
             * \brief initialize FST-PSO specific parameters
             */
            void initializeFSTPSO()
            {
                // Calculate delta_max (diagonal of search space)
                T sum_sq = static_cast<T>(0);
                for (int i = 0; i < dim; ++i)
                {
                    T diff = super::b[i] - super::a[i];
                    sum_sq += diff * diff;
                }
                delta_max = std::sqrt(sum_sq);
            }

            /**
             * \brief Core FSTPSO: Update fuzzy parameters for a particle.
             * \param particle Pointer to the FuzzyParticle to update.
             * \param global_best_pos Current global best position.
             *
             * Applies the complete FST-PSO fuzzy inference chain:
             * - Calculate delta
             * - Calculate phi
             * - Calculate membership degrees
             * - Apply 15 fuzzy rules (3 per output variable)
             * - Defuzzify using Sugeno method to obtain crisp values
             * - Set particle's individual parameters (w, c_cog, c_soc, L, U)
             *
             * \see computeDelta
             * \see computePhi
             * \see computeInertia
             * \see computeSocial
             * \see computeCognitive
             * \see computeLowerClamp
             * \see computeUpperClamp
             */
            void updateParticleFuzzyParameters(FuzzyParticle<T, dim> *particle,
                                               const Point<T, dim> &global_best_pos)
            {
                T delta = computeDelta(particle->getPosition(), global_best_pos);

                T phi = computePhi(particle);

                T w = computeInertia(delta, phi);
                T c_soc = computeSocial(delta, phi);
                T c_cog = computeCognitive(delta, phi);
                T L = computeLowerClamp(delta, phi);
                T U = computeUpperClamp(delta, phi);

                particle->setFuzzyParameters(w, c_soc, c_cog, L, U);
            }

            /**
             * \brief Compute Euclidean distance between two points.
             * \param pos Current particle position.
             * \param g_best Global best position.
             * \return Euclidean distance
             *
             * This distance is normalized by delta_max when used in fuzzy inference.
             */
            T computeDelta(const Point<T, dim> &pos, const Point<T, dim> &g_best)
            {
                T sum = static_cast<T>(0);
                for (int i = 0; i < dim; ++i)
                {
                    T diff = pos[i] - g_best[i];
                    sum += diff * diff;
                }
                return std::sqrt(sum);
            }

            /**
             * \brief Compute normalized fitness incremental factor (Equation 7 from paper).
             * \param particle Pointer to the FuzzyParticle.
             * \return Normalized fitness incremental factor
             */
            T computePhi(FuzzyParticle<T, dim> *particle)
            {
                Point<T, dim> current_pos = particle->getPosition();
                auto prev_eval = particle->getPreviousEval();
                Point<T, dim> prev_pos = prev_eval.point;
                T prev_fit = prev_eval.value;
                T current_fit = super::fitness_function->evaluate(current_pos);

                // Distance moved (normalized by delta_max)
                T dist_moved = computeDelta(current_pos, prev_pos);
                T normalized_dist = dist_moved / delta_max;

                // Fitness variation (clamped to worst_fitness, normalized)
                T fit_current_clamped = std::min(current_fit, worst_fitness);
                T fit_prev_clamped = std::min(prev_fit, worst_fitness);
                T fit_variation = (fit_current_clamped - fit_prev_clamped) / std::abs(worst_fitness);

                // Phi = normalized_distance × fitness_variation
                T phi = normalized_dist * fit_variation;

                // Clamp to [-1, 1]
                return std::max(static_cast<T>(-1.0), std::min(static_cast<T>(1.0), phi));
            }

            /**
             * \brief Fuzzy rules for Inertia (Rules 1, 2, 3 from Table 1).
             * \param delta Distance from global best (normalized).
             * \param phi Normalized fitness incremental factor.
             * \return Defuzzified inertia weight w ∈ [0.3, 1.0]
             *
             * Increase inertia when particle is improving or far from global best
             * to encourage exploration. Decrease when degrading or close to global best.
             */
            T computeInertia(T delta, T phi)
            {
                using namespace FuzzyLogic;

                // Calculate membership degrees for delta
                T delta_same = DeltaMembership<T>::membershipSame(delta, delta_max);
                T delta_near = DeltaMembership<T>::membershipNear(delta, delta_max);
                T delta_far = DeltaMembership<T>::membershipFar(delta, delta_max);

                // Calculate membership degrees for phi
                T phi_better = PhiMembership<T>::membershipBetter(phi);
                T phi_same = PhiMembership<T>::membershipSame(phi);
                T phi_worse = PhiMembership<T>::membershipWorse(phi);

                // Rule 1: if (phi is Worse OR delta is Same) then (Inertia is Low)
                T low_degree = std::max(phi_worse, delta_same);

                // Rule 2: if (phi is Same OR delta is Near) then (Inertia is Medium)
                T med_degree = std::max(phi_same, delta_near);

                // Rule 3: if (phi is Better OR delta is Far) then (Inertia is High)
                T high_degree = std::max(phi_better, delta_far);

                // Defuzzify using Sugeno method
                return sugenoDefuzzify<T>(low_degree, static_cast<T>(0.3),
                                          med_degree, static_cast<T>(0.5),
                                          high_degree, static_cast<T>(1.0));
            }

            /**
             * \brief Fuzzy rules for Social factor (Rules 4, 5, 6 from Table 1).
             * \param delta Distance from global best (normalized).
             * \param phi Normalized fitness incremental factor.
             * \return Defuzzified social factor c_soc ∈ [1.0, 3.0]
             *
             * Increase social attraction when particle is degrading or far
             * to bring it closer to swarm. Decrease when improving or close.
             */
            T computeSocial(T delta, T phi)
            {
                using namespace FuzzyLogic;

                T delta_same = DeltaMembership<T>::membershipSame(delta, delta_max);
                T delta_near = DeltaMembership<T>::membershipNear(delta, delta_max);
                T delta_far = DeltaMembership<T>::membershipFar(delta, delta_max);

                T phi_better = PhiMembership<T>::membershipBetter(phi);
                T phi_same = PhiMembership<T>::membershipSame(phi);
                T phi_worse = PhiMembership<T>::membershipWorse(phi);

                // Rule 4: if (phi is Better OR delta is Near) then (Social is Low)
                T low_degree = std::max(phi_better, delta_near);

                // Rule 5: if (phi is Same OR delta is Same) then (Social is Medium)
                T med_degree = std::max(phi_same, delta_same);

                // Rule 6: if (phi is Worse OR delta is Far) then (Social is High)
                T high_degree = std::max(phi_worse, delta_far);

                return sugenoDefuzzify<T>(low_degree, static_cast<T>(1.0),
                                          med_degree, static_cast<T>(2.0),
                                          high_degree, static_cast<T>(3.0));
            }

            /**
             * \brief Fuzzy rules for Cognitive factor (Rules 7, 8, 9 from Table 1).
             * \param delta Distance from global best (normalized).
             * \param phi Normalized fitness incremental factor.
             * \return Defuzzified cognitive factor c_cog ∈ [0.1, 3.0]
             *
             * Increase cognitive attraction when improving to exploit
             * local search around personal best. Decrease when far from global best.
             */
            T computeCognitive(T delta, T phi)
            {
                using namespace FuzzyLogic;

                T delta_same = DeltaMembership<T>::membershipSame(delta, delta_max);
                T delta_near = DeltaMembership<T>::membershipNear(delta, delta_max);
                T delta_far = DeltaMembership<T>::membershipFar(delta, delta_max);

                T phi_better = PhiMembership<T>::membershipBetter(phi);
                T phi_same = PhiMembership<T>::membershipSame(phi);
                T phi_worse = PhiMembership<T>::membershipWorse(phi);

                // Rule 7: if (delta is Far) then (Cognitive is Low)
                T low_degree = delta_far;

                // Rule 8: if (phi is Worse OR phi is Same OR delta is Same OR delta is Near)
                //         then (Cognitive is Medium)
                T med_degree = std::max({phi_worse, phi_same, delta_same, delta_near});

                // Rule 9: if (phi is Better) then (Cognitive is High)
                T high_degree = phi_better;

                return sugenoDefuzzify<T>(low_degree, static_cast<T>(0.1),
                                          med_degree, static_cast<T>(1.5),
                                          high_degree, static_cast<T>(3.0));
            }

            /**
             * \brief Fuzzy rules for Lower velocity clamp (Rules 10, 11, 12 from Table 1).
             * \param delta Distance from global best (normalized).
             * \param phi Normalized fitness incremental factor.
             * \return Defuzzified lower clamp factor L ∈ [0.0, 0.01]
             *
             * Increase minimum velocity when degrading to force exploration.
             * Decrease when improving or far to allow local exploitation.
             */
            T computeLowerClamp(T delta, T phi)
            {
                using namespace FuzzyLogic;

                T delta_same = DeltaMembership<T>::membershipSame(delta, delta_max);
                T delta_near = DeltaMembership<T>::membershipNear(delta, delta_max);
                T delta_far = DeltaMembership<T>::membershipFar(delta, delta_max);

                T phi_better = PhiMembership<T>::membershipBetter(phi);
                T phi_same = PhiMembership<T>::membershipSame(phi);
                T phi_worse = PhiMembership<T>::membershipWorse(phi);

                // Rule 10: if (phi is Same OR phi is Better OR delta is Far) then (L is Low)
                T low_degree = std::max({phi_same, phi_better, delta_far});

                // Rule 11: if (delta is Same OR delta is Near) then (L is Medium)
                T med_degree = std::max(delta_same, delta_near);

                // Rule 12: if (phi is Worse) then (L is High)
                T high_degree = phi_worse;

                return sugenoDefuzzify<T>(low_degree, static_cast<T>(0.0),
                                          med_degree, static_cast<T>(0.001),
                                          high_degree, static_cast<T>(0.01));
            }

            /**
             * \brief Fuzzy rules for Upper velocity clamp (Rules 13, 14, 15 from Table 1).
             * \param delta Distance from global best (normalized).
             * \param phi Normalized fitness incremental factor.
             * \return Defuzzified upper clamp factor U ∈ [0.1, 0.2]
             *
             * Increase maximum velocity when degrading or far to allow
             * larger jumps for exploration. Decrease when close to global best.
             */
            T computeUpperClamp(T delta, T phi)
            {
                using namespace FuzzyLogic;

                T delta_same = DeltaMembership<T>::membershipSame(delta, delta_max);
                T delta_near = DeltaMembership<T>::membershipNear(delta, delta_max);
                T delta_far = DeltaMembership<T>::membershipFar(delta, delta_max);

                T phi_better = PhiMembership<T>::membershipBetter(phi);
                T phi_same = PhiMembership<T>::membershipSame(phi);
                T phi_worse = PhiMembership<T>::membershipWorse(phi);

                // Rule 13: if (delta is Same) then (U is Low)
                T low_degree = delta_same;

                // Rule 14: if (phi is Same OR phi is Better OR delta is Near) then (U is Medium)
                T med_degree = std::max({phi_same, phi_better, delta_near});

                // Rule 15: if (phi is Worse OR delta is Far) then (U is High)
                T high_degree = std::max(phi_worse, delta_far);

                return sugenoDefuzzify<T>(low_degree, static_cast<T>(0.1),
                                          med_degree, static_cast<T>(0.15),
                                          high_degree, static_cast<T>(0.2));
            }
        };

    };
}

ENABLE_REDUCTION(float, 2);

#endif
