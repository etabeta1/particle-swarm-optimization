#ifndef _SWARM_HPP_
#define _SWARM_HPP_

#include <omp.h>

#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

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
namespace Swarm
{
    /**
     * \brief Structure representing a point and its evaluated value.
     * \tparam T The type used to store the coordinates (defaults to `float`).
     * \tparam dim The number of dimensions for the vector (defaults to 2).
     *
     * This structure holds a point in the search space along with its corresponding fitness value.
     */
    template <typename T = float, int dim = 2>
    struct EvaluatedPoint
    {
        /**
         * \brief The point in the search space.
         */
        Point<T, dim> point;
        /**
         * \brief The value of the fitness function at the point.
         */
        float value;
    };
}

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
         * \param initial_position The initial position for the normal particles.
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
                        const Point<T, dim> &initial_position,
                        const Point<T, dim> &_a,
                        const Point<T, dim> &_b,
                        size_t num_normal_particles,
                        size_t num_chaotic_particles,
                        const ChaosMap<T, dim> &chaos_map,
                        IterationType _max_iterations,
                        bool save_on_file)
            : particles(0),
              fitness_function(std::move(p)),
              global_best(initial_position, fitness_function->evaluate(initial_position)),
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

            if (!initial_position.isInsideBox(a, b))
            {
                throw std::invalid_argument("Initial position is not inside global search space");
            }

            for (size_t i = 0; i < num_normal_particles; i++)
            {
                std::unique_ptr<Particle<T, dim>> particle =
                    std::make_unique<NormalParticle<T, dim>>();
                particle->reinit(initial_position, *fitness_function);
                particles.emplace_back(std::move(particle));
            }

            for (size_t i = 0; i < num_chaotic_particles; i++)
            {
                std::unique_ptr<Particle<T, dim>> particle =
                    std::make_unique<ChaoticParticle<T, dim>>(chaos_map);
                particle->reinit(Point<T, dim>([this](size_t i)
                                               { return generate_random(a[i], b[i]); }),
                                 *fitness_function);
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
#pragma omp parallel default(private)
            {
                while (current_iteration < max_iterations)
                {
                    updateEveryone();
                }
            }
        }

        /**
         * \brief Updates the positions of all particles in the swarm and, if enabled, logs their positions to a file.
         */
        void updateEveryone()
        {
#pragma omp for schedule(static) reduction(findBestPoint : global_best)
            for (auto &particle : particles)
            {
                particle->updatePosition(global_best.point, a, b, current_iteration, max_iterations, constraints);
                // We cannot check here whether a particle is inside or outside the constraints
                // because each type of particles has a different behavior, so we make each particle
                // update itself.
                particle->updatePersonalBest(*fitness_function, constraints);

                float fitness = particle->getPersonalBestValue();
                global_best = {particle->getPosition(), fitness};

                if (positions_file.is_open())
                {
                    for (size_t i = 0; i < dim; i++)
                    {
                        positions_file << particle->getPosition()[i] << " ";
                    }
                    positions_file << fitness_function->evaluate(particle->getPosition()) << " " // Z
                                   << current_iteration << " "                                   // Iterazione (k)
                                   << particle->getType() << "\n";                               // Type (0=Normal, 1=Chaotic)
                }
            }

            ++current_iteration;
        }

        /**
         * \brief Returns the best evaluated point found by the swarm.
         * \return The best evaluated point.
         */
        EvaluatedPoint<T, dim> getGlobalBest() const { return global_best; }
    };
}

ENABLE_REDUCTION(float, 2);

#endif
