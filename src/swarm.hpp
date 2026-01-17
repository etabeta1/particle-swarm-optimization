#ifndef _SWARM_HPP_
#define _SWARM_HPP_

#include "point.hpp"
#include "chaosmap.hpp"
#include "function.hpp"
#include "particle.hpp"
#include "funcs.hpp"
#include "utils.hpp"

#include <cmath>
#include <vector>
#include <memory>
#include <limits>
#include <fstream>
#include <omp.h>

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
 * This macro defines an OpenMP reduction operation to find the best evaluated point among particles in the swarm.
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
        Point<T, dim> point;
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
    class Swarm
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
        EvaluatedPoint<T, dim> global_best{Point<T, dim>(T(0)), std::numeric_limits<T>::infinity()};

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
         * \brief File stream to log particle positions.
         */
        std::ofstream positions_file;

    public:
        /**
         * \brief Constructs a `Swarm` with the given parameters.
         * \param p A unique pointer to the fitness function.
         * \param _a The lower bounds of the search space.
         * \param _b The upper bounds of the search space.
         * \param _max_iterations The maximum number of iterations for the swarm.
         *
         * This constructor initializes the swarm with the provided fitness function, search space bounds, and maximum iterations.
         */
        Swarm(std::unique_ptr<ObjectiveFunction<T, dim>> &p, const Point<T, dim> &_a, const Point<T, dim> &_b, IterationType _max_iterations) : particles(0), fitness_function(std::move(p)), a(_a), b(_b), current_iteration(1), max_iterations(_max_iterations)
        {
            positions_file.open("points_xy.txt");
            if (!positions_file.is_open())
            {
                std::cerr << "Can't open the points file" << std::endl;
            }
        }

        /**
         * \brief Destructor for the `Swarm` class.
         *
         * Closes the positions file if it is open.
         */
        ~Swarm()
        {
            if (positions_file.is_open())
            {
                positions_file.close();
            }
        }

        /**
         * \brief Adds a particle to the swarm.
         * \param p A unique pointer to the particle to add.
         *
         * This function initializes the particle's position and personal best using random values within the search space bounds, and then adds it to the swarm.
         */
        void addParticle(std::unique_ptr<Particle<T, dim>> p)
        {
            p->reinit(Point<T, dim>([this](size_t i)
                                    { return generate_random(a[i], b[i]); }),
                      *fitness_function);
            particles.emplace_back(std::move(p));
        }

        /**
         * \brief Runs the swarm optimization algorithm.
         *
         * This function iteratively updates the positions of all particles in the swarm until the maximum number of iterations is reached.
         */
        void run()
        {
#pragma omp parallel default(private)
            {
                for (IterationType it = 0; it < max_iterations; ++it)
                {
                    updateEveryone();
                }
            }
        }

        /**
         * \brief Updates the positions of all particles in the swarm.
         */
        void updateEveryone()
        {
#pragma omp for schedule(static) reduction(findBestPoint : global_best)
            for (auto &particle : particles)
            {
                particle->updatePosition(global_best.point, a, b, current_iteration, max_iterations);
                bool update = particle->updatePersonalBest(*fitness_function);

                float fitness = particle->getPersonalBestValue();
                global_best = {particle->getPosition(), fitness};

                // if (positions_file.is_open())
                // {

                //     /*
                //     positions_file << pos[0] << " "                // X
                //                    << pos[1] << " "                // Y
                //                    << val << " "                   // Z
                //                    << current_iteration << " "     // Iterazione (k)
                //                    << particle->getType() << "\n"; // Type (0=Normal, 1=Chaotic)
                //                    */
                // }
            }

            ++current_iteration;
        }

        /**
         * \brief Returns the best evaluated point found by the swarm.
         * \return The best evaluated point.
         */
        EvaluatedPoint<T, dim> getGlobalBest() const
        {
            return global_best;
        }
    };
}

ENABLE_REDUCTION(float, 2);

#endif
