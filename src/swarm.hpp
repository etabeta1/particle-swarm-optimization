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

#define _TO_STRING(x) #x
#define TO_STRING(x) _TO_STRING(x)

#define ENABLE_REDUCTION(T, dim)                                                                           \
    _Pragma(TO_STRING(                                                                                     \
        omp declare reduction(                                                                             \
            findBestPoint : Swarm::EvaluatedPoint<T, dim> : betterPointReduction<T, dim>(omp_out, omp_in)) \
            initializer(omp_priv = {                                                                       \
                            Swarm::EvaluatedPoint<T, dim>({Swarm::Point<T, dim>(T(0)),                     \
                                                           std::numeric_limits<T>::infinity()})})))
namespace Swarm
{
    template <typename T = float, int dim = 2>
    struct EvaluatedPoint
    {
        Point<T, dim> point;
        float value;
    };
}

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
    // Class Declarations
    /*
        Definition of a swarm optimization algorithm.
        - `T` is the type used to store the coordinates (defaults to `float`)
        - `dim` is the number of dimensions for the vector (defaults to 2)

        - addParticle: adds a particle to the swarm
        - findGlobalBest: finds the best position among all particles
        - updateEveryone: updates the position of all particles in the swarm
    */
    template <typename T = float, int dim = 2>
    class Swarm
    {
    private:
        std::vector<std::unique_ptr<Particle<T, dim>>> particles; // stores both normal and chaotic particles
        // const ObjectiveFunction<T, dim> &fitness_function;
        std::unique_ptr<ObjectiveFunction<T, dim>> fitness_function;

        // Point<T, dim> global_best;
        // float global_best_value;

        EvaluatedPoint<T, dim> global_best{Point<T, dim>(T(0)), std::numeric_limits<T>::infinity()};

        Point<T, dim> a;
        Point<T, dim> b;

        IterationType current_iteration;
        IterationType max_iterations;

        std::ofstream positions_file;

    public:
        Swarm(std::unique_ptr<ObjectiveFunction<T, dim>> &p, const Point<T, dim> &_a, const Point<T, dim> &_b, IterationType _max_iterations) : particles(0), fitness_function(std::move(p)), a(_a), b(_b), current_iteration(1), max_iterations(_max_iterations)
        {
            positions_file.open("points_xy.txt");
            if (!positions_file.is_open())
            {
                std::cerr << "Can't open the points file" << std::endl;
            }
        }

        ~Swarm()
        {
            if (positions_file.is_open())
            {
                positions_file.close();
            }
        }

        void addParticle(std::unique_ptr<Particle<T, dim>> p)
        {
            p->reinit(Point<T, dim>([this](size_t i)
                                    { return generate_random(a[i], b[i]); }),
                      *fitness_function);
            particles.emplace_back(std::move(p));
        }

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

        /*
            Description of the function updateEveryone
            - parameters: none

            the function updates the position of all particles, updates the personal best for normal particles,
            finds the global best, and increments the current iteration
        */
        void updateEveryone()
        {
#pragma omp for schedule(static) reduction(findBestPoint : global_best)
            for (auto &particle : particles)
            {
                particle->updatePosition(global_best.point, a, b, current_iteration, max_iterations);
                bool update = particle->updatePersonalBest(*fitness_function);

                // Point<T, dim> pos = particle->getPosition();
                // float val = fitness_function->evaluate(pos);

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

            // findGlobalBest();

            ++current_iteration;
        }

        // Point<T, dim> getGlobalBest() const
        // {
        //     return global_best;
        // }

        // float getGlobalBestValue() const
        // {
        //     return global_best_value;
        // }

        EvaluatedPoint<T, dim> getGlobalBest() const
        {
            return global_best;
        }
    };
}

ENABLE_REDUCTION(float, 2);

#endif
