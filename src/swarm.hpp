#ifndef _SWARM_HPP_
#define _SWARM_HPP_
#include "point.hpp"
#include "chaosmap.hpp"
#include "function.hpp"
#include "particle.hpp"
#include "../lib/funcs.hpp"

#include <cmath>
#include <vector>
#include <memory>
#include <limits>

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
        // const Function<T, dim> &fitness_function;
        std::unique_ptr<Function<T, dim>> fitness_function;

        Point<T, dim> global_best;
        float global_best_value;

        Point<T, dim> a;
        Point<T, dim> b;

        int current_iteration;
        int max_iterations;

    public:
        Swarm(std::unique_ptr<Function<T, dim>> &p, const Point<T, dim> &_a, const Point<T, dim> &_b) : particles(0), fitness_function(std::move(p)), global_best(0.f), global_best_value(std::numeric_limits<T>::infinity()), a(_a), b(_b), current_iteration(0), max_iterations(0) {}

        void addParticle(std::unique_ptr<Particle<T, dim>> p)
        {
            p->reinit(Point<T, dim>([this](size_t i)
                                    { return generate_random(a[i], b[i]); }));
            particles.emplace_back(std::move(p));
        }

        /*
            Description of the function findGlobalBest
            - parameters: none

            the function iterates through all particles to find the one with the best fitness value
        */
        void findGlobalBest()
        {
            for (const auto &particle : particles)
            {
                float fitness = fitness_function->evaluate(particle->getPosition());
                if (fitness < global_best_value)
                {
                    global_best_value = fitness;
                    global_best = particle->getPosition();
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
            for (auto &particle : particles)
            {
                std::cout << particle->getPosition() << std::endl;
                particle->updatePosition(this->global_best, this->a, this->b, current_iteration, max_iterations);
                particle->updatePersonalBest(*fitness_function);
            }

            findGlobalBest();

            ++current_iteration;
        }

        Point<T, dim> getGlobalBest() const
        {
            return global_best;
        }

        float getGlobalBestValue() const
        {
            return global_best_value;
        }
    };
}
#endif
