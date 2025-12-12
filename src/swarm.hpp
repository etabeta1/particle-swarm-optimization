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

        IterationType current_iteration;
        IterationType max_iterations;

        std::ofstream positions_file;

    public:
        Swarm(std::unique_ptr<Function<T, dim>> &p, const Point<T, dim> &_a, const Point<T, dim> &_b, IterationType _max_iterations) : particles(0), fitness_function(std::move(p)), global_best(T(0)), global_best_value(std::numeric_limits<T>::infinity()), a(_a), b(_b), current_iteration(1), max_iterations(_max_iterations)
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

        /*
            Description of the function findGlobalBest
            - parameters: none

            the function iterates through all particles to find the one with the best fitness value
        */
        void findGlobalBest()
        {
            for (const auto &particle : particles)
            {
                // float fitness = fitness_function->evaluate(particle->getPosition());
                float fitness = particle->getPersonalBestValue();
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
                particle->updatePosition(this->global_best, this->a, this->b, current_iteration, max_iterations);
                particle->updatePersonalBest(*fitness_function);

                if (positions_file.is_open())
                {
                    Point<T, dim> pos = particle->getPosition();
                    float val = fitness_function->evaluate(pos);

                    positions_file << pos[0] << " "                // X
                                   << pos[1] << " "                // Y
                                   << val << " "                   // Z
                                   << current_iteration << " "     // Iterazione (k)
                                   << particle->getType() << "\n"; // Type (0=Normal, 1=Chaotic)
                }
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
