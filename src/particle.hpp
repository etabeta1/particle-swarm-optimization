#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include "point.hpp"
#include "../lib/funcs.hpp"

namespace Swarm
{

    /*
        Definition of a generic particle in the swarm.
        - `T` is the type used to store the coordinates (defaults to `float`)
        - `dim` is the number of dimensions for the vector (defaults to 2)

        - updatePosition: updates the position of the particle
    */
    template <typename T, int dim>
    class Particle
    {
    protected:
        Point<T, dim> position;
        Point<T, dim> personal_best;
        float personal_best_value;

    public:
        Particle() : position(T(0)), personal_best(T(0)), personal_best_value(T(0)) {}

        void reinit(const Point<T, dim> &initial_position, const Function<T, dim> &func)
        {
            position = initial_position;
            personal_best = initial_position;
            personal_best_value = func.evaluate(initial_position);
        }

        virtual void updatePosition(const Point<T, dim> &global_best, const Point<T, dim> &a, const Point<T, dim> &b, int current_iteration, int max_iterations) = 0;
        virtual ~Particle() {}

        /*
            Description of the function updatePersonalBest
            - parameters:
                - func: the fitness function to evaluate the particle's position
                - current_iteration: the current iteration of the swarm

            the function updates the personal best position and value of the normal particle if the current position is better
        */
        void updatePersonalBest(const Function<T, dim> &func)
        {
            float current_value = func.evaluate(this->position);
            if (current_value < this->personal_best_value)
            {
                this->personal_best_value = current_value;
                this->personal_best = this->position;
            }
        }

        Point<T, dim> &getPosition()
        {
            return position;
        }

        Point<T, dim> &getPersonalBest()
        {
            return personal_best;
        }

        float getPersonalBestValue()
        {
            return personal_best_value;
        }
    };

    /*
        Definition of the first kind of particle
        The normal particle follows the PSO logic

        - `T` is the type used to store the coordinates (defaults to `float`)
        - `dim` is the number of dimensions for the vector (defaults to 2)

        the constructor initializes the position, speed, and personal best
        constants c1 and c2 control the influence of personal and global best positions and are fixed based on the CHOPSO implementation

        - updatePosition: updates the position of the particle
        - updatePersonalBest: updates the personal best position and value of the particle
        - updateSpeed: updates the speed of the particle

        the last two functions are specific to the normal particle

    */
    template <typename T = float, int dim = 2>
    class NormalParticle : public Particle<T, dim>
    {
    private:
        Point<T, dim> speed;

        float c1 = 2.5f;
        float c2 = 2.5f;

        /*
            Description of the function updateSpeed
            - parameters:
                - global_best: the best position found by the swarm
                - current_iteration: the current iteration of the swarm
                - max_iterations: the maximum number of iterations for the swarm

            the function updates the speed of the normal particle based on the CHOPSO velocity update formula
        */
        void updateSpeed(const Point<T, dim> &global_best, int current_iteration, int max_iterations)
        {
            // CHOPSO velocity update formula
            // v(t+1) = w*v(t) + k1*(pBest - x(t)) + k2*(gBest - x(t))
            // k1 = rand[0,1] * 2.5
            // k2 = rand[0,1] * 2.5
            // w = 0.9 - (0.5 * iteration) / maxIterations) is the dynamic inertia weight

            float w = 0.9f - (0.5f * static_cast<float>(current_iteration) / static_cast<float>(max_iterations));

            float k1 = generate_random(0.0f, 1.0f) * c1;
            float k2 = generate_random(0.0f, 1.0f) * c2;

            Point<T, dim> inertia = speed * w;
            Point<T, dim> cognitive = (this->personal_best - this->position) * k1;
            Point<T, dim> social = (global_best - this->position) * k2;

            this->speed = inertia + cognitive + social;
        }

    public:
        NormalParticle() : Particle<T, dim>(), speed(0.f) {}

        /*
            Description of the function updatePosition
            - parameters:
                - global_best: the best position found by the swarm
                - a: the minimum boundary for the position
                - b: the maximum boundary for the position
                - current_iteration: the current iteration of the swarm
                - max_iterations: the maximum number of iterations for the swarm

            the function updates the position of the normal particle based on its speed and clamps it within the boundaries
        */
        void updatePosition(const Point<T, dim> &global_best, const Point<T, dim> &a, const Point<T, dim> &b, int current_iteration, int max_iterations) override
        {
            updateSpeed(global_best, current_iteration, max_iterations);

            this->position = (this->position + this->speed).clamp(a, b);
        }
    };

    /*
        Definition of the chaotic particle
        - `T` is the type used to store the coordinates (defaults to `float`)
        - `dim` is the number of dimensions for the vector (defaults to 2)

        - updatePosition: updates the position of the particle based on a chaotic map
    */
    template <typename T = float, int dim = 2>
    class ChaoticParticle : public Particle<T, dim>
    {
    private:
        const ChaosMap<T, dim> &chaosMap;

    public:
        ChaoticParticle(const ChaosMap<T, dim> &map) : Particle<T, dim>(), chaosMap(map) {}

        /*
            Description of the function updatePosition
            - parameters:
                - global_best: the best position found by the swarm
                - a: the minimum boundary for the position
                - b: the maximum boundary for the position
                - current_iteration: the current iteration of the swarm
                - max_iterations: the maximum number of iterations for the swarm

            the function updates the position of the chaotic particle based on a chaotic map
        */
        void updatePosition(const Point<T, dim> &global_best, const Point<T, dim> &a, const Point<T, dim> &b, int current_iteration, int max_iterations) override
        {
            this->position = chaosMap.getPoint(this->position, a, b, current_iteration);
        }
    };
};

#endif