#ifndef _PARTICLE_HPP_
#define _PARTICLE_HPP_

#include "point.hpp"
#include "funcs.hpp"

namespace Swarm
{
    /**
     * \brief Abstract base class representing a particle in the swarm.
     * \tparam T The type used to store the coordinates (defaults to `float`)
     * \tparam dim The number of dimensions for the vector (defaults to 2)
     *
     * This class serves as a base for different types of particles in the swarm optimization algorithm.
     */
    template <typename T, int dim>
    class Particle
    {
    protected:
        /**
         * \brief The current position of the particle.
         */
        Point<T, dim> position;

        /**
         * \brief The personal best position of the particle.
         */
        Point<T, dim> personal_best;

        /**
         * \brief The personal best value of the particle.
         */
        float personal_best_value;

    public:
        /**
         * \brief Constructs a `Particle` with default position and personal best.
         * The position and personal best are initialized to zero.
         */
        Particle() : position(T(0)), personal_best(T(0)), personal_best_value(T(0)) {}

        /**
         * \brief Initializes the particle as the CHOPSOOptimizer likes.
         * \param initial_position The initial position of the particle.
         * \param func The fitness function to evaluate the particle's position.
         *
         * The function sets the position, personal best, and personal best value of the particle.
         */
        void reinit(const Point<T, dim> &initial_position, const ObjectiveFunction<T, dim> &func)
        {
            position = initial_position;
            personal_best = initial_position;
            personal_best_value = func.evaluate(initial_position);
        }

        /**
         * \brief Updates the position of the particle.
         * \param global_best The best position found by the swarm.
         * \param a The minimum boundary for the position.
         * \param b The maximum boundary for the position.
         * \param current_iteration The current iteration of the swarm.
         * \param max_iterations The maximum number of iterations for the swarm.
         */
        virtual void updatePosition(const Point<T, dim> &global_best, const Point<T, dim> &a, const Point<T, dim> &b, IterationType current_iteration, IterationType max_iterations, const std::vector<Constraint<T, dim>> &constraints) = 0;

        /**
         * \brief Returns the type of the particle.
         * \return An integer representing the type of the particle.
         * The type can be used to distinguish between different particle behaviors.
         */
        virtual int getType() const = 0;

        /**
         * \brief Virtual destructor for the `Particle` class.
         *
         * Ensures proper cleanup of derived classes.
         */
        virtual ~Particle() {}

        /**
         * \brief Updates the personal best position and value of the particle.
         * \param func The fitness function to evaluate the particle's position.
         * \return `true` if the personal best was updated, `false` otherwise.
         */
        virtual bool updatePersonalBest(const ObjectiveFunction<T, dim> &func, const std::vector<Constraint<T, dim>> &constraints)
        {
            // If the particle does not respect the constraint
            for (auto &c : constraints)
            {
                if (!(c(this->position)))
                {
                    return handleOutOfBounds(c, constraints);
                }
            }

            float current_value = func.evaluate(this->position);
            if (current_value < this->personal_best_value)
            {
                this->personal_best_value = current_value;
                this->personal_best = this->position;

                return true;
            }

            return false;
        }

        virtual bool handleOutOfBounds(const Constraint<T, dim> &first_failed_constraint, const std::vector<Constraint<T, dim>> &all_constraints) = 0;

        /**
         * \brief Returns a reference to the current position of the particle.
         * \return A reference to the `Point` representing the particle's position.
         * \see Point
         */
        Point<T, dim> &getPosition()
        {
            return position;
        }

        /**
         * \brief Returns a reference to the personal best position of the particle.
         * \return A reference to the `Point` representing the particle's personal best position.
         * \see Point
         */
        Point<T, dim> &getPersonalBest()
        {
            return personal_best;
        }

        /**
         * \brief Returns the personal best value of the particle.
         * \return A float representing the particle's personal best value.
         */
        float getPersonalBestValue()
        {
            return personal_best_value;
        }
    };

    /**
     * \brief Definition of the normal particle.
     * \tparam T The type used to store the coordinates (defaults to `float`).
     * \tparam dim The number of dimensions for the vector (defaults to 2).
     *
     * This class represents a normal particle that follows the PSO logic.
     */
    template <typename T = float, int dim = 2>
    class NormalParticle : public Particle<T, dim>
    {
    private:
        /**
         * \brief The speed of the particle.
         */
        Point<T, dim> speed;

        /**
         * \brief Cognitive coefficient.
         */
        float c1 = 2.5f;

        /**
         * \brief Social coefficient.
         */
        float c2 = 2.5f;

        /**
         * \brief Updates the speed of the normal particle.
         * \param global_best The best position found by the swarm.
         * \param current_iteration The current iteration of the swarm.
         * \param max_iterations The maximum number of iterations for the swarm.
         *
         * This function updates the speed of the normal particle based on the CHOPSO velocity update formula.
         */
        void updateSpeed(const Point<T, dim> &global_best, IterationType current_iteration, IterationType max_iterations)
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
        /**
         * \brief Constructs a `NormalParticle` with default position and speed.
         *
         * The position and speed are initialized to zero.
         */
        NormalParticle() : Particle<T, dim>(), speed(T(0)) {}

        /**
         * \brief Updates the position of the normal particle.
         * \param global_best The best position found by the swarm.
         * \param a The minimum boundary for the position.
         * \param b The maximum boundary for the position.
         * \param current_iteration The current iteration of the swarm.
         * \param max_iterations The maximum number of iterations for the swarm.
         *
         * This function updates the position of the normal particle based on its speed and clamps it within the boundaries.
         */
        void updatePosition(const Point<T, dim> &global_best, const Point<T, dim> &a, const Point<T, dim> &b, IterationType current_iteration, IterationType max_iterations, const std::vector<Constraint<T, dim>> &) override
        {
            updateSpeed(global_best, current_iteration, max_iterations);

            this->position = (this->position + this->speed).clamp(a, b);
        }

        bool handleOutOfBounds(const Constraint<T, dim> &, const std::vector<Constraint<T, dim>> &)
        {
            // If we are OOB, just do not update the personal best, the velocity will eventually point towards a valid point
            return false;
        }

        /**
         * \brief Returns the type of the normal particle.
         * \return An integer representing the type of the particle (0 for normal particles).
         *
         * The type can be used to distinguish between different particle behaviors.
         */
        int getType() const override
        {
            return 0; // Normal particle type
        }
    };

    /**
     * \brief Definition of the chaotic particle.
     * \tparam T The type used to store the coordinates (defaults to `float`)
     * \tparam dim The number of dimensions for the vector (defaults to 2)
     *
     * This class represents a chaotic particle that updates its position based on a chaotic map.
     */
    template <typename T = float, int dim = 2>
    class ChaoticParticle : public Particle<T, dim>
    {
    private:
        /**
         * \brief The chaos map used to update the particle's position.
         */
        const ChaosMap<T, dim> &chaosMap;

    public:
        /**
         * \brief Constructs a `ChaoticParticle` with the given chaos map.
         * \param map The chaos map used to update the particle's position.
         */
        ChaoticParticle(const ChaosMap<T, dim> &map) : Particle<T, dim>(), chaosMap(map) {}

        /**
         * \brief Updates the position of the chaotic particle.
         * \param global_best The best position found by the swarm.
         * \param a The minimum boundary for the position.
         * \param b The maximum boundary for the position.
         * \param current_iteration The current iteration of the swarm.
         * \param max_iterations The maximum number of iterations for the swarm.
         *
         * This function updates the position of the chaotic particle based on the chaos map.
         */
        void updatePosition(const Point<T, dim> &, const Point<T, dim> &a, const Point<T, dim> &b, IterationType current_iteration, IterationType max_iterations, const std::vector<Constraint<T, dim>> &) override
        {
            this->position = chaosMap.getPoint(this->position, a, b, current_iteration);
        }

        bool handleOutOfBounds(const Constraint<T, dim> &, const std::vector<Constraint<T, dim>> &)
        {
            // If we are OOB, just do not update the personal best, the velocity will eventually point towards a valid point
            return false;
        }

        /**
         * \brief Returns the type of the chaotic particle.
         * \return An integer representing the type of the particle (1 for chaotic particles).
         */
        int getType() const override
        {
            return 1; // Chaotic particle type
        }
    };
};

#endif