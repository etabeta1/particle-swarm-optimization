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

        EvaluatedPoint<T, dim> personal_best;

    public:
        /**
         * \brief Constructs a `Particle` with default position and personal best.
         *
         * The position and personal best are initialized to zero.
         */
        Particle() : position(T(0)), personal_best(T(0), T(0)) {}

        /**
         * \brief Initializes the particle as the optimizer likes.
         * \param initial_position The initial position of the particle.
         * \param func The fitness function to evaluate the particle's position.
         *
         * The function sets the position, personal best, and personal best value of the particle.
         * The idea is to split the instantiation from the initialization so that the same class can be easily reused for multiple different optimizers that may need to initialize/reinitialize it at different stages and/or with different parameters.
         */
        void reinit(const Point<T, dim> &initial_position, const ObjectiveFunction<T, dim> &func)
        {
            reinit(initial_position, initial_position, func);
        }

        /**
         * \brief Initializes the particle as the optimizer likes.
         * \param initial_position The initial position of the particle.
         * \param local_best The personal best position of the particle.
         * \param func The fitness function to evaluate the particle's position.
         *
         * The function sets the position, personal best, and personal best value of the particle.
         * The idea is to split the instantiation from the initialization so that the same class can
         * be easily reused for multiple different optimizers that may need to
         * initialize/reinitialize it at different stages and/or with different parameters.
         *
         * `initial_position` and `local_best` can be different to allow more complex initialization
         * strategies where the initial position does not necessarily coincide with the initial
         * personal best.
         */
        virtual void reinit(const Point<T, dim> &initial_position, const Point<T, dim> &local_best, const ObjectiveFunction<T, dim> &func)
        {
            this->position = initial_position;
            this->personal_best = {local_best, func.evaluate(local_best)};
        }

        /**
         * \brief Updates the position of the particle.
         * \param global_best The best position found by the swarm.
         * \param a The minimum boundary for the position.
         * \param b The maximum boundary for the position.
         * \param current_iteration The current iteration of the swarm.
         * \param max_iterations The maximum number of iterations for the swarm.
         * \param constraints A vector of constraints that the particle position must satisfy.
         * \see Constraint
         */
        virtual void updatePosition(const Point<T, dim> &global_best, const Point<T, dim> &a, const Point<T, dim> &b, IterationType current_iteration, IterationType max_iterations, const std::vector<Constraint<T, dim>> &constraints) = 0;

        /**
         * \brief Returns the type of the particle.
         * \return An integer representing the type of the particle.
         *
         * The type can be used to distinguish between different particle behaviors.
         * This method is just used to save the type of the particle in the log files.
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
         * \param constraints A vector of constraints that the particle position must satisfy.
         * \return `true` if the personal best was updated, `false` otherwise.
         *
         * \warning This method calls `handleOutOfBounds` **instead** of actually performing the update if the particle happens to be out of bounds.
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

            // Minimization
            if (current_value < this->personal_best.value)
            {
                this->personal_best = {this->position, current_value};

                return true;
            }

            return false;
        }

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
        Point<T, dim> &getPersonalBestPosition()
        {
            return personal_best.point;
        }

        /**
         * \brief Returns the personal best value of the particle.
         * \return A float representing the particle's personal best value.
         */
        float getPersonalBestValue()
        {
            return personal_best.value;
        }

    protected:
        /**
         * \brief Called when the particle is out of bounds.
         * \param first_failed_constraint The first constraint that the particle failed to satisfy.
         * \param all_constraints A vector of all constraints.
         * \return `true` if the personal best was updated, `false` otherwise.
         * \see Constraint
         * \see updatePersonalBest
         *
         * This function is called when the particle's position violates any of the constraints. It allows derived classes to define specific behaviors for handling out-of-bounds situations.
         *
         * \warning This method is called by `updatePersonalBest` when the particle is out of bounds **instead** of actually performing the update.
         */
        virtual bool handleOutOfBounds(const Constraint<T, dim> &first_failed_constraint, const std::vector<Constraint<T, dim>> &all_constraints) = 0;
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

        float r1, r2;

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

            // float k1 = generate_random(0.0f, 1.0f) * c1;
            // float k2 = generate_random(0.0f, 1.0f) * c2;

            float k1 = c1 * r1;
            float k2 = c2 * r2;

            Point<T, dim> inertia = speed * w;
            Point<T, dim> cognitive = (this->personal_best.point - this->position) * k1;
            Point<T, dim> social = (global_best - this->position) * k2;

            this->speed = inertia + cognitive + social;
        }

    public:
        /**
         * \brief Constructs a `NormalParticle` with default position and speed.
         *
         * The position and speed are initialized to zero.
         */
        NormalParticle() : Particle<T, dim>(), speed(T(0)), r1(generate_random(0.f, 1.f)), r2(generate_random(0.f, 1.f)) {}

        /**
         * \copydoc Particle::updatePosition
         *
         * This function updates the position of the normal particle based on its speed and clamps it within the boundaries.
         */
        void updatePosition(const Point<T, dim> &global_best, const Point<T, dim> &a, const Point<T, dim> &b, IterationType current_iteration, IterationType max_iterations, const std::vector<Constraint<T, dim>> &) override
        {
            updateSpeed(global_best, current_iteration, max_iterations);

            this->position = (this->position + this->speed).clamp(a, b);
        }

        /**
         * \copydoc Particle::getType
         */
        int getType() const override
        {
            return 0; // Normal particle type
        }

    protected:
        /**
         * \copydoc Particle::handleOutOfBounds
         */
        bool handleOutOfBounds(const Constraint<T, dim> &, const std::vector<Constraint<T, dim>> &) override
        {
            // If we are OOB, just do not update the personal best, the velocity will eventually point towards a valid point
            return false;
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
         * \copydoc Particle::updatePosition
         *
         * This function updates the position of the chaotic particle based on the chaos map.
         */
        void updatePosition(const Point<T, dim> &, const Point<T, dim> &a, const Point<T, dim> &b, IterationType current_iteration, IterationType, const std::vector<Constraint<T, dim>> &) override
        {
            this->position = chaosMap.getPoint(this->position, a, b, current_iteration);
        }

        /**
         * \copydoc Particle::getType
         */
        int getType() const override
        {
            return 1; // Chaotic particle type
        }

    protected:
        /**
         * \copydoc Particle::handleOutOfBounds
         */
        bool handleOutOfBounds(const Constraint<T, dim> &, const std::vector<Constraint<T, dim>> &) override
        {
            // If we are OOB, just do not update the personal best, the velocity will eventually point towards a valid point
            return false;
        }
    };

    /**
     * \brief Definition of the fuzzy particle for FST-PSO algorithm.
     * \tparam T The type used to store the coordinates (defaults to float).
     * \tparam dim The number of dimensions for the vector (defaults to 2).
     *
     * Each particle independently adjusts its inertia, cognitive factor, social factor,
     * and velocity bounds at each iteration using fuzzy inference based on:
     * - delta: distance from the global best position
     * - phi: normalized fitness incremental factor
     */
    template <typename T = float, int dim = 2>
    class FuzzyParticle : public Particle<T, dim>
    {
    private:
        /**
         * \brief The velocity of the particle.
         */
        Point<T, dim> speed;

        /**
         * \brief Previous evaluated position for phi calculation.
         *
         * Stores both the position and fitness value from the previous iteration.
         * This is used to compute the normalized fitness incremental factor (phi).
         *
         * \see FSTPSOOptimizer::computePhi
         */
        EvaluatedPoint<T, dim> previous;

        /**
         * \brief Individual inertia weight (calculated by fuzzy inference).
         *
         * Controls how much the particle "remembers" its previous velocity.
         * Calculated using Rules 1-3:
         * - Rule 1: if (phi is Worse OR delta is Same) then (Inertia is Low = 0.3)
         * - Rule 2: if (phi is Same OR delta is Near) then (Inertia is Medium = 0.5)
         * - Rule 3: if (phi is Better OR delta is Far) then (Inertia is High = 1.0)
         *
         * Defuzzified range: [0.3, 1.0]
         */
        T w;

        /**
         * \brief Individual cognitive factor (calculated by fuzzy inference).
         *
         * Controls attraction toward the particle's personal best position.
         * Calculated using Rules 7-9:
         * - Rule 7: if (delta is Far) then (Cognitive is Low = 0.1)
         * - Rule 8: if (phi is Worse OR phi is Same OR delta is Same OR delta is Near)
         *           then (Cognitive is Medium = 1.5)
         * - Rule 9: if (phi is Better) then (Cognitive is High = 3.0)
         *
         * Defuzzified range: [0.1, 3.0]
         */
        T c1;

        /**
         * \brief Individual social factor (calculated by fuzzy inference).
         *
         * Controls attraction toward the global best position.
         * Calculated using Rules 4-6:
         * - Rule 4: if (phi is Better OR delta is Near) then (Social is Low = 1.0)
         * - Rule 5: if (phi is Same OR delta is Same) then (Social is Medium = 2.0)
         * - Rule 6: if (phi is Worse OR delta is Far) then (Social is High = 3.0)
         *
         * Defuzzified range: [1.0, 3.0]
         */
        T c2;

        /**
         * \brief Lower velocity clamp factor (calculated by fuzzy inference).
         *
         * Multiplied by the search space range to obtain the minimum velocity per dimension:
         * v_min[i] = L × (b[i] - a[i])
         *
         * Calculated using Rules 10-12:
         * - Rule 10: if (phi is Same OR phi is Better OR delta is Far) then (L is Low = 0.0)
         * - Rule 11: if (delta is Same OR delta is Near) then (L is Medium = 0.001)
         * - Rule 12: if (phi is Worse) then (L is High = 0.01)
         *
         * Defuzzified range: [0.0, 0.01]
         */
        T L;

        /**
         * \brief Upper velocity clamp factor (calculated by fuzzy inference).
         *
         * Multiplied by the search space range to obtain the maximum velocity per dimension:
         * v_max[i] = U × (b[i] - a[i])
         *
         * Calculated using Rules 13-15:
         * - Rule 13: if (delta is Same) then (U is Low = 0.1)
         * - Rule 14: if (phi is Same OR phi is Better OR delta is Near) then (U is Medium = 0.15)
         * - Rule 15: if (phi is Worse OR delta is Far) then (U is High = 0.2)
         *
         * Defuzzified range: [0.1, 0.2]
         */
        T U;

        float r1, r2;

        /**
         * \brief Updates the velocity of the fuzzy particle.
         * \param global_best The best position found by the swarm.
         * \param current_iteration Current iteration (unused in FST-PSO).
         * \param max_iterations Maximum iterations (unused in FST-PSO).
         *
         * Implements the FST-PSO velocity update formula
         */
        void updateSpeed(const Point<T, dim> &global_best, IterationType current_iteration, IterationType max_iterations)
        {
            float k1 = c1 * r1;
            float k2 = c2 * r2;

            Point<T, dim> inertia = speed * w;
            Point<T, dim> cognitive = (this->personal_best.point - this->position) * k1;
            Point<T, dim> social = (global_best - this->position) * k2;

            this->speed = inertia + cognitive + social;
        }

    public:
        /**
         * \brief Constructs a `FuzzyParticle` with default Medium values.
         *
         * Initializes all fuzzy parameters to their "Medium" linguistic values:
         *
         * These values will be dynamically updated by fuzzy inference starting
         * from the first iteration via setFuzzyParameters().
         */
        FuzzyParticle()
            : Particle<T, dim>(),
              speed(T(0)),
              previous({Point<T, dim>(T(0)), T(0)}),
              w(0.5),
              c1(1.5),
              c2(2.0),
              L(0.001),
              U(0.15),
              r1(generate_random(0.f, 1.f)),
              r2(generate_random(0.f, 1.f)) {}

        /**
         * \brief Initializes the fuzzy particle with real initial position and fitness.
         * \param initial_position The initial position of the particle in the search space.
         * \param local_best The personal best position (typically same as initial_position).
         * \param func The fitness function to evaluate the particle's position.
         *
         * This override ensures that the previous evaluation is initialized with the
         * real initial position and fitness value. This is critical
         * for correct phi calculation starting from the first iteration.
         */
        void reinit(const Point<T, dim> &initial_position, const Point<T, dim> &local_best, const ObjectiveFunction<T, dim> &func) override
        {
            // Call base class reinit to set position and personal_best
            Particle<T, dim>::reinit(initial_position, local_best, func);

            // Initialize previous with REAL initial position and fitness
            previous = {initial_position, func.evaluate(initial_position)};
        }

        /**
         * \brief Updates the position of the fuzzy particle with velocity clamping.
         * \param global_best The best position found by the swarm.
         * \param a The minimum boundary for the position.
         * \param b The maximum boundary for the position.
         * \param current_iteration Current iteration (unused in FST-PSO).
         * \param max_iterations Maximum iterations (unused in FST-PSO).
         * \param constraints Constraints (unused in this implementation).
         */
        void updatePosition(const Point<T, dim> &global_best, const Point<T, dim> &a, const Point<T, dim> &b, IterationType current_iteration, IterationType max_iterations, const std::vector<Constraint<T, dim>> &) override
        {
            // Save current position for next iteration's phi calculation
            previous.point = this->position;

            // Update velocity using fuzzy-calculated parameters
            updateSpeed(global_best, current_iteration, max_iterations);

            // Calculate velocity bounds based on fuzzy L and U factors
            Point<T, dim> v_min = (b - a) * L;
            Point<T, dim> v_max = (b - a) * U;

            // Clamp velocity component-wise, preserving sign
            this->speed = Point<T, dim>([this, &v_min, &v_max](size_t i)
                                        {
                T vel = this->speed[i];
                T abs_vel = std::abs(vel);

                // If velocity magnitude is below minimum, set to minimum
                if (abs_vel < v_min[i])
                    return (vel < 0) ? -v_min[i] : v_min[i];

                // If velocity magnitude exceeds maximum, clamp to maximum
                if (abs_vel > v_max[i])
                    return (vel < 0) ? -v_max[i] : v_max[i];

                return vel; });

            // Update position and clamp to search space boundaries
            this->position = (this->position + this->speed).clamp(a, b);
        }

        /**
         * \copydoc Particle::getType
         *
         * \return 2 for FuzzyParticle (0=Normal, 1=Chaotic, 2=Fuzzy)
         */
        int getType() const override
        {
            return 2; // Fuzzy particle type
        }

        /**
         * \brief Set individual fuzzy parameters calculated by FSTPSOOptimizer.
         * \param _w Inertia weight (defuzzified from Rules 1-3).
         * \param _c_soc Social factor (defuzzified from Rules 4-6).
         * \param _c_cog Cognitive factor (defuzzified from Rules 7-9).
         * \param _L Lower velocity clamp factor (defuzzified from Rules 10-12).
         * \param _U Upper velocity clamp factor (defuzzified from Rules 13-15).
         *
         * These parameters are used immediately in the next updatePosition() call
         * to calculate velocity and position updates.
         */
        void setFuzzyParameters(T _w, T _c_soc, T _c_cog, T _L, T _U)
        {
            w = _w;
            c2 = _c_soc;
            c1 = _c_cog;
            L = _L;
            U = _U;
        }

        /**
         * \brief Get previous evaluated position for phi calculation.
         * \return Const reference to previous evaluation containing position and fitness.
         */
        const EvaluatedPoint<T, dim> &getPreviousEval() const
        {
            return previous;
        }

        /**
         * \brief Update the fitness value of the previous evaluation.
         * \param fitness The fitness value to store for the current position.
         *
         * Called by FSTPSOOptimizer after each fitness evaluation to keep
         * previous.value synchronized with the actual fitness of the current position.
         */
        void updatePreviousFitness(float fitness)
        {
            previous.value = fitness;
        }

    protected:
        /**
         * \copydoc Particle::handleOutOfBounds
         */
        bool handleOutOfBounds(const Constraint<T, dim> &, const std::vector<Constraint<T, dim>> &) override
        {
            return false;
        }
    };
};

#endif