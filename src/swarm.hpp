#ifndef _SWARM_HPP_
#define _SWARM_HPP_
#include "point.hpp"
#include "chaosmap.hpp"
#include "function.hpp"
#include "../lib/funcs.hpp"

#include <cmath>
#include <vector>
#include <memory>

namespace Swarm
{
    //Class Declarations
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
        std::vector<std::unique_ptr<Particle<T, dim>>> particles; //stores both normal and chaotic particles
        const Function<T, dim>& fitness_function;

        Point<T, dim> global_best;
        float global_best_value;

        Point<T, dim> a;
        Point<T, dim> b;

        int current_iteration;
        int max_iterations;

    public:
        Swarm(const Function<T, dim>& func) : fitness_function(func) {}

        void addParticle(std::unique_ptr<Particle<T, dim>> p) {
            particles.emplace_back(std::move(p));
    }

        void findGlobalBest();
        void updateEveryone();

        Point<T, dim> getGlobalBest() const {
            return global_best;
        }

        float getGlobalBestValue() const {
            return global_best_value;
        }
    };

    /*
        Definition of a generic particle in the swarm.
        - `T` is the type used to store the coordinates (defaults to `float`)
        - `dim` is the number of dimensions for the vector (defaults to 2)

        - updatePosition: updates the position of the particle
    */
    template <typename T = float, int dim = 2>
    class Particle
    {
    protected:
        Point<T, dim> position;
    public:
        virtual void updatePosition(const Point<T, dim>& global_best,const Point<T, dim>& a, const Point<T, dim>& b,int current_iteration, int max_iterations) = 0;
        virtual ~Particle() {}
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
        
        Point<T, dim> personal_best;
        float personal_best_value; // only normal particle can store the personal best

        void updateSpeed(const Point<T, dim>& global_best, int current_iteration, int max_iterations);

    public:
        NormalParticle(const Point<T, dim>& position_, const Point<T, dim>& speed_, const Function<T, dim>& func) : speed(speed_), personal_best(position_), personal_best_value(func.evaluate(position_)) {
            this->position = position_;
        } 

        void updatePosition(const Point<T, dim>& global_best,const Point<T, dim>& a, const Point<T, dim>& b,int current_iteration, int max_iterations) override;
        void updatePersonalBest(const Function<T, dim>& func, int current_iteration);
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
        const ChaosMap<T, float, dim>& chaosMap;
    public:
        ChaoticParticle(const ChaosMap<T, float, dim>& map) : chaosMap(map) {}
    
        void updatePosition(const Point<T, dim>& global_best,const Point<T, dim>& a, const Point<T, dim>& b,int current_iteration, int max_iterations) override;
    };

    //Class functions
    /*
        Description of the function findGlobalBest
        - parameters: none

        the function iterates through all particles to find the one with the best fitness value
    */
    template <typename T, int dim>
    void Swarm<T, dim>::findGlobalBest()
    {
        for (const auto& particle : particles)
        {
            float fitness = fitness_function.evaluate(particle->position);
            if (fitness < global_best_value)
            {
                global_best_value = fitness;
                global_best = particle->position;
            }
        }
    }
    /*
        Description of the function updateEveryone
        - parameters: none

        the function updates the position of all particles, updates the personal best for normal particles,
        finds the global best, and increments the current iteration
    */
    template <typename T, int dim>
    void Swarm<T, dim>::updateEveryone(){
        for(auto& particle : particles){
            particle->updatePosition(this->global_best, this->a, this->b, current_iteration, max_iterations);

            if (auto* normalParticle = dynamic_cast<NormalParticle<T, dim>*>(particle.get())) {
                normalParticle->updatePersonalBest(fitness_function, current_iteration);
            }
        }
        
        findGlobalBest();

        ++current_iteration;
    }
    /*
        Description of the function updateSpeed
        - parameters:
            - global_best: the best position found by the swarm
            - current_iteration: the current iteration of the swarm
            - max_iterations: the maximum number of iterations for the swarm

        the function updates the speed of the normal particle based on the CHOPSO velocity update formula
    */
    template <typename T, int dim>
    void NormalParticle<T, dim>::updateSpeed(const Point<T, dim>& global_best, int current_iteration, int max_iterations)
    {
        // CHOPSO velocity update formula
        // v(t+1) = w*v(t) + k1*(pBest - x(t)) + k2*(gBest - x(t))
        // k1 = rand[0,1] * 2.5
        // k2 = rand[0,1] * 2.5
        // w = 0.9 - (0.5 * iteration) / maxIterations) is the dynamic inertia weight

        float w = 0.9f - ((0.5f * static_cast<float>(current_iteration)) / static_cast<float>(max_iterations));

        float k1 = Funcs::randFloat(0.0f, 1.0f) * c1;
        float k2 = Funcs::randFloat(0.0f, 1.0f) * c2;

        Point<T, dim> inertia = speed * w;
        Point<T, dim> cognitive = (this->personal_best - this->position) * k1;
        Point<T, dim> social = (global_best - this->position) * k2;

        this->speed = inertia + cognitive + social;
    }

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
    template <typename T, int dim>
    void NormalParticle<T, dim>::updatePosition(const Point<T, dim>& global_best,const Point<T, dim>& a, const Point<T, dim>& b,int current_iteration, int max_iterations)
    {
        updateSpeed(global_best, current_iteration, max_iterations);

        this->position = this->position + this->speed;

        this-> position = this->position.clamp(a,b);
    }
    /*
        Description of the function updatePersonalBest
        - parameters:
            - func: the fitness function to evaluate the particle's position
            - current_iteration: the current iteration of the swarm

        the function updates the personal best position and value of the normal particle if the current position is better
    */
    template <typename T, int dim>
    void NormalParticle<T, dim>::updatePersonalBest(const Function<T, dim>& func, int current_iteration)
    {
        float current_value = func.evaluate(this->position);
        if(current_value < this->personal_best_value){
            this->personal_best_value = current_value;
            this->personal_best = this->position;
        }
    }
    
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
    template <typename T, int dim>
    void ChaoticParticle<T, dim>::updatePosition(const Point<T, dim>& global_best,const Point<T, dim>& a, const Point<T, dim>& b,int current_iteration, int max_iterations)
    {
        this->position = chaosMap.getPoint(this->position, a, b);
    }
}
#endif
