#ifndef _SWARM_HPP_
#define _SWARM_HPP_
#include "point.hpp"
#include "chaosmap.hpp"
#include "function.hpp"
#include "../lib/funcs.hpp"

#include <cmath>
#include <vector>

namespace Swarm
{
    //Class Declaration
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

        void findGlobalBest();
        void updateEveryone();

        Point<T, dim> getGlobalBest() const {
            return global_best;
        }

        float getGlobalBestValue() const {
            return global_best_value;
        }
    };

    template <typename T = float, int dim = 2>
    class Particle
    {
    protected:
        Point<T, dim> position;
        Point<T, dim> personal_best;

    public:
        virtual void updatePosition(const Point<T, dim>& global_best,const Point<T, dim>& a, const Point<T, dim>& b,int current_iteration, int max_iterations) = 0;
        virtual ~Particle() {}
    };

    template <typename T = float, int dim = 2>
    class NormalParticle : public Particle<T, dim>
    {
    private:
        Point<T, dim> speed;
        float c1 = 2.5f;
        float c2 = 2.5f;
        float personal_best_value; // solo la normal particle

        void updateSpeed(const Point<T, dim>& global_best, int current_itertion, int max_iterations);

    public:
        void updatePosition(const Point<T, dim>& global_best,const Point<T, dim>& a, const Point<T, dim>& b,int current_iteration, int max_iterations) override;
        void updatePersonalBest(const Function<T, dim>& func, int current_iteration);
    };

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

    template <typename T, int dim>
    void Swarm<T, dim>::updateEveryone(){
        for(auto& particle : particles){
            particle->updatePosition(this->global_best, this->a, this->b, current_iteration, max_iterations);
            // allora praticamente faccio tutto in update position
            // chaotic semplicemente calcola la posizione
            // normal calcola la velocità, update posizione, e update personal best
            // non è importante il personal best per il chaotic
            // dio caro ho sonno
            particle->updatePersonalBest(fitness_function, current_iteration); // da sistemare per chaotic
        }
        ++current_iteration;
    }

    template <typename T, int dim>
    void NormalParticle<T, dim>::updateSpeed(const Swarm& swarm, int current_iteration, int max_iterations)
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
        Point<T, dim> social = (swarm.getGlobalBest() - this->position) * k2;

        speed = inertia + cognitive + social;
    }

    template <typename T, int dim>
    void NormalParticle<T, dim>::updatePosition(const Swarm& swarm,const Point<T, dim>& a, const Point<T, dim>& b,int current_iteration, int max_iterations)
    {
        updateSpeed(swarm, current_iteration, max_iterations);
        this->position = this->position + this->speed; // clamp strategy is needed
    }

    template <typename T, int dim>
    void NormalParticle<T, dim>::updatePersonalBest(const Function<T, dim>& func, int current_iteration)
    {
        if(current_iteration==0){
            this->personal_best_value = func.evaluate(this->position);
            this->personal_best = this->position;
        }
        else{
            float current_value = func.evaluate(this->position);
            if(current_value < this->personal_best_value){
                this->personal_best_value = current_value;
                this->personal_best = this->position;
            }
        }
    }

    template <typename T, int dim>
    void ChaoticParticle<T, dim>::updatePosition(const Swarm& swarm,const Point<T, dim>& a, const Point<T, dim>& b,int current_iteration, int max_iterations)
    {
        this->position = chaosMap.getPoint(this->position, a, b);
    }
}
#endif
