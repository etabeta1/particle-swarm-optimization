#ifndef _SWARM_HPP_
#define _SWARM_HPP_
#include "point.hpp"
#include "chaosmap.hpp"
#include "../lib/funcs.hpp"

#include <cmath>
#include <vector>

namespace Swarm
{
    // Class declarations
    template <typename T = float, int dim = 2>
    class Function
    {
    public:
        virtual float evaluate(const Point<T, dim>& p) const = 0;
        virtual ~Function() = default;
    };

    template <typename T = float, int dim = 2>
    class SphereFunction : public Function<T, dim>
    {
    public:
        float evaluate(const Point<T, dim>& p) const override {
            float sum = 0.0f;  
            for (int i = 0; i < dim; ++i) {
                sum += p[i] * p[i]; 
            }
            return sum;
        }
    };

    template <typename T = float, int dim = 2>
    class EllipsoidFunction : public Function<T, dim>
    {
    public:
        float evaluate(const Point<T, dim>& p) const override {
            float sum = 0.0f;  
            for (int i = 0; i < dim; ++i) {
                sum += i*p[i] * p[i]; 
            }
            return sum;
        }
    };

    template <typename T = float, int dim = 2>
    class QuinticFunction : public Function<T, dim>
    {
    public:
        float evaluate(const Point<T, dim>& p) const override {
            float sum = 0.0f;  
            for (int i = 0; i < dim; ++i) {
                sum += abs(pow(p[i],5) - 3*pow(p[i],4) + 4*pow(p[i],3) + 2*pow(p[i],2) -10*p[i] -4); 
            }
            return sum;
        }
    };

    template <typename T = float, int dim = 2>
    class DropwaveFunction : public Function<T, dim>
    {
    public:
        float evaluate(const Point<T, dim>& p) const override {
            float sum = 0.0f;
            float part_sum = 0.0f;
            for (int i = 0; i < dim; ++i) {
                part_sum = p[i]*p[i];
            }
            sum = 1 - (1 + cos(12 * sqrt(part_sum))) / (0.5 * part_sum + 2);
            return sum;
        }
    };

    template <typename T = float, int dim = 2>
    class Alpine1Function : public Function<T, dim>
    {
    public:
        float evaluate(const Point<T, dim>& p) const override {
            float sum = 0.0f;

            for (int i = 0; i < dim; ++i) {
                sum += abs(p[i] * sin(p[i]) + 0.1 * p[i]);
            }
            return sum;
        }
    };

    template <typename T = float, int dim = 2>
    class AckleyFunction : public Function<T, dim>
    {
    public:
        float evaluate(const Point<T, dim>& p) const override {
            float sum = 0.0f;
            float part_sum_1 = 0.0f;
            float part_sum_2 = 0.0f;
            for (int i = 0; i < dim; ++i) {
               part_sum_1 += p[i]*p[i];
                part_sum_2 += cos(2 * M_PI * p[i]); 
            }
            sum = -20 * exp(-0.2 * sqrt(part_sum_1 / dim)) - exp(part_sum_2 / dim) + 20 + exp(1);
            return sum;
        }
    };

    





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
    //protected:
        // float global

    public:
        Swarm(const Function<T, dim>& func) : fitness_function(func) {}

        void findGlobalBest();
        void updateEveryone();
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

        void updateSpeed(const Point<T, dim>& global_best, int current_itertions, int max_iterations);

    public:
        virtual void updatePosition(const Point<T, dim>& global_best,const Point<T, dim>& a, const Point<T, dim>& b,int current_iteration, int max_iterations);
        void updatePersonalBest(const Function<T, dim>& func, int current_iterations);
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
            particle->updatePosition(this->global_best, this->a, this->b, current_iteration,max_iterations);
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
    void NormalParticle<T, dim>::updateSpeed(const Point<T, dim>& global_best, int current_iterations, int max_iterations)
    {
        // CHOPSO velocity update formula
        // v(t+1) = w*v(t) + k1*(pBest - x(t)) + k2*(gBest - x(t))
        // k1 = rand[0,1] * 2.5
        // k2 = rand[0,1] * 2.5
        // w = 0.9 - (0.5 * iteration) / maxIterations) is the dynamic inertia weight

        float w = 0.9f - ((0.5f * static_cast<float>(current_iterations)) / static_cast<float>(max_iterations));

        float k1 = Funcs::randFloat(0.0f, 1.0f) * c1;
        float k2 = Funcs::randFloat(0.0f, 1.0f) * c2;

        Point<T, dim> inertia = speed * w;
        Point<T, dim> cognitive = (this->personal_best - this->position) * k1;
        Point<T, dim> social = (global_best - this->position) * k2;

        speed = inertia + cognitive + social;
    }

    template <typename T, int dim>
    void NormalParticle<T, dim>::updatePosition(const Point<T, dim>& global_best,const Point<T, dim>& a, const Point<T, dim>& b,int current_iteration, int max_iterations)
    {
        updateSpeed(global_best, current_iteration, max_iterations);
        this->position = this->position + this->speed; // clamp strategy is needed
    }

    template <typename T, int dim>
    void NormalParticle<T, dim>::updatePersonalBest(const Function<T, dim>& func, int current_iterations)
    {
        if(current_iterations==0){
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
    void ChaoticParticle<T, dim>::updatePosition(const Point<T, dim>& global_best,const Point<T, dim>& a, const Point<T, dim>& b,int current_iteration, int max_iterations)
    {
        this->position = chaosMap.getPoint(this->position, a, b);
    }
}
#endif
