#ifndef _SWARM_HPP_
#define _SWARM_HPP_
#include "point.hpp"
#include "chaosmap.hpp"
#include "../lib/funcs.hpp"

#include <vector>

namespace Swarm
{
    template <typename T = float, int dim = 2>
    class Function
    {
    public:
        float evaluate(Point<T, dim>& p) const;
        virtual ~Function() = default;
    };

    template <typename T = float, int dim = 2>
    class Swarm
    {
    private:
        std::vector<Particle> particles;
        Point<T, dim> gBest;
        float gBest_value;
        Point<T, dim> a;
        Point<T, dim> b;
        int current_iteration;
        int max_iterations;

    public:
        void findGBest();
        void updateEveryone();
    };

    template <typename T = float, int dim = 2>
    class Particle
    {
    protected:
        Point<T, dim> position;
        Point<T, dim> lBest;

    public:
        void updatePosition();
        virtual ~Particle() {}
    };

    template <typename T = float, int dim = 2>
    class NormalParticle : public Particle<T, dim>
    {
    private:
        Point<T, dim> speed;
        float c1 = 2.5f;
        float c2 = 2.5f;
        // float w; dynamic inertia weight
        float lBest_value;

        void updateSpeed(const Point<T, dim>& gBest, int it, int maxiter);

    public:
        void updatelBest();
    };

    template <typename T = float,int dim = 2>
    class ChaoticParticle : public Particle<T, dim>
    {
        const ChaosMap &chaosMap;
    };

    template <typename T, int dim>
    void NormalParticle<T, dim>::updateSpeed(const Point<T, dim>& global_best, int it, int maxiter)
    {
        // CHOPSO velocity update formula
        // v(t+1) = w*v(t) + k1*(pBest - x(t)) + k2*(gBest - x(t))
        // k1 = rand[0,1] * 2.5
        // k2 = rand[0,1] * 2.5
        // w = 0.9 - (0.5 * iteration) / maxIterations) is the dynamic inertia weight

        float w = 0.9f - ((0.5f * static_cast<float>(it)) / static_cast<float>(maxiter));

        float k1 = Funcs::randFloat(0.0f, 1.0f) * c1;
        float k2 = Funcs::randFloat(0.0f, 1.0f) * c2;

        Point<T, dim> inertia = speed * w;
        Point<T, dim> cognitive = (this->personal_best - this->position) * k1;
        Point<T, dim> social = (global_best - this->position) * k2;

        speed = inertia + cognitive + social;
    }

    template <typename T, int dim>
    void Swarm<T, dim>::findGbest()
    {
        for (const auto& particle : particles)
        {
            float fitness = /* evaluate fitness of particle.position */;
            if (fitness < gBest_value)
            {
                gBest_value = fitness;
                gBest = particle.position;
            }
        }
    }

    template <typename T, int dim>
    void Swarm<T, dim>::updateEveryone(){
        for(auto& particle : particles){
            particle.updateSpeed(this->gBest, current_iteration, max_iterations);
            particle.updatePosition();
            particle.updatelBest(); //maybe is useful to update lBest directly in the particle updatePosition method
        }
    }

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

    
}
#endif
