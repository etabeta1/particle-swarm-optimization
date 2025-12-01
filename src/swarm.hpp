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

    public:
        void findGBest();
        void updateEveryone();
    };

    template <typename T = float, int dim = 2>
    class Particle
    {
    private:
        Point<T, dim> position;
        Point<T, dim> lBest;

    public:
        void updatePosition();
    };

    template <typename T = float, int dim = 2>
    class NormalParticle : Particle<T>
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
    void NormalParticle<T, dim>::updateSpeed(const Point<T, dim>& gBest, int it, int maxiter)
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
        Point<T, dim> cognitive = (this->lBest - this->position) * k1;
        Point<T, dim> social = (gBest - this->position) * k2;

        speed = inertia + cognitive + social;
    }
}
#endif
