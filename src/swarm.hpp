#ifndef _SWARM_HPP_
#define _SWARM_HPP_
#include "point.hpp"
#include "chaosmap.hpp"
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
        std::Vector<Particle> particles;
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
        float c1;
        float c2;
        float w;
        float lBest_value;
        void updateSpeed();

    public:
        void updatelBest();
    };

    template <typename T = float,int dim = 2>
    class ChaoticParticle : public Particle<T, dim>
    {
        const ChaosMap &chaosMap;
    };
}

#endif
