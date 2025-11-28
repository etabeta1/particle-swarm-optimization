#ifndef _CHAOSMAP_HPP_
#define _CHAOSMAP_HPP_

#include "point.hpp"

#include <functional>

namespace Swarm
{
    template <typename T = float, typename U = float, int dim = 2>
    class ChaosMap
    {
    private:
        Point<T, dim> a;
        Point<T, dim> b;
        std::function<U(Point<T, dim> &)> map;
        inline Point<T, dim> mapBetweenDomains(const Point<T, dim> &);

    public:
        ChaosMap();
        inline Point<T, dim> toLocalDomain(const Point<T, dim> &);
        inline Point<T, dim> toGlobalDomain(const Point<T, dim> &);
    };
};

#endif