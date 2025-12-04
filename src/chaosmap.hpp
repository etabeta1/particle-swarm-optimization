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
        Point<T, dim> a_min;
        Point<T, dim> b_max;
        std::function<Point<T,dim>(Point<T, dim> &)> map;
        inline Point<T, dim> mapBetweenDomains(const Point<T, dim> &point, const Point<T, dim> &first_range_start, const Point<T, dim> &first_range_end, const Point<T, dim> &min_range, const Point<T, dim> &max_range)
        {
            // codice per trasformare le coordinate
            return (point - first_range_start) / (first_range_end - first_range_start) * (max_range - min_range) + min_range;
        };
        inline Point<T, dim> generatePoint(const Point<T, dim> &point)
        {
            // map evaluate point map(Point),per ogni dimensione genera per ogni coordinata
            return map(point);
        };

    public:
        ChaosMap(std::function<U(Point<T, dim> &)> lambda_fun, Point<T, dim> min, Point<T, dim> max)
        {
            map = lambda_fun;
            a_min = min;
            b_max = max;
        };
        inline Point<T, dim> toLocalDomain(const Point<T, dim> &point, dim > &min_range, const Point<T, dim> &max_range)
        {
            return mapBetweenDomains(point, min_range, max_range, a_min, b_max);
        };
        inline Point<T, dim> toGlobalDomain(const Point<T, dim> &point, const Point<T, dim> &min_range, const Point<T, dim> &max_range)
        {
            return mapBetweenDomains(point, a_min, b_max, min_range, max_range);
        };
        inline Point<T, dim> getPoint(const Point<T, dim> &point, const Point<T, dim> &min_global, const Point<T, dim> &max_global)
        {
            Point<T, dim> new_point;

            new_point = toLocalDomain(point, min_global, max_global);
            new_point = generatePoint(new_point);
            new_point = toGlobalDomain(new_point, min_global, max_global);

            return new_point;
        };
    };
};

#endif