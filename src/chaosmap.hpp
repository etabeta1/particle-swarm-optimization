#ifndef _CHAOSMAP_HPP_
#define _CHAOSMAP_HPP_

#include "point.hpp"

#include <functional>
// Chaosmap is used to generate Random coordinates for Chaotic Particles.
namespace Swarm

{
    template <typename T = float, int dim = 2>
    using MapFunction = std::function<Point<T, dim>(const Point<T, dim> &, int)>;

    template <typename T = float, int dim = 2>
    class ChaosMap
    {
    private:
        MapFunction<T, dim> map;

        Point<T, dim> a_min;
        Point<T, dim> b_max;

        // Esegue la trasformazione delle coordiante del punto in un'altro dominio
        inline Point<T, dim> mapBetweenDomains(const Point<T, dim> &point, const Point<T, dim> &first_range_start, const Point<T, dim> &first_range_end, const Point<T, dim> &min_range, const Point<T, dim> &max_range) const
        {
            // codice per trasformare le coordinate
            return (point - first_range_start) / (first_range_end - first_range_start) * (max_range - min_range) + min_range;
        };
        // genera un punto secondo la lambda function
        inline Point<T, dim> generatePoint(const Point<T, dim> &point, int iter) const
        {
            // map evaluate point map(Point),per ogni dimensione genera per ogni coordinata
            return map(point, iter);
        };

    public:
        // ChaosMap Construct
        ChaosMap(const MapFunction<T, dim> &lambda_fun, const Point<T, dim> &min, const Point<T, dim> &max) : map(lambda_fun), a_min(min), b_max(max) {};
        // change coordinates of the Point from Global domain to local domain
        inline Point<T, dim> toLocalDomain(const Point<T, dim> &point, const Point<T, dim> &min_range, const Point<T, dim> &max_range) const
        {
            return mapBetweenDomains(point, min_range, max_range, a_min, b_max);
        };
        // change coordinates of the Point from Local domain to Global domain
        inline Point<T, dim> toGlobalDomain(const Point<T, dim> &point, const Point<T, dim> &min_range, const Point<T, dim> &max_range) const
        {
            return mapBetweenDomains(point, a_min, b_max, min_range, max_range);
        };
        // return a New Point from ChaosFunction
        inline Point<T, dim> getPoint(const Point<T, dim> &point, const Point<T, dim> &min_global, const Point<T, dim> &max_global, int iter) const
        {
            Point<T, dim> new_point = toLocalDomain(point, min_global, max_global);
            new_point = generatePoint(new_point, iter);
            new_point = toGlobalDomain(new_point, min_global, max_global);

            return new_point;
        };
    };
};

#endif