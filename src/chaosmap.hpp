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
        std::function<U(Point<T, dim> &)> map;
        inline Point<T, dim> mapBetweenDomains(const Point<T, dim> &point,const Point<T, dim> & min_range,const Point<T, dim> & max_range){
            //codice per trasformare le coordinate
            Point<T,dim> new_point;
            for(int i=0;i<dim;i++){
                new_point[i] = 
            }
        };
        inline Point<T,dim> generatePoint(const Point<T, dim> &point){
            //map evaluate point map(Point),per ogni dimensione genera per ogni coordinata
            Point<T,dim> new_point;
            for(int i==0;i<dim;i++){
                new_point[i]=map(point[i]);
            }
            return new_point;
        };

    public:
        ChaosMap(std::function<U(Point<T, dim> &)> lambda_fun,Point<T,dim> min,Point<T,dim> max){
            map = lambda_fun;
            a_min = min;
            b_max = max;
        };
        inline Point<T, dim> toLocalDomain(const Point<T, dim> &point){
            return mapBetweenDomains(point,a_min,b_max);
        };
        inline Point<T, dim> toGlobalDomain(const Point<T, dim> &point,const Point<T, dim> & min_range,const Point<T, dim> & max_range){
            return mapBetweenDomains(point,min_range,max_range);
        };
        inline Point<T,dim> getPoint(const Point<T, dim> &point,const Point<T, dim> & min_global,const Point<T, dim> & max_global){
            Point<T,dim> new_point;

            new_point = toLocalDomain(point);
            new_point = generatePoint(new_point);
            new_point = toGlobalDomain(new_point,min_global,max_global);

            return new_point;

        };
    };
};

#endif