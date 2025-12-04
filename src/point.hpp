#ifndef _POINT_HPP_
#define _POINT_HPP_

#include <vector>

namespace Swarm
{
    template<typename T, typename U>
    concept Addable = requires(T a, U b) {
        a + b;
    };

    template<typename T, typename U>
    concept Subtractable = requires(T a, U b) {
        a - b;
    };

    template<typename T, typename U>
    concept Multipliable = requires(T a, U b) {
        a * b;
    };

    template<typename T, typename U, typename V>
    concept TriComparable = requires(T a, U b, V c) {
        a <= b <= c;
    };

    template <typename T = float, int dim = 2>
    class Point
    {
    public:
        Point(std::vector<T>& data) {
            for(size_t i = 0; i < dim; i++) {
                coordinates[i] = data[i];
            }
        }

        Point(T data[dim]) {
            for(size_t i = 0; i < dim; i++) {
                coordinates[i] = data[i];
            }
        }

        Point(T data) {
            for(size_t i = 0; i < dim; i++) {
                coordinates[i] = data;
            }
        }

        Point(Point<T, dim>& data) {
            for(size_t i = 0; i < dim; i++) {
                coordinates[i] = data.coordinates[i];
            }
        }

        Point(std::function<T(size_t)> generator) {
            for(size_t i = 0; i < dim; i++) {
                coordinates[i] = generator(i);
            }
        }

        template<typename U> requires Addable<T, U>
        Point<T, dim> operator+(Point<U, dim>& other) const {}

        template<typename U> requires Subtractable<T, U>
        Point<T, dim> operator-(Point<U, dim>& other) const {}

        template<typename U> requires Multipliable<T, U>
        Point<T, dim> operator*(U other) const {}

        template<typename U> requires Multipliable<T, U>
        Point<T, dim> operator*(Point<U, dim>& other) const {}

        T& operator[](size_t index) const {}

        template<typename U, typename V> requires TriComparable<U, T, V>
        Point<T, dim> clamp(Point<U, dim>& a, Point<V, dim>& b) {}
    
    private:
        T coordinates[dim];
    };
};

#endif