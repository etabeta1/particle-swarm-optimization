#ifndef _POINT_HPP_
#define _POINT_HPP_

#include <vector>
#include <functional>
#include <ostream>

namespace Swarm
{
    template <typename T, typename U>
    concept Addable = requires(T a, U b) {
        a + b;
    };

    template <typename T, typename U>
    concept Subtractable = requires(T a, U b) {
        a - b;
    };

    template <typename T, typename U>
    concept Multipliable = requires(T a, U b) {
        a * b;
    };

    template <typename T, typename U>
    concept Divisible = requires(T a, U b) {
        a / b;
    };

    template <typename T, typename U, typename V>
    concept TriComparable = requires(T a, U b, V c) {
        a <= b <= c;
    };

    template <typename T = float, int dim = 2>
    class Point
    {
    public:
        Point(std::vector<T> &data)
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = data[i];
            }
        }

        Point(T data[dim])
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = data[i];
            }
        }

        Point(T data)
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = data;
            }
        }

        Point(const Point<T, dim> &data)
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = data.coordinates[i];
            }
        }

        Point(const std::function<T(size_t)> &generator)
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = generator(i);
            }
        }

        template <typename U>
            requires Addable<T, U>
        Point<T, dim> operator+(const Point<U, dim> &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] + other.coordinates[index]; });
        }

        template <typename U>
            requires Subtractable<T, U>
        Point<T, dim> operator-(const Point<U, dim> &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] - other.coordinates[index]; });
        }

        template <typename U>
            requires Multipliable<T, U>
        Point<T, dim> operator*(const U &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] * other; });
        }

        template <typename U>
            requires Multipliable<T, U>
        Point<T, dim> operator*(Point<U, dim> &other) const
        {
        }

        template <typename U>
            requires Divisible<T, U>
        Point<T, dim> operator/(Point<U, dim> &other) const
        {
        }

        template <typename U>
            requires Divisible<T, U>
        Point<T, dim> operator/(U other) const
        {
        }

        T &operator[](size_t index) const {}

        constexpr int getDimensions() const
        {
            return dim;
        }

        template <typename U, typename V>
            requires TriComparable<U, T, V>
        Point<T, dim> clamp(Point<U, dim> &a, Point<V, dim> &b)
        {
        }

        friend std::ostream &operator<<(std::ostream &os, const Point<T, dim> &p)
        {
            os << "[ ";
            for (size_t i = 0; i < dim; i++)
            {
                os << p.coordinates[i] << ", ";
            }
            os << "]";

            return os;
        }

    private:
        T coordinates[dim];
    };
};

#endif