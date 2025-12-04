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

    // Representation of an n-dimensional vector.
    //
    // - `T` is the type used to store the coordinated (defaults to `float`).
    // - `dim` is the number of dimensions for the vector (defaults to 2).
    template <typename T = float, int dim = 2>
    class Point
    {
    public:
        // Builds a `Point` using the first `dim` elements of an `std::vector`.
        Point(std::vector<T> &data)
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = data[i];
            }
        }

        // Builds a `Point` using the elements from an array.
        Point(T data[dim])
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = data[i];
            }
        }

        // Builds a `Point` where all the members are set to the same constant.
        Point(T data)
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = data;
            }
        }

        // Copy constructor
        Point(const Point<T, dim> &data)
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = data.coordinates[i];
            }
        }

        // Builds a `Point` using the returned values from a given lambda.
        //
        // The lambda is called `dim` times, with an incremental index as paramenter that starts at 0 and then it is incremented by one after each iteration.
        Point(const std::function<T(size_t)> &generator)
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = generator(i);
            }
        }

        // Returns a new `Point` whose coordinates are obtained through the elementwise sum of the coordinates of two `Point`s.
        template <typename U>
            requires Addable<T, U>
        Point<T, dim> operator+(const Point<U, dim> &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] + other.coordinates[index]; });
        }

        // Returns a new `Point` whose coordinates are obtained through the elementwise subtraction of the coordinates of two `Point`s.
        template <typename U>
            requires Subtractable<T, U>
        Point<T, dim> operator-(const Point<U, dim> &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] - other.coordinates[index]; });
        }

        // Returns a new `Point` whose coordinates are obtained through the elementwise product of the coordinates by a scalar.
        template <typename U>
            requires Multipliable<T, U>
        Point<T, dim> operator*(const U &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] * other; });
        }

        // Returns a new `Point` whose coordinates are obtained through the elementwise product of the coordinates of two `Point`s.
        template <typename U>
            requires Multipliable<T, U>
        Point<T, dim> operator*(const Point<U, dim> &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] * other.coordinates[index]; });
        }

        // Returns a new `Point` whose coordinates are obtained through the elementwise division of the coordinates by a scalar.
        template <typename U>
            requires Divisible<T, U>
        Point<T, dim> operator/(const U &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] / other; });
        }

        // Returns a new `Point` whose coordinates are obtained through the elementwise division of the coordinates of two `Point`s.
        template <typename U>
            requires Divisible<T, U>
        Point<T, dim> operator/(const Point<U, dim> &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] / other.coordinates[index]; });
        }

        // Returns a new `Point` whose coordinates are obtained through the elementwise clamping between the coordinates of two other `Point`s defining an n-dimensional box.
        template <typename U, typename V>
            requires TriComparable<U, T, V>
        Point<T, dim> clamp(const Point<U, dim> &a, const Point<V, dim> &b) const
        {
            return Point<T, dim>([this, &a, &b](size_t index)
                                 { return std::max(a[index], std::min(coordinates[index], b[index])); });
        }

        // Returns a copy of the coordinate at the selected index.
        T operator[](size_t index) const
        {
            return coordinates[index];
        }

        // Returns a reference to the coordinate at the selected index.
        T &operator[](size_t index)
        {
            return coordinates[index];
        }

        // Returns the number of dimensions.
        constexpr int getDimensions() const
        {
            return dim;
        }

        // Makes `Point` printable with the standard `<<` operator.
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
        // Stores the value of each coordinate.
        T coordinates[dim];
    };
};

#endif