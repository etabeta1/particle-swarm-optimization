#ifndef _POINT_HPP_
#define _POINT_HPP_

#include <vector>
#include <functional>
#include <ostream>
#include <cmath>

#include "utils.hpp"

namespace Swarm
{
    // Representation of an n-dimensional vector.
    //
    // - `T` is the type used to store the coordinated (defaults to `float`).
    // - `dim` is the number of dimensions for the vector (defaults to 2).
    template <typename T = float, int dim = 2>
        requires Addable<T, T> ||
                 Multipliable<T, T>
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

        // Returns a new `Point` whose coordinates are obtained through the elementwise application of `std::sin`.
        Point<T, dim> sin() const
        {
            return Point<T, dim>([this](size_t index)
                                 { return std::sin(coordinates[index]); });
        }

        // Returns a new `Point` whose coordinates are obtained through the elementwise application of `std::asin`.
        Point<T, dim> arcsin() const
        {
            return Point<T, dim>([this](size_t index)
                                 { return std::asin(coordinates[index]); });
        }

        // Returns a new `Point` whose coordinates are obtained through the elementwise application of `std::cos`.
        Point<T, dim> cos() const
        {
            return Point<T, dim>([this](size_t index)
                                 { return std::cos(coordinates[index]); });
        }

        // Returns a new `Point` whose coordinates are obtained through the elementwise application of `std::acos`.
        Point<T, dim> arccos() const
        {
            return Point<T, dim>([this](size_t index)
                                 { return std::acos(coordinates[index]); });
        }

        // Returns a new `Point` whose coordinates are obtained through the elementwise exponentiation to a positive integer.
        Point<T, dim> pow(int exponent) const
        {
            Point<T, dim> remainder(1);
            Point<T, dim> accum(*this);

            while (exponent > 1)
            {
                if (exponent & 1)
                {
                    remainder = remainder * accum;
                }

                accum = accum * accum;
                exponent >>= 1;
            }

            return accum * remainder;
        }

        // Returns the 1-norm of the point.
        T norm1() const
        {
            T sum = 0;
            for (T c : coordinates)
            {
                sum += c;
            }
            return sum;
        }

        // Returns the square 2-norm of the point.
        T squareNorm2() const
        {
            T sum = 0;
            for (T c : coordinates)
            {
                sum += c * c;
            }
            return sum;
        }

        // Returns the 2-norm of the point.
        T norm2() const
        {
            return std::sqrt(squareNorm2());
        }

        // Returns a new `Point` whose coordinates are obtained through the elementwise application of `std::abs`.
        Point<T, dim> abs() const
        {
            return Point<T, dim>([this](size_t index)
                                 { return std::abs(coordinates[index]); });
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