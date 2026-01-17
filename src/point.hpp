#ifndef _POINT_HPP_
#define _POINT_HPP_

#include <vector>
#include <functional>
#include <ostream>
#include <cmath>

#include "utils.hpp"

namespace Swarm
{
    /**
     * \brief Class representing an n-dimensional vector.
     * \tparam T The data type for the coordinates (default is float).
     * \tparam dim The dimensionality of the point (default is 2).
     * \see Addable
     * \see Multipliable
     */
    template <typename T = float, int dim = 2>
        requires Addable<T, T> ||
                 Multipliable<T, T>
    class Point
    {
    public:
        /**
         * \brief Builds a `Point` using the elements from a vector.
         * \param data The vector containing the coordinate values.
         * \see std::vector
         */
        Point(std::vector<T> &data)
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = data[i];
            }
        }

        /**
         * \brief Builds a `Point` using the elements from a raw array.
         * \param data The array containing the coordinate values.
         */
        Point(T data[dim])
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = data[i];
            }
        }

        /**
         * \brief Builds a `Point` with all coordinates set to the same value.
         * \param data The value to set for all coordinates.
         */
        Point(T data)
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = data;
            }
        }

        /**
         * \brief Builds a `Point` as a copy of another `Point`.
         * \param data The `Point` to copy.
         */
        Point(const Point<T, dim> &data)
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = data.coordinates[i];
            }
        }

        /**
         * \brief Builds a `Point` using a generator function.
         * \param generator A function that generates the coordinate values based on their index.
         * \see std::function
         *
         * The lambda is called `dim` times, with an incremental index as paramenter that starts at 0 and then it is incremented by one after each iteration.
         */
        Point(const std::function<T(size_t)> &generator)
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = generator(i);
            }
        }

        /**
         * \brief Returns a new `Point` whose coordinates are obtained through the elementwise addition of the coordinates of two `Point`s.
         * \param other The other `Point` to add.
         * \return A new `Point` representing the elementwise sum.
         * \see Addable
         */
        template <typename U>
            requires Addable<T, U>
        Point<std::common_type_t<T, U>, dim> operator+(const Point<U, dim> &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] + other.coordinates[index]; });
        }

        /**
         * \brief Returns a new `Point` whose coordinates are obtained through the elementwise subtraction of the coordinates of two `Point`s.
         * \param other The other `Point` to subtract.
         * \return A new `Point` representing the elementwise difference.
         * \see Subtractable
         */
        template <typename U>
            requires Subtractable<T, U>
        Point<std::common_type_t<T, U>, dim> operator-(const Point<U, dim> &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] - other.coordinates[index]; });
        }

        /**
         * \brief Returns a new `Point` whose coordinates are obtained through the elementwise multiplication of the coordinates by a scalar.
         * \param other The scalar to multiply with.
         * \return A new `Point` representing the elementwise product.
         * \see Multipliable
         */
        template <typename U>
            requires Multipliable<T, U>
        Point<std::common_type_t<T, U>, dim> operator*(const U &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] * other; });
        }

        /**
         * \brief Returns a new `Point` whose coordinates are obtained through the elementwise multiplication of the coordinates of two `Point`s.
         * \param other The other `Point` to multiply with.
         * \return A new `Point` representing the elementwise product.
         * \see Multipliable
         */
        template <typename U>
            requires Multipliable<T, U>
        Point<std::common_type_t<T, U>, dim> operator*(const Point<U, dim> &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] * other.coordinates[index]; });
        }

        /**
         * \brief Returns a new `Point` whose coordinates are obtained through the elementwise division of the coordinates by a scalar.
         * \param other The scalar to divide by.
         * \return A new `Point` representing the elementwise quotient.
         * \see Divisible
         */
        template <typename U>
            requires Divisible<T, U>
        Point<std::common_type_t<T, U>, dim> operator/(const U &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] / other; });
        }

        /**
         * \brief Returns a new `Point` whose coordinates are obtained through the elementwise division of the coordinates of two `Point`s.
         * \param other The other `Point` to divide by.
         * \return A new `Point` representing the elementwise quotient.
         * \see Divisible
         */
        template <typename U>
            requires Divisible<T, U>
        Point<std::common_type_t<T, U>, dim> operator/(const Point<U, dim> &other) const
        {
            return Point<T, dim>([this, &other](size_t index)
                                 { return coordinates[index] / other.coordinates[index]; });
        }

        /**
         * \brief Returns a new `Point` whose coordinates are clamped between the corresponding coordinates of two other `Point`s.
         * \param a The minimum `Point`.
         * \param b The maximum `Point`.
         * \return A new `Point` with clamped coordinates.
         * \see TriComparable
         */
        template <typename U, typename V>
            requires TriComparable<U, T, V>
        Point<T, dim> clamp(const Point<U, dim> &a, const Point<V, dim> &b) const
        {
            return Point<T, dim>([this, &a, &b](size_t index)
                                 { return std::max(a[index], std::min(coordinates[index], b[index])); });
        }

        /**
         * \brief Returns a new `Point` whose coordinates are obtained through the elementwise application of `std::sin`.
         * \return A new `Point` with the sine of each coordinate.
         * \see std::sin
         */
        Point<T, dim> sin() const
        {
            return Point<T, dim>([this](size_t index)
                                 { return std::sin(coordinates[index]); });
        }

        /**
         * \brief Returns a new `Point` whose coordinates are obtained through the elementwise application of `std::asin`.
         * \return A new `Point` with the arcsine of each coordinate.
         * \see std::asin
         */
        Point<T, dim> arcsin() const
        {
            return Point<T, dim>([this](size_t index)
                                 { return std::asin(coordinates[index]); });
        }

        /**
         * \brief Returns a new `Point` whose coordinates are obtained through the elementwise application of `std::cos`.
         * \return A new `Point` with the cosine of each coordinate.
         * \see std::cos
         */
        Point<T, dim> cos() const
        {
            return Point<T, dim>([this](size_t index)
                                 { return std::cos(coordinates[index]); });
        }

        /**
         * \brief Returns a new `Point` whose coordinates are obtained through the elementwise application of `std::acos`.
         * \return A new `Point` with the arccosine of each coordinate.
         * \see std::acos
         */
        Point<T, dim> arccos() const
        {
            return Point<T, dim>([this](size_t index)
                                 { return std::acos(coordinates[index]); });
        }

        /**
         * \brief Returns a new `Point` whose coordinates are obtained through the elementwise exponentiation to a positive integer.
         * \param exponent The exponent to raise each coordinate to.
         * \return A new `Point` with the coordinates raised to the power of `exponent`.
         */
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

        /**
         * \brief Returns the 1-norm of the point.
         * \return The 1-norm value.
         */
        T norm1() const
        {
            T sum = 0;
            for (T c : coordinates)
            {
                sum += c;
            }
            return sum;
        }

        /**
         * \brief Returns the squared 2-norm of the point.
         * \return The squared 2-norm value.
         */
        T squareNorm2() const
        {
            T sum = 0;
            for (T c : coordinates)
            {
                sum += c * c;
            }
            return sum;
        }

        /**
         * \brief Returns the 2-norm of the point.
         * \return The 2-norm value.
         */
        T norm2() const
        {
            return std::sqrt(squareNorm2());
        }

        /**
         * \brief Returns a new `Point` whose coordinates are the absolute values of the original coordinates.
         * \return A new `Point` with the absolute values of each coordinate.
         * \see std::abs
         */
        Point<T, dim> abs() const
        {
            return Point<T, dim>([this](size_t index)
                                 { return std::abs(coordinates[index]); });
        }

        /**
         * \brief Returns the coordinate at the selected index.
         * \param index The index of the coordinate to retrieve.
         * \return The value of the coordinate at the specified index.
         */
        T operator[](size_t index) const
        {
            return coordinates[index];
        }

        /**
         * \brief Returns a reference to the coordinate at the selected index.
         * \param index The index of the coordinate to retrieve.
         * \return A reference to the value of the coordinate at the specified index.
         */
        T &operator[](size_t index)
        {
            return coordinates[index];
        }

        /**
         * \brief Returns the dimensionality of the point.
         * \return The number of dimensions.
         */
        constexpr int getDimensions() const
        {
            return dim;
        }

        /**
         * \brief Makes `Point` printable with the standard `<<` operator.
         */
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

        /**
         * \brief Assignment operator for `Point`.
         * \param other The other `Point` to assign from.
         * \return A reference to the assigned `Point`.
         */
        Point<T, dim> &operator=(Point<T, dim> other)
        {
            for (size_t i = 0; i < dim; i++)
            {
                coordinates[i] = other[i];
            }

            return *this;
        }

        /**
         * \brief Checks if the point is inside the box defined by two other points.
         * \param a One corner of the box.
         * \param b The opposite corner of the box.
         * \return `true` if the point is inside the box, `false` otherwise
         */
        bool isInsideBox(const Point<T, dim> &a, const Point<T, dim> &b) const
        {
            for (size_t i = 0; i < dim; i++)
            {
                if (coordinates[i] < std::min(a[i], b[i]) || coordinates[i] > std::max(a[i], b[i]))
                {
                    return false;
                }
            }
            return true;
        }

    private:
        /**
         * \brief The array storing the coordinates of the point.
         */
        T coordinates[dim];
    };
};

#endif