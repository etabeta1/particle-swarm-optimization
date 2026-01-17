#ifndef _CONSTRAINT_HPP_
#define _CONSTRAINT_HPP_

#include <functional>
#include "point.hpp"

namespace Swarm
{
    /**
     * \brief Type definition for a generic constraint.
     * \tparam T The The data type for the coordinates (default is float).
     * \tparam dim The dimensionality of the points (default is 2)
     *
     * A constraint is nothing more than an std::function that takes a Point and returns a bool.
     * Constraints must return true if they are respected and false otherwise.
     */
    template <typename T = float, int dim = 2>
    using Constraint = std::function<bool(const Point<T, dim> &)>;

    /**
     * \brief Creates a box constraint defined by two points.
     * \tparam T The The data type for the coordinates (default is float).
     * \tparam dim The dimensionality of the points (default is 2)
     * \param a The minimum point defining the box.
     * \param b The maximum point defining the box.
     * \return A Constraint that checks if a point is inside the defined box.
     * \see Constraint
     *
     * The box constraint checks if all coordinates of a point are within the ranges defined by points `a` and `b`.
     */
    template <typename T = float, int dim = 2>
    Constraint<T, dim> makeBoxConstraint(const Point<T, dim> &a, const Point<T, dim> &b)
    {
        return [&a, &b](const Point<T, dim> &p)
        {
            return p.isInsideBox(a, b);
        }
    }

    /**
     * \brief Negates a given constraint.
     * \tparam T The The data type for the coordinates (default is float).
     * \tparam dim The dimensionality of the points (default is 2)
     * \param constraint The constraint to negate.
     * \return A Constraint that represents the negation of the given constraint.
     * \see Constraint
     *
     * The negated constraint returns true when the original constraint returns false, and vice versa.
     */
    template <typename T = float, int dim = 2>
    Constraint<T, dim> negateConstraint(const Constraint<T, dim> &constraint)
    {
        return [&constraint](const Point<T, dim> &p)
        {
            return !(constraint(p));
        }
    }

    /**
     * \brief Creates a maximum distance constraint from a center point.
     * \tparam T The The data type for the coordinates (default is float).
     * \tparam dim The dimensionality of the points (default is 2)
     * \param center The center point.
     * \param distance The maximum allowed distance from the center.
     * \return A Constraint that checks if a point is within the specified distance from the center.
     * \see Constraint
     *
     * The maximum distance constraint checks if the Euclidean distance between a point and the center is less than or equal to the specified distance.
     */
    template <typename T = float, int dim = 2>
    Constraint<T, dim> makeMaxDistanceConstraint(const Point<T, dim> &center, T distance)
    {
        return [&center, distance](const Point<T, dim> &p)
        {
            return (center - p).norm2() <= distance;
        }
    }

    /**
     * \brief Creates a elliptic constraint from 2 foci point and a tollerance.
     * \tparam T The The data type for the coordinates (default is float).
     * \tparam dim The dimensionality of the points (default is 2)
     * \param focus1  The first focus of the ellipse.
     * \param focus2  The second focus of the ellipse.
     * \param constant The length of the maximus axes of the ellipse.
     * \return A Constraint that checks if a point is within the ellipse.
     * \see Constraint
     *
     * The elliptic constraint checks if the a point is inside of the ellipse.
     */
    template <typename T = float, int dim = 2>
    Constraint<T, dim> makeEllipticConstraint(const Point<T, dim> &focus1, const Point<T, dim> &focus2, T constant)
    {
        return [&focus1, &focus2, constant](const Point<T, dim> &p)
        {
            return (focus1 - p).norm2() + (focus2 - p).norm2() <= constant;
        }
    }

}

#endif