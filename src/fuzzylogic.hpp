#ifndef _FUZZYLOGIC_HPP_
#define _FUZZYLOGIC_HPP_

#include <algorithm>
#include <cmath>

namespace Swarm
{
    namespace FuzzyLogic
    {
        /**
         * \brief Membership functions for Delta (distance from global best)
         * \tparam T The type used for calculations (defaults to float)
         *
         * Implements trapezoidal and triangular membership functions for the
         * linguistic variable "Distance from global best" with three fuzzy sets:
         * - Same: particle is very close to global best
         * - Near: particle is at medium distance
         * - Far: particle is far from global best
         */
        template <typename T = float>
        struct DeltaMembership
        {
            /**
             * \brief "Same" is a trapezoidal membership function
             * \param delta Current distance value
             * \param delta_max Maximum possible distance for the linguistic value "Same"
             * \return Membership value in [0,1]
             */
            static T membershipSame(T delta, T delta_max)
            {
                T d1 = static_cast<T>(0.2) * delta_max;
                T d2 = static_cast<T>(0.4) * delta_max;

                if (delta < d1)
                    return static_cast<T>(1.0);
                if (delta >= d2)
                    return static_cast<T>(0.0);
                return (d2 - delta) / (d2 - d1);
            }

            /**
             * \brief "Near" is a triangular membership function
             * \param delta Current distance value
             * \param delta_max Maximum possible distance for the linguistic value "Near"
             * \return Membership value in [0,1]
             */
            static T membershipNear(T delta, T delta_max)
            {
                T d1 = static_cast<T>(0.2) * delta_max;
                T d2 = static_cast<T>(0.4) * delta_max;
                T d3 = static_cast<T>(0.6) * delta_max;

                if (delta < d1 || delta >= d3)
                    return static_cast<T>(0.0);
                if (delta < d2)
                    return (delta - d1) / (d2 - d1);
                return (d3 - delta) / (d3 - d2);
            }

            /**
             * \brief "Far" is a trapezoidal membership function
             * \param delta Current distance value
             * \param delta_max Maximum possible distance for the linguistic value "Far"
             * \return Membership value in [0,1]
             */
            static T membershipFar(T delta, T delta_max)
            {
                T d2 = static_cast<T>(0.4) * delta_max;
                T d3 = static_cast<T>(0.6) * delta_max;

                if (delta >= d3)
                    return static_cast<T>(1.0);
                if (delta < d2)
                    return static_cast<T>(0.0);
                return (delta - d2) / (d3 - d2);
            }
        };

        /**
         * \brief Membership functions for Phi (normalized fitness incremental factor)
         * \tparam T The type used for calculations (defaults to float)
         *
         * Implements triangular membership functions for the
         * linguistic variable "Normalized fitness incremental factor" with three fuzzy sets:
         * - Better: particle is very close to global best
         * - Same: particle is at medium distance
         * - Worse: particle is far from global best
         */
        template <typename T = float>
        struct PhiMembership
        {
            /**
             * \brief Triangular membership function for "Better"
             * \param phi Normalized fitness incremental factor in [-1, 1]
             * \return Membership degree in [0, 1]
             */
            static T membershipBetter(T phi)
            {
                if (phi <= static_cast<T>(-1.0))
                    return static_cast<T>(1.0);
                if (phi >= static_cast<T>(0.0))
                    return static_cast<T>(0.0);
                return -phi;
            }

            /**
             * \brief Triangular membership function for "Same"
             * \param phi Normalized fitness incremental factor in [-1, 1]
             * \return Membership degree in [0, 1]
             */
            static T membershipSame(T phi)
            {
                return static_cast<T>(1.0) - std::abs(phi);
            }

            /**
             * \brief Triangular membership function for "Worse"
             * \param phi Normalized fitness incremental factor in [-1, 1]
             * \return Membership degree in [0, 1]
             */
            static T membershipWorse(T phi)
            {
                if (phi <= static_cast<T>(0.0))
                    return static_cast<T>(0.0);
                if (phi >= static_cast<T>(1.0))
                    return static_cast<T>(1.0);
                return phi;
            }
        };

        /**
         * \brief Sugeno defuzzification method
         * \tparam T The type used for calculations (defaults to `float`)
         * \param low_deg Membership degree for "Low" linguistic value
         * \param low_val Crisp value for "Low" (singleton)
         * \param med_deg Membership degree for "Medium" linguistic value
         * \param med_val Crisp value for "Medium" (singleton)
         * \param high_deg Membership degree for "High" linguistic value
         * \param high_val Crisp value for "High" (singleton)
         * \return Defuzzified crisp output value (weighted average)
         *
         * Computes the weighted average of the consequents based on their
         * activation degrees. This is the standard Sugeno inference method.
         */
        template <typename T = float>
        T sugenoDefuzzify(T low_deg, T low_val,
                          T med_deg, T med_val,
                          T high_deg, T high_val)
        {
            T sum_degrees = low_deg + med_deg + high_deg;

            // Avoid division by zero - return medium value as fallback
            if (sum_degrees < static_cast<T>(1e-10))
                return med_val;

            return (low_deg * low_val + med_deg * med_val + high_deg * high_val) / sum_degrees;
        }

    }
}
#endif