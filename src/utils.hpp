#ifndef _UTILS_HPP_
#define _UTILS_HPP_

#include <cstdint>

namespace Swarm
{
    /**
     * \brief Type alias for iteration count.
     */
    using IterationType = uint64_t;

    /**
     * \brief Concept requiring that the selected types supports addition.
     * \tparam T The first type.
     * \tparam U The second type.
     */
    template <typename T, typename U>
    concept Addable = requires(T a, U b) {
        a + b;
    };

    /**
     * \brief Concept requiring that the selected types supports subtraction.
     * \tparam T The first type.
     * \tparam U The second type.
     */
    template <typename T, typename U>
    concept Subtractable = requires(T a, U b) {
        a - b;
    };

    /**
     * \brief Concept requiring that the selected types supports multiplication.
     * \tparam T The first type.
     * \tparam U The second type.
     */
    template <typename T, typename U>
    concept Multipliable = requires(T a, U b) {
        a * b;
    };

    /**
     * \brief Concept requiring that the selected types supports division.
     * \tparam T The first type.
     * \tparam U The second type.
     */
    template <typename T, typename U>
    concept Divisible = requires(T a, U b) {
        a / b;
    };

    /**
     * \brief Concept requiring that the selected types supports comparison.
     * \tparam T The first type.
     * \tparam U The second type.
     */
    template <typename T, typename U>
    concept Comparable = requires(T a, U b) {
        a <= b;
    };

    /**
     * \brief Concept requiring that the selected types supports triple comparison.
     * \tparam T The first type.
     * \tparam U The second type.
     * \tparam V The third type.
     */
    template <typename T, typename U, typename V>
    concept TriComparable = requires(T a, U b, V c) {
        a <= b <= c;
    };
};

#endif