#ifndef _UTILS_HPP_
#define _UTILS_HPP_

namespace Swarm
{
    using IterationType = uint64_t;

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
};

#endif