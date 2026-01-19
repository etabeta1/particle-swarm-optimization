#ifndef _FUNCS_HPP_
#define _FUNCS_HPP_

#include <random>
#include <chrono>

/**
 * \brief Generates a random number within the specified range [a, b].
 * \tparam T The data type of the random number (can be integral or floating-point).
 * \param a The lower bound of the range.
 * \param b The upper bound of the range.
 * \return A random number of type T within the range [a, b].
 * \note This function uses the Mersenne Twister engine for random number generation.
 */
template <typename T = float>
inline T generate_random(T a, T b)
{
    static std::mt19937 generator(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    if constexpr (std::is_integral_v<T>)
    {
        std::uniform_int_distribution<T> distribution(a, b);
        return distribution(generator);
    }
    else
    {
        std::uniform_real_distribution<T> distribution(a, b);
        return distribution(generator);
    }
}

#endif