#ifndef _FUNCS_HPP_
#define _FUNCS_HPP_

#include <random>
#include <chrono>

template<typename T = float>
inline T generate_random(T a, T b)
{
    static std::mt19937 generator(std::random_device{}());
    
    if constexpr (std::is_integral_v<T>) {
        std::uniform_int_distribution<T> distribution(a, b);
        return distribution(generator);
    } else {
        std::uniform_real_distribution<T> distribution(a, b);
        return distribution(generator);
    }
}

#endif