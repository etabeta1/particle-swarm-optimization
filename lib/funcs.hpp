#ifndef _FUNCS_HPP_
#define _FUNCS_HPP_

#include <random>
#include <chrono>

inline float generate_random(float a, float b)
{
    static std::mt19937 generator(std::random_device{}());
    std::uniform_real_distribution<float> distribution(a, b);
    return distribution(generator);
}

#endif