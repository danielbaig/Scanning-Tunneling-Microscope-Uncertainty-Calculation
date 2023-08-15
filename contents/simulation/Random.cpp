#include "Random.h"

#include<iostream>
#include<random>



double Random::uniform()
{
    /*
    Generate a random number between 0 and 1.

    Outputs:
        - double m_uniform: Uniform random number.
    */

    static std::mt19937 mt{ static_cast<std::mt19937::result_type>
    (std::time(nullptr)) };
    static std::uniform_real_distribution<double> m_uniform{ 0,1 };
    return m_uniform(mt);
}

double Random::randn()
{
    /*
    Generate a random number based on a standard normal
    distribution.

    Outputs:
        - double m_randn: Normal random number.
    */

    static std::mt19937 mt{ static_cast<std::mt19937::result_type>
    (std::time(nullptr)) };
    static std::normal_distribution<double> m_randn{ 0, 1 };
    return m_randn(mt);
}

int Random::intUniform(int lower, int upper)
{
    /*
    Generates a random integer between given bounds.

    Inputs:
        - int lower: Lower bound.
        - int upper: Upper bound.

    Outputs:
        - int value: Random number.
    */

    static std::mt19937 mt{ static_cast<std::mt19937::result_type>
    (std::time(nullptr)) };

    if (lower >= upper)
    {
        std::cout << "lower >= upper in intUniform\n";
        std::cout << "Lower: " << lower << "\nUpper: " 
            << upper << '\n';
        return -1;
    }
    std::uniform_int_distribution<int> value{ lower, upper-1 };
    return value(mt);
}
