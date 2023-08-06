#include "MathematicaConversions.h"



#include <boost/math/special_functions/airy.hpp> // airy functions


template<typename T, typename U>
double Power(T a, U b)
{
    return pow(a, b);
}

template<typename T>
T AiryAi(T x)
{
    return boost::math::airy_ai(x);
}

template<typename T>
T AiryBi(T x)
{
    return boost::math::airy_bi(x);
}

template<typename T>
T AiryAiPrime(T x)
{
    return boost::math::airy_ai_prime(x);
}

template<typename T>
T AiryBiPrime(T x)
{
    return boost::math::airy_bi_prime(x);
}
