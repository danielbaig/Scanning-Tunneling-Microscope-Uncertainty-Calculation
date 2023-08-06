#ifndef RANDOM_H
#define RANDOM_H

class Random
{
public:
    static double uniform();

    static double randn();

    static int intUniform(int lower, int upper);
};

#endif