#ifndef MATHEMATICAROOTSOLVE_H
#define MATHEMATICAROOTSOLVE_H



double mathematicaCalculatedFunction(const double w, const double r,
    const double mass,
    const double energy, const double bias,
    const double workFunction, const double k1, const double k3, const double gamma);



double newton_raphson(double (*f)(const double, const double, const double,
    const double, const double, const double, const double, const double, const double),
    double guess,
    const unsigned int numIterations, const double r,
    const double mass,
    const double energy, const double bias,
    const double workFunction, const double k1, const double k3, const double gamma = 100.);







#endif
