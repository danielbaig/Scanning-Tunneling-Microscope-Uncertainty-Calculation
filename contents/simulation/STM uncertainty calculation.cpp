#include "Random.h"
#include "Timer.h"
#include "MathematicaRootSolve.h"


#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <fstream> // Saving to file
#include <algorithm> // std::remove, std::fill_n
#include <numeric> // std::iota
#include <tuple> // std::tuple, std::make_tuple
#include <complex> // std::complex

#include <math.h> // exp
#include <stdlib.h> // abs
#include <boost/math/special_functions/airy.hpp> // airy functions
#include <boost/math/constants/constants.hpp> // constants



constexpr double ELECTRONCHARGE{ 1.60217663e-19 }; // [C]
constexpr double ELECTRONMASS{ 9.1093837e-31 }; // [kg]
constexpr double HBAR{ 1.054571817e-34 }; // [Js]
const double PI{ boost::math::constants::pi<double>() };
constexpr std::complex<double> I(0, 1.);





std::vector<double> linspace(double start, double end, int numPoints)
{
    std::vector<double> arr(numPoints);

    const double dx{ (end - start) / numPoints };
    arr[0] = start;

    for (int i{ 1 }; i < numPoints; ++i)
    {
        arr[i] = arr[i - 1] + dx;
    }

    if (abs(arr[numPoints] - end) < 1e-3)
        return arr;
    else
    {
        std::cout << "Error occured\n";
        return arr;
    }


}

template<typename T>
T pythonModulo(T x, T M)
{
    T output{ fmod(x, M) };

    if (M > 0 && output < 0)
    {
        return M + output;
    }
    else if (M < 0 && output > 0)
    {
        return output - M;
    }
    
    return output;
}


class STM
{
    /*
    Creates a system to model electrons tunneling in a Scanning Tunneling
    Microscope by using the wavefunction to generate a PDF and then using
    this to determine the transmittance to estimate the height of the prob.
    
    Member variables:

    */

private:
    const unsigned int m_numReadings{};


    struct SamplePoint
    {
        double x{};
        double y{};

        double* tunnelProb{};
        double* estimatedTransmittance{};
    };

    const double m_GAMMA{ 538 };
    const unsigned int m_numSamples{};
    SamplePoint* m_samples{ new SamplePoint[m_numSamples] {} };

    const unsigned int m_maxNumElectrons{};

    // System properties
    double m_workFunction{5}; // [eV]
    double m_bias{0.3}; // V
    double m_energy{2.5}; // [eV]

    const double m_k1{};
    const double m_k3{};


    double m_actualWidth{};
    std::complex<double> m_fluxRef{};
    std::complex<double> m_fluxOut{};
    double m_ai0{};
    double m_aih{};
    double m_bi0{};
    double m_bih{};
    double m_M_det{};



public:
    STM(const unsigned int numSamples, double lowerBound, double upperBound,
        const unsigned int maxNumElectrons, const unsigned int numReadings) 
        : m_numSamples{numSamples},
        m_maxNumElectrons{maxNumElectrons},
        m_numReadings{numReadings},
        m_k1{ sqrt(2 * ELECTRONMASS * m_energy*ELECTRONCHARGE) / HBAR },
        m_k3{ sqrt(2 * ELECTRONMASS * (m_energy + m_bias)*ELECTRONCHARGE) / HBAR }
    {
        /*
        Constructor for the class.

        Inputs:
            - const unsigned int numSamples: The number of samples to take.
            - double lowerBound: The lower bound of the system.
            - double upperBound: The upper bound of the system.
            - const unsigned int maxNumElectrons: Maximum number of electrons
            to try to tunnel.
            - const unsigned int numReadings: The number of transmittance readings
            to take.
        */

        // Sets random positions to the samples.
        for (int i{ 0 }; i < m_numSamples; ++i)
        {
            m_samples[i].x = (upperBound - lowerBound)*Random::uniform() + lowerBound;
            m_samples[i].y = ((upperBound - lowerBound)*Random::uniform() + lowerBound)* sin(PI / 3);
            m_samples[i].tunnelProb = new double[m_numReadings] {0};
            m_samples[i].estimatedTransmittance = new double[m_numReadings] {0};
        }

        // Convert to J
        m_workFunction *= ELECTRONCHARGE;
        m_bias *= ELECTRONCHARGE;
        m_energy *= ELECTRONCHARGE;

    }


    ~STM()
    {
        /*
        Destructor for the class.
        Deletes dynamic arrays.
        */
        for (int i{ 0 }; i < m_numSamples; ++i)
        {
            delete[] m_samples[i].tunnelProb;
            m_samples[i].tunnelProb = nullptr;
            delete[] m_samples[i].estimatedTransmittance;
            m_samples[i].estimatedTransmittance = nullptr;
        }


        delete[] m_samples;
        m_samples = nullptr;
    }


    void sampleLoop()
    {
        /*
        Samples all of the points and saves them.
        
        */


        for (int i{ 0 }; i < m_numSamples; ++i)
        {
            
            if ((i + 1) % (m_numSamples / 5) == 0)
            {
                std::cout << "Progress: " << 100 * (i + 1) / m_numSamples << "%\n";
            }

            samplePoint(&(m_samples[i]));

                
            for (int j{ 0 }; j < m_numReadings; ++j)
            {
                
                m_samples[i].estimatedTransmittance[j]
                    = newton_raphson(
                        mathematicaCalculatedFunction, 4.5e-10, 100,
                        m_samples[i].tunnelProb[j], ELECTRONMASS, m_energy, m_bias, m_workFunction,
                        m_k1, m_k3, m_GAMMA);
            }

        }
        




        saveSamples();
        savePositions();

    }

    void savePositions() const
    {
        std::ofstream myfile;
        myfile.open("positions.txt");

        for (int i{ 0 }; i < m_numSamples; ++i)
        {
            myfile << m_samples[i].x << ',' << m_samples[i].y << '\n';
        }
        myfile.close();
    }

    void saveSamples() const
    {
        std::ofstream myfileSamples;
        myfileSamples.open("samples.txt");
        std::ofstream myfileProp;
        myfileProp.open("proportionTunnelled.txt");

        for (int i{ 0 }; i< m_numSamples; ++i)
        {
            for (int j{ 0 }; j < m_numReadings; ++j)
            {
                myfileSamples << m_samples[i].estimatedTransmittance[j] << ',';
                myfileProp << m_samples[i].tunnelProb[j] << ',';
            }
            myfileSamples << '\n';
            myfileProp << '\n';
        }
        myfileSamples.close();
        myfileProp.close();
    }


    double height(const double x, const double y) const
    {
        /*
        Determines the actual height of the sample.
        Width between (4-5)e-10m

        Inputs:
            - const double x: x-coord of the tip.
            - const double y: y-coord of the tip.
        Outputs:
            - double z: Actual height of the sample.
        */

        static constexpr double verticalOffset{ 0.5 };
        static constexpr double horizontalOffset{ 1 };
        static constexpr unsigned int gap{ 4 };
        static constexpr double particleRadius{ 1e-10 };
        static const double verticalSeparation{ sin(PI / 3) };
        static const double tipHeight{ 5 * particleRadius };

        // Distance between tip and particle centre in even rows.
        double r{ hypot(pythonModulo(y / particleRadius - verticalOffset, gap * verticalSeparation) - 1.,
            pythonModulo(x / particleRadius - horizontalOffset, 2.) - 1.) };
        double z { -particleRadius };

        if (r < 1)
        {
            z = particleRadius * sqrt(1 - r * r);
            return tipHeight - z;
        }

        // Distance between tip and particle centre in odd rows.
        r = hypot(pythonModulo(y / particleRadius - verticalOffset 
            + 2 * verticalSeparation, gap * verticalSeparation) - 1.,
            pythonModulo(x / particleRadius - horizontalOffset - 1., 2.) - 1.);

        if (r < 1)
            z = particleRadius * sqrt(1 - r * r);

        return tipHeight - z;


    }

    double prob_psi1(std::complex<double> x) const
    {
        /*
        Gets the probability density of a given point for the incident
        wave.

        Inputs:
            - std::complex<double> x: The value to evaluate the probability density function at.
        
        Outputs:
            - double x: The probability density evaluated at x.
        */
        return (m_fluxRef * std::conj(m_fluxRef)
            + m_fluxRef * std::exp(-2. * I * m_k1 * x)
            + std::conj(m_fluxRef) * std::exp(2. * I * m_k1 * x) + 1.).real();
    }


    void samplePoint(SamplePoint* point)
    {
        /*
        Determines what proportion of electrons can tunnel across the gap.

        Inputs:
            - SamplePoint point: The point to sample.
        */

        m_actualWidth = height(point->x, point->y);
        const double zeta{ -pow((2. * m_bias * ELECTRONMASS
            / (m_actualWidth * HBAR * HBAR)), 1. / 3.) };
        const double alpha_0{ zeta * (m_energy * m_actualWidth
            - m_actualWidth * m_workFunction) / m_bias };
        const double alpha_h{ zeta * (m_energy * m_actualWidth
            - m_actualWidth * m_workFunction
            + m_bias * m_actualWidth) / m_bias };


        using namespace boost::math;
        m_ai0 = airy_ai(alpha_0);
        m_aih = airy_ai(alpha_h);
        m_bi0 = airy_bi(alpha_0);
        m_bih = airy_bi(alpha_h);


        const double ai0_der{ airy_ai_prime(alpha_0) };
        const double aih_der{ airy_ai_prime(alpha_h) };
        const double bi0_der{ airy_bi_prime(alpha_0) };
        const double bih_der{ airy_bi_prime(alpha_h) };

        m_M_det = m_ai0 * m_bih - m_aih * m_bi0;



        const double kappa11{ (ai0_der * m_bih - m_aih * bi0_der)
            / m_M_det * zeta };
        const double kappa12{ (m_ai0 * bi0_der - ai0_der * m_bi0)
            / m_M_det * zeta };
        const double kappa21{ (aih_der * m_bih - m_aih * bih_der)
            / m_M_det * zeta };
        const double kappa22{ (m_ai0 * bih_der - aih_der * m_bi0)
            / m_M_det * zeta };


        // Calculate flux reflected and transmitted
        m_fluxRef = -(kappa12 * kappa21 + I * m_k1 * kappa22
            - kappa11 * kappa22 + m_k1 * m_k3 + I * kappa11 * m_k3)
            / (kappa12 * kappa21 - I * m_k1 * kappa22
                - kappa11 * kappa22 - m_k1 * m_k3 + I * kappa11 * m_k3);
        m_fluxOut = -2. * m_k1 * kappa21
            * std::exp(-I * m_k3 * m_actualWidth)
            / (I * kappa12 * kappa21 + m_k1 * kappa22
                - I * kappa11 * kappa22 - I * m_k1 * m_k3 - kappa11 * m_k1);

        const double fluxMag{ sqrt((m_fluxRef * std::conj(m_fluxRef)
            + m_fluxOut * std::conj(m_fluxOut)).real()) };

        m_fluxRef /= fluxMag;
        m_fluxOut /= fluxMag;

        // Determine the maximum probability density.
        const double maxProb{pow(1 + abs(m_fluxRef), 2.) };
        
        // Sample electrons
        unsigned int numTunnelled{ 0 };
        // intended integer division
        unsigned int samplingPeriod{ m_maxNumElectrons / m_numReadings };
        
        for (int i{ 0 }; i < m_maxNumElectrons; ++i)
        {
            if (measureElectronPosition(maxProb) > m_actualWidth)
            {
                ++numTunnelled;
            }

            if ((i+1) % samplingPeriod == 0 && i!=0)
            {
                point->tunnelProb[((i+1) / samplingPeriod) - 1] = static_cast<double>(
                    numTunnelled) / static_cast<double>(i+1);

            }
        }

 
    }


    double measureElectronPosition(double maxProb) const
    {
        /*
        Measures the system by collapsing the wavefunction to a single point.

        Inputs:
            - double maxProb: The peak probability density.
        */

        const double lowerBound{ -PI/m_k1 };
        const double upperBound{ m_GAMMA * m_actualWidth };
        double trial{};

        while (true)
        {
            trial = (upperBound - lowerBound) * Random::uniform() + lowerBound;

            if (maxProb * Random::uniform() < probWavefunction(trial))
            {
                return trial;
            }
        }

    }


    std::complex<double> wavefunction(const double x) const
    {
        /*
        The value of the wavefunction evaluated at x.

        Inputs:
            - const double x: The spatial point to evaluate the wavefunction at.

        Outputs:
            - std::complex<double> psi: The wavefunction value.
        */


        // Approaching barrier.
        if (x < 0)
        {
            return std::exp(I * m_k1 * x) + m_fluxRef * std::exp(-I * m_k1 * x);
        }

        // After barrier
        else if (x > m_actualWidth)
        {
            return m_fluxOut * std::exp(I * m_k3 * x);
        }

        // Inside barrier
        else
        {
            double alpha{ -pow((2. * m_bias * ELECTRONMASS / (m_actualWidth * HBAR*HBAR)), 1. / 3.)
                * (m_energy * m_actualWidth - m_actualWidth * m_workFunction + m_bias * x) / m_bias };

            const double ai{ boost::math::airy_ai(alpha) };
            const double bi{ boost::math::airy_bi(alpha) };

            const std::complex<double> c1{ ((1. + m_fluxRef) * m_bih - m_bi0 * m_fluxOut 
                * std::exp(I * m_k1 * m_actualWidth)) / m_M_det };
            const std::complex<double> c2{ (m_ai0 * m_fluxOut
                * std::exp(I * m_k1 * m_actualWidth) - (1. + m_fluxRef) * m_aih) / m_M_det };

            return c1 * ai + c2 * bi;
        }
    }

    double probWavefunction(double x) const
    {
        /*
        The probability density of the wavefunction for a given x.

        Inputs:
            - double x: The spatial point to evaluate at.

        Outputs:
            - double P_x: The probability density.
        */


        return pow(std::abs(wavefunction(x)),2.);
    }

    double get_workFunction() const
    {
        return m_workFunction;
    }
    double get_bias() const
    {
        return m_bias;
    }
    double get_energy() const
    {
        return m_energy;
    }



};



void saveSystemConditions(const unsigned int numSamples, const unsigned int maxNumElectrons,
    const unsigned int numReadings, const double lowerBound, const double upperBound,
    const double workFunction, const double bias, const double energy)
{
    std::ofstream myfile;
    myfile.open("systemConditions.txt");

    myfile << "numSamples:" << numSamples;
    myfile << "\nmaxNumElectrons:" << maxNumElectrons;
    myfile << "\nnumReadings:" << numReadings;
    myfile << "\nlowerBound[m]:" << lowerBound;
    myfile << "\nupperBound[m]:" << upperBound;
    myfile << "\nworkFunction[eV]:" << workFunction;
    myfile << "\nbias[V]:" << bias;
    myfile << "\nenergy[eV]:" << energy;
    

    myfile.close();
}


int main()
{
    std::cout << "START\n";

    constexpr unsigned int numSamples{ 10 };
    constexpr unsigned int maxNumElectrons{ static_cast<unsigned int>(1e+4) };
    constexpr unsigned int numReadings{30};

    // Set system size
    constexpr double lowerBound{ 0. };
    constexpr double upperBound{ 2e-10 };


    STM sys{ numSamples, lowerBound, upperBound, maxNumElectrons, numReadings };

    Timer t{};
    sys.sampleLoop();
    std::cout << "Time taken: " << t.elapsed()
        << " seconds\n";


    saveSystemConditions(numSamples, maxNumElectrons,
        numReadings, lowerBound, upperBound, sys.get_workFunction(), sys.get_bias(), sys.get_energy());


    std::cout << "END\n";


    return 0;
}

