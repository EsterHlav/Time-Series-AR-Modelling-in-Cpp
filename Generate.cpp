#include <iostream>
#include <iterator>
#include <random>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <fstream>

using namespace std;

/*
 Code for generating the random noise of normal distribution
 Box-Muller transform is a random sampling method for generating pairs of independent, standard, normally distributed (zero expectation, unit variance) random numbers, i.e. transforming random numbers given uniform distribution to normal distribution.
 */
double generateGaussianNoise(const double MEAN, const double SIGMA)
{
    static const double noise = numeric_limits<double>::min();
    static const double two_pi = 2.0*3.14159265358979323846;
    
    double u1, u2;
    do
    {
        u1 = rand() * (1.0 / RAND_MAX);
        u2 = rand() * (1.0 / RAND_MAX);
    }
    while ( u1 <= noise );
    
    double z;
    z = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    return z * SIGMA + MEAN;
}


/*
 Function "generate" for generating aritificial data important for checking the validity and quality of function, i.e. if we have data generated from an AR process, then our function should fit this AR process and we should be able to find back the phis and order
 
 Phi is the the list of all the coefficients in the AR model -- the number of coefficients is stored in the variable ORDER.
 x is an array of the first ORDER data points. N is the number of x to generate (x1,..,xN)
 
 Example -> AR(3) $X_t = \phi_1 * X_t-1 + \phi_2 * X_t-2 + \phi_3 * X_t-3$
 */
double* generate(double* phis, double* x, const int ORDER, const int N, const double SIGMA)
{
    double* arrXs = new double[N]; // dynamically allocating memory for array of Xs
    
    for (int index = 0; index < ORDER; index++)
    {
        *(arrXs+index) = *(x+index);
    }
    
    // computing X of N+j, i.e. all Xs given by N
    for (int t = ORDER; t < N; t++)
    {
        double sum = 0;
        
        // here computing X of j
        for (int i=0; i < ORDER; i++)
        {
            sum += ( *(phis+i) * (*(arrXs+t-i-1)) );
        }
        
        // adding random noise from normal distribution with mean 0 and standard deviation SIGMA
        // SIGMA is an input variable that specifies how much noise there should be in the data
        // I am using the library for normal distribution to generate the noise based on the choice of SIGMA
        
        // Noise epsilon for Gaussian distribution
        double epsilon = generateGaussianNoise(0.0, SIGMA);
        
        sum += epsilon;
        *(arrXs+t) = sum;
    }
    
    
    return arrXs;
}

