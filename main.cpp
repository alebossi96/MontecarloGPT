#include <iostream>
#include <random>
#include <functional>
#include <cmath>
#include <fstream>

#define PI 3.14159265358979323846

// Structure to represent a 3D vector
struct Vector {
    double x;
    double y;
    double z;
};

// Structure to represent a photon
struct Photon {
    Vector position;
    Vector direction;
};

// Function to generate a random number between 0 and 1 using a uniform distribution
double rand01(std::mt19937& rng) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}

// Function to generate a random number with a Gaussian distribution
double randGauss(std::mt19937& rng) {
    std::normal_distribution<double> dist;
    return dist(rng);
}


double henyey_greenstein_F(double theta)
{
    // Calculate the CDF for the Henyey-Greenstein phase function
    // g is the asymmetry parameter
    double g = 0.9;
    //Martelli's book 9.12 old version
    return (1-g*g)/(2*g)*(1/(1-g) -1/sqrt(1+g*g-2*g*cos(theta)));
}

double inverse_transform_sampling(std::function<double(double)> cdf,    // CDF for the distribution
                                    std::mt19937& rng)    // random number generator
{
    std::uniform_real_distribution<double> uniform(0, 1);    // uniform distribution from 0 to 1
    double u = uniform(rng);    // generate a random number from the uniform distribution
    double x = 0;    // initial value for x
    double cdf_val = cdf(x);    // initial value for CDF(x)

    // Find the value x that satisfies F(x) = u
    while (cdf_val < u) {
        x += 0.01;    // increment x
        cdf_val = cdf(x);    // update CDF(x)
    }

    return x;
}



int main() {//eulero stocastico
    // Create a random number generator
    std::ofstream out_file("output.txt");
    if (!out_file) {
        std::cerr << "Error opening file" << std::endl;
        return 1;
    }
    std::mt19937 rng(12345);
    for(int i = 0; i<1e5; ++i){
        double n = inverse_transform_sampling(henyey_greenstein_F, rng);
        //std::cout<<n*180/3.14<<std::endl;
        out_file<<n<<std::endl;
    }
    out_file.close();
    return 0;
}
