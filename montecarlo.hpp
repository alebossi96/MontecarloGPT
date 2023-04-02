#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include <iostream>
#include <random>
#include <functional>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <curand.h>
#include <curand_kernel.h>

#define PI 3.14159265358979323846
#define C_LIGHT 29.9  //TODO speed of light cm/ns

#define CH_PER_UNIT int(1e3)
#define NUM_PHOTONS 1e7
#define SIZE_LIST_ANGLE 1000
#define TIME_LIMIT 3
#define PHOTON_INTEGRATION 1e5
// Structure to represent a 3D vector
class Vector 
    {
    public:
        __host__ __device__ Vector(const double &x,const double &y, const double &z);
        __device__ Vector(const Vector &v);
        double x;
        double y;
        double z;
        __device__ double operator *(const Vector &a);
        __device__ Vector operator *(const double &a);
        __device__ Vector operator +(const Vector &a);
        __device__ Vector& operator +=(const Vector &a);
        
        //friend std::ostream& operator<<(std::ostream& os, const Vector& v);
    };
    

class Photon
    {
    public:

        Vector position;
        Vector direction;
        __device__ Photon(const Vector &position_,const Vector &direction_, const double &mu_s_);
        double length, time;
        __device__ void propagatePhoton(curandState_t &state, double *deflectionAngleArray);//propagate the photon
        __device__ void generateRandomDirection(curandState_t &state, double *deflectionAngleArray);
    private:
        const double mu_s;
        // generate a random direction of scattering

    };
class Detector
    {
    public:
        __host__ __device__ Detector(const Vector &position, const double radius);
        __device__ bool is_recorded(const Photon &photon, const Vector &previous_position);
    private:
        Vector position;
        double radius;
    };
class Results
    {
    public:
        __host__ Results();
        std::array<int, TIME_LIMIT*CH_PER_UNIT> tcspc;
        std::array<int, SIZE_LIST_ANGLE+1> cos_angle;
    };
// Function to generate a random number between 0 and 1 using a uniform distribution
__device__ int sign(const double &x);
__device__ void find_v1(Vector &v1, const Vector &v0, double theta, double phi);
// Calculate the CDF for the Henyey-Greenstein phase function
__host__ double henyey_greenstein_F(const double &theta, const double &g);
// fill array of angles of scattering
__host__ std::array<double, SIZE_LIST_ANGLE> inverse_transform_sampling(std::function<double( const double &, const double &)> cdf, const double &g);
__global__ void propagation(double mu_s, double g, double * deflectionAngleArray, int *tcspc, Detector *detector);
Results simulate(const double &g, const double &mu_s, Detector &detector);
#endif
