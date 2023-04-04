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
#define SIZE_LIST_ANGLE 1000
#define TIME_LIMIT 3
// Structure to represent a 3D vector
class Vector 
    {
    public:
        __host__ __device__ Vector(const float &x,const float &y, const float &z);
        __device__ Vector(const Vector &v);
        float x;
        float y;
        float z;
        __device__ float operator *(const Vector &a);
        __device__ Vector operator *(const float &a);
        __device__ Vector operator +(const Vector &a);
        __device__ Vector& operator +=(const Vector &a);
        
        //friend std::ostream& operator<<(std::ostream& os, const Vector& v);
    };
    

class Photon
    {
    public:

        Vector position;
        Vector direction;
        __device__ Photon(const Vector &position_,const Vector &direction_, const float &mu_s_);
        float length, time;
        __device__ void propagatePhoton(curandState_t &state);//propagate the photon
        __device__ void generateRandomDirection(curandState_t &state);
    private:
        const float mu_s;
        // generate a random direction of scattering

    };
class Detector
    {
    public:
        __host__ __device__ Detector(const Vector &position, const float radius);
        __device__ bool is_recorded(const Photon &photon, const Vector &previous_position);
    private:
        Vector position;
        float radius;
    };
class Results
    {
    public:
        __host__ Results();
        std::array<int, TIME_LIMIT*CH_PER_UNIT> tcspc;
        std::array<int, SIZE_LIST_ANGLE+1> cos_angle;
    };
// Function to generate a random number between 0 and 1 using a uniform distribution
__device__ int sign(const float &x);
__device__ void find_v1(Vector &v1, const Vector &v0, float theta, float phi);
// Calculate the CDF for the Henyey-Greenstein phase function
__host__ float henyey_greenstein_F(const float &theta, const float &g);
// fill array of angles of scattering
__host__ std::array<float, SIZE_LIST_ANGLE> inverse_transform_sampling(std::function<float( const float &, const float &)> cdf, const float &g);
__global__ void propagation(float mu_s, float g, int *tcspc, Detector *detector, int NUM_PHOTONS, int PHOTON_INTEGRATION);
Results simulate(const float &g, const float &mu_s, Detector &detector, int NUM_PHOTONS, int PHOTON_INTEGRATION);
#endif
