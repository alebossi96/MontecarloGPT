#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include <iostream>
#include <random>
#include <functional>
#include <cmath>
#include <fstream>
#include <cstdlib>

#define PI 3.14159265358979323846
#define C_LIGHT 29.9  //TODO speed of light cm/ns

#define NUM_PHOTONS 1e9
#define SIZE_LIST_ANGLE 1000
#define TIME_LIMIT 1
#define PHOTON_INTEGRATION 1e5
// Structure to represent a 3D vector
class Vector 
    {
    public:
        Vector(const double &x,const double &y, const double &z);
        Vector(const Vector &v);
        double x;
        double y;
        double z;
        Vector operator *(const double &a);
        Vector operator +(const Vector &a);
        Vector& operator +=(const Vector &a);
        friend std::ostream& operator<<(std::ostream& os, const Vector& v);
    };
    
class VectorSphericalCoordinate
    {
    public:
        VectorSphericalCoordinate(const Vector &v);
        VectorSphericalCoordinate(const double &r, const double &theta, const double &phi);
        double r, theta, phi;
        Vector to_cartesian_coordinates();// convert spherical coordinates to a direction versor
    };

// Structure to represent a photon
class Photon
    {
    public:

        Vector position;
        Vector direction;
        VectorSphericalCoordinate direction_spherical;
        Photon(const Vector &position_,const Vector &direction_, const double &mu_s_);
        double length, time;
        void propagatePhoton(std::mt19937& rng, const std::array<double, SIZE_LIST_ANGLE>& deflectionAngleArray);//propagate the photon
    private:
        const double mu_s;
        void computeOutputVersor(const double &deflectionAngle, const double &azimuthAngleDeflection);
        // generate a random direction of scattering
        void generateRandomDirection(std::mt19937& rng, const std::array<double, SIZE_LIST_ANGLE> &deflectionAngleArray);
    };
class Detector
    {
    public:
        Detector(const Vector &position, const double radius);
        bool is_recorded(const Photon &photon, const Vector &previous_position);
    private:
        Vector position;
        double radius;
    };
// Function to generate a random number between 0 and 1 using a uniform distribution
double rand01(std::mt19937& rng);
// Calculate the CDF for the Henyey-Greenstein phase function
double henyey_greenstein_F(const double &theta, const double &g);
// fill array of angles of scattering
std::array<double, SIZE_LIST_ANGLE> inverse_transform_sampling(std::function<double( const double &, const double &)> cdf, const double &g);

std::vector<int> simulate(const double &g, const double &mu_s);
#endif
