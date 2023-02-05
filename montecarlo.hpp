#ifndef MONTECARLO_HPP
#define MONTECARLO_HPP

#include <iostream>
#include <random>
#include <functional>
#include <cmath>
#include <fstream>
#include <cstdlib>

#define PI 3.14159265358979323846
#define C_LIGHT 29.97925  //TODO speed of light cm/ns

#define CH_PER_UNIT 1e3
#define NUM_PHOTONS 1e9
#define SIZE_LIST_ANGLE 10000
#define TIME_LIMIT 2
#define PHOTON_INTEGRATION int(1e4)
// Structure to represent a 3D vector
class Vector 
    {
    public:
        Vector(const double &x,const double &y, const double &z);
        Vector(const Vector &v);
        double x;
        double y;
        double z;
        double operator *(const Vector &a);
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
class Results
    {
    public:
        std::array<int, PHOTON_INTEGRATION*TIME_LIMIT> tcspc;
        std::array<int, SIZE_LIST_ANGLE> cos_angle;
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
        void generateRandomDirection(std::mt19937& rng, const std::array<double, SIZE_LIST_ANGLE> &deflectionAngleArray);
    private:
        const double mu_s;
        void computeOutputVersor(const double &deflectionAngle, const double &azimuthAngleDeflection);
        // generate a random direction of scattering

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
int sign(const double &x);
void find_v1(Vector &v1, const Vector &v0, double theta, double phi);
double rand01(std::mt19937& rng);
// Calculate the CDF for the Henyey-Greenstein phase function
double henyey_greenstein_F(const double &theta, const double &g);
// fill array of angles of scattering
std::array<double, SIZE_LIST_ANGLE> inverse_transform_sampling(std::function<double( const double &, const double &)> cdf, const double &g);
std::vector<double> test_angle(const double &g, const int &num_sct);
std::vector<double> test_mus(const double &mu_s, const int &num_sct);
Results simulate(const double &g, const double &mu_s, Detector &detector);
#endif
