#include "montecarlo.hpp"

Vector::Vector(const double &x,const double &y, const double &z)
    {
    this->x = x;
    this->y = y;
    this->z = z;
    }
Vector::Vector(const Vector &v)
    {
    this->x = v.x;
    this->y = v.y;
    this->z = v.z;
    }
Vector Vector::operator *(const double &a)
    {
    Vector v(a*this->x, a*this->y, a*this->z);
    return v;
    }
double Vector::operator *(const Vector &a)
    {
    return a.x*this->x + a.y*this->y + a.z*this->z;
    }
Vector Vector::operator +(const Vector &a)
    {
    Vector v(a.x+this->x, a.y+this->y, a.z+this->z);
    return v;
    }
Vector& Vector::operator +=(const Vector &a)
    {
    this->x += a.x;
    this->y += a.y;
    this->z += a.z;
    return *this;
    }
VectorSphericalCoordinate::VectorSphericalCoordinate(const Vector &v)
    {
    this->r = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    this->theta = std::acos(v.z / r);
    this->phi = std::atan2(v.y, v.x);
    //double z = r * std::cos(theta);
    }
VectorSphericalCoordinate::VectorSphericalCoordinate(const double &r_, const double &theta_, const double &phi_):
    r{r_},
    theta{theta_},
    phi{phi_}{}
Vector VectorSphericalCoordinate::to_cartesian_coordinates()
    {
    double x = r * std::sin(theta) * std::cos(phi);
    double y = r * std::sin(theta) * std::sin(phi);
    double z = r * std::cos(theta);
    //std::cout<<" th  = "<<theta<<" z = "<<z<<std::endl;
    Vector v(x,y,z);
    return v;
    }

std::ostream& operator<<(std::ostream& os, const Vector& v)
    {
    os <<v.x<<" "<<v.y<<" "<<v.z;
    return os;
    }

Photon::Photon(const Vector &position_,const Vector &direction_, const double &mu_s_):
    position(position_),
    direction(direction_),
    direction_spherical(direction_),
    mu_s(mu_s_)
    {
    
    length = 0;
    time = 0;
    }
Detector::Detector(const Vector &position_, const double radius_):
    position(position_),
    radius{radius_}{}
bool Detector::is_recorded(const Photon &photon, const Vector &previous_position)
    {
    if (previous_position.z < 0 && photon.position.z>0){
        double t =  -photon.position.z/photon.direction.z;
        double x_intersection = photon.direction.x*t+ photon.position.x;
        double y_intersection = photon.direction.y*t+ photon.position.y;
        double x = x_intersection - this->position.x;
        double y = y_intersection - this->position.y;
        //std::cout<<x_intersection << " "<<x_det<<" "<<y_intersection << " "<<y_det<<std::endl;
        if (x*x+y*y< this->radius*this->radius)
            return true;
        }
    return false;
    }
double rand01(std::mt19937& rng) 
    {// Function to generate a random number between 0 and 1 using a uniform distribution
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
    }
double henyey_greenstein_F(const double &theta, const double &g)
    {
    if(g == 0)
        return (1-cos(theta))/2;
    // Calculate the CDF for the Henyey-Greenstein phase function - Martelli's book 9.12 old version
    return (1-g*g)/(2*g)*(1/(1-g) -1/sqrt(1+g*g-2*g*cos(theta)));
    }
std::array<double, SIZE_LIST_ANGLE> inverse_transform_sampling(std::function<double( const double &, const double &)> cdf, const double &g)
    {
    std::array<double, SIZE_LIST_ANGLE> deflectionAngleArray;
    for(std::size_t i = 0; i<SIZE_LIST_ANGLE; ++i)
        {
        double x = 0;    // initial value for x
        double cdf_val = cdf(x, g);    // initial value for CDF(x)

        // Find the value x that satisfies F(x) = u
        while (cdf_val < double(i)/SIZE_LIST_ANGLE)
            {
            x += 1e-4;    // increment x
            cdf_val = cdf(x, g);    // update CDF(x)
            }
       deflectionAngleArray[i] = x;
        }
    return deflectionAngleArray;
    }


void Photon::computeOutputVersor(const double &deflectionAngle, const double &azimuthAngleDeflection)
    {// compute the output versor given the input direction versor, the deflection angle, and the azimuth angle
    // convert input direction versor to spherical coordinates
    this->direction_spherical.phi+=azimuthAngleDeflection;
    // modify the polar angle to reflect the deflection angle
    this->direction_spherical.theta += deflectionAngle;
    this->direction  = this->direction_spherical.to_cartesian_coordinates();
    //std::cout<<" d_a = "<<deflectionAngle;
    // convert modified spherical coordinates to output direction versor
    }
void Photon::generateRandomDirection(std::mt19937& rng, const std::array<double, SIZE_LIST_ANGLE> &deflectionAngleArray)
    {// generate a random direction of scattering
    // Generate the scattering angle using the inverse CDF of the provided probability function
    std::uniform_int_distribution<> uniform_int(0, SIZE_LIST_ANGLE); 
    double deflectionAngle = deflectionAngleArray[std::uniform_int_distribution<std::size_t>(0,SIZE_LIST_ANGLE-1)(rng)];
    std::cout<<cos(deflectionAngle)<<" ";
    this->computeOutputVersor(deflectionAngle, 2 * PI * rand01(rng));
    }
void Photon::propagatePhoton(std::mt19937& rng, const std::array<double, SIZE_LIST_ANGLE>& deflectionAngleArray)
    {
    double dl = -log(1-rand01(rng))/this->mu_s;
    //double dt = dl/C_LIGHT;
    this->length+=dl;
    generateRandomDirection(rng, deflectionAngleArray);
    this->position += this->direction*dl;    // Propagate the photon in new direction
    this->time = length/C_LIGHT; // in ns
    }
std::vector<int> simulate(const double &g, const double &mu_s, Detector &detector)
    {
    std::array<double, SIZE_LIST_ANGLE> deflectionAngleArray;
    // Create a random number generator
    std::mt19937 rng(12345);
    deflectionAngleArray = inverse_transform_sampling(henyey_greenstein_F, g);
    //int numScatteringEvents = 3;
    std::vector<int> tcspc(TIME_LIMIT*CH_PER_UNIT);
    std::vector<double> cos_angle_sct;
    int tot = 0;
    for( int j = 0; j<NUM_PHOTONS; ++j)
        {
        Vector position_start(0,0,0);
        Vector direction_start(0,0,-1);
        Photon photon(position_start, direction_start, mu_s);
        // Propagate the photon through the medium
        while (photon.time < TIME_LIMIT) 
            {
            Vector position_previous(photon.position);// Store the current position in a temporary variable
            Vector direction_previous(photon.direction);
            photon.propagatePhoton(rng, deflectionAngleArray);
            Vector direction_new(photon.direction);
            double cos_angle = direction_previous*direction_new;
            std::cout<< cos_angle<<std::endl;
            cos_angle_sct.emplace_back(cos_angle);
            if (detector.is_recorded(photon, position_previous))
                {
                std::cout<<double(tot)/PHOTON_INTEGRATION<<std::endl;
                if(int(photon.time*CH_PER_UNIT)>TIME_LIMIT*CH_PER_UNIT) break;
                ++tcspc[int(photon.time*CH_PER_UNIT)];
                ++tot;
                if(tot > PHOTON_INTEGRATION) goto exit;//return tcspc;  
                break;
                }  
            }
        }
    exit:
    std::ofstream out_file("cos_ang.txt");
    if (!out_file) 
        {
        std::cerr << "Error opening file" << std::endl;
        return tcspc;
        }
    for(std::size_t i = 0; i<cos_angle_sct.size(); ++i)
        {
        out_file<<cos_angle_sct[i]<<std::endl;
        }
    out_file.close();
    std::cout<<"tot:"<<tot<<std::endl;
    return tcspc;
    }
