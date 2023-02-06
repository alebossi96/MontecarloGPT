#include "montecarlo.hpp"

int sign(const double &x)
    {
    if(x>0)
        return 1;
    return -1;
    }
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
Results::Results()
    {
    for(int i = 0; i<CH_PER_UNIT*TIME_LIMIT; ++i)
        tcspc[i] = 0;
    for(int i = 0; i<SIZE_LIST_ANGLE; ++i)
        cos_angle[i] = 0;
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
void find_v1(Vector &v1, const Vector &v0, double theta, double phi)
    {
    if(fabs(v0.z)>0.99)
        {
        v1.x= sin(theta)*cos(phi);
        v1.y = sin(theta)*sin(phi);
        v1.z = sign(v0.z)*cos(theta);
        return;
        }
    double c0 = sin(theta)/sqrt(1.-v0.z*v0.z);    //preso da mcx_core.cu riga 931
    v1.x = c0*(v0.x*v0.z*cos(phi)-v0.y*sin(phi))+v0.x*cos(theta);
    v1.y = c0*(v0.y*v0.z*cos(phi)+v0.x*sin(phi))+v0.y*cos(theta);
    v1.z =  -c0*(1.-v0.z*v0.z)*cos(phi)+v0.z*cos(theta);
    //std::cout<< v1*v0<<" "<<cos(theta)<<std::endl;
    }
void Photon::generateRandomDirection(std::mt19937& rng, const std::array<double, SIZE_LIST_ANGLE> &deflectionAngleArray)
    {// generate a random direction of scattering
    // Generate the scattering angle using the inverse CDF of the provided probability function
    std::uniform_int_distribution<> uniform_int(0, SIZE_LIST_ANGLE); 
    double theta = deflectionAngleArray[std::uniform_int_distribution<std::size_t>(0,SIZE_LIST_ANGLE-1)(rng)];
    //std::cout<<cos(theta)<<" ";
    double phi = 2 * PI * rand01(rng);
    Vector u(this->direction);
    find_v1(this->direction, u, theta, phi);
    //this->computeOutputVersor(theta, azimuth);
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
std::vector<double> test_angle(const double &g, const int &num_sct)
    {
    std::mt19937 rng(12345);
    std::array<double, SIZE_LIST_ANGLE> deflectionAngleArray = inverse_transform_sampling(henyey_greenstein_F, g);
    Vector position_start(0,0,0);
    Vector direction_start(0,0,-1);
    Photon photon(position_start, direction_start, 10);
    std::vector<double> cos_angle;
    for(int i = 0; i<num_sct; ++i)
        {
        /*
        Vector prev(photon.direction);
        photon.generateRandomDirection(rng, deflectionAngleArray);
        Vector succ(photon.direction);
        */
        std::uniform_int_distribution<> uniform_int(0, SIZE_LIST_ANGLE); 
        double theta = deflectionAngleArray[std::uniform_int_distribution<std::size_t>(0,SIZE_LIST_ANGLE-1)(rng)];
        cos_angle.emplace_back(cos(theta));
        }
    return cos_angle;
    }
std::vector<double> test_mus(const double &mu_s, const int &num_sct)
    {
    std::mt19937 rng(12345);
    std::array<double, SIZE_LIST_ANGLE> deflectionAngleArray = inverse_transform_sampling(henyey_greenstein_F, 0);
    Vector position_start(0,0,0);
    Vector direction_start(0,0,-1);
    Photon photon(position_start, direction_start, mu_s);
    std::vector<double> dl;
    double previous_length;
    for(int i= 0; i<num_sct; ++i)
        {
        previous_length = photon.length;
        photon.propagatePhoton(rng, deflectionAngleArray);
        dl.emplace_back(photon.length-previous_length);
        }
    return dl;
    }
Results simulate(const double &g, const double &mu_s, Detector &detector)
    {
    std::array<double, SIZE_LIST_ANGLE> deflectionAngleArray;
    // Create a random number generator
    std::mt19937 rng(12345);
    deflectionAngleArray = inverse_transform_sampling(henyey_greenstein_F, g);
    //int numScatteringEvents = 3;
    Results res;
    //std::vector<double> cos_angle_sct;
    int tot = 0;
    for( int j = 0; j<NUM_PHOTONS && tot < PHOTON_INTEGRATION ; ++j)
        {
        Vector position_start(0,0,0);
        Vector direction_start(0,0,-1);
        Photon photon(position_start, direction_start, mu_s);
        // Propagate the photon through the medium
        while (photon.time < TIME_LIMIT) 
            {
            Vector position_previous(photon.position);// Store the current position in a temporary variable
            Vector direction_previous(photon.direction);
            double time_prev = photon.time;
            photon.propagatePhoton(rng, deflectionAngleArray);
            Vector direction_new(photon.direction);
            double cos_angle = direction_previous*direction_new;
            ++res.cos_angle[int(cos_angle*SIZE_LIST_ANGLE)];
            //std::cout<< cos_angle<<std::endl;
            //cos_angle_sct.emplace_back(cos_angle);
            if (detector.is_recorded(photon, position_previous))
                {
                std::cout<<double(tot)/PHOTON_INTEGRATION<<std::endl;
                
                double dt = - position_previous.z/direction_new.z;
                double time_at_detector = time_prev + dt;
                if(int(time_at_detector*CH_PER_UNIT)>TIME_LIMIT*CH_PER_UNIT) break;
                ++res.tcspc[int(time_at_detector*CH_PER_UNIT)];
                ++tot;
                break;
                }  
            }
        }
    for(int i = 0; i< TIME_LIMIT*CH_PER_UNIT;++i)
        std::cout<<res.tcspc[i]<<std::endl; 
    std::cout<<"tot:"<<tot<<std::endl;
    return res;
    }
