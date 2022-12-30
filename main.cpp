#include <iostream>
#include <random>
#include <functional>
#include <cmath>
#include <fstream>

#define PI 3.14159265358979323846
#define C_LIGHT 29.9  //speed of light cm/ns
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


double henyey_greenstein_F(double theta)
{
    // Calculate the CDF for the Henyey-Greenstein phase function
    // g is the asymmetry parameter
    double g = 0.1;
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



// convert a direction versor to spherical coordinates
void toSpherical(const Vector& v, double& r, double& theta, double& phi) {
  r = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
  theta = std::acos(v.z / r);
  phi = std::atan2(v.y, v.x);
}

// convert spherical coordinates to a direction versor
Vector fromSpherical(double r, double theta, double phi) {
  Vector v;
  v.x = r * std::sin(theta) * std::cos(phi);
  v.y = r * std::sin(theta) * std::sin(phi);
  v.z = r * std::cos(theta);
  return v;
}

// compute the output versor given the input direction versor, the deflection angle, and the azimuth angle
Vector computeOutputVersor(const Vector& v, double deflectionAngle, double azimuthAngleDeflection) {
  // convert input direction versor to spherical coordinates
  double r, theta, phi;
  toSpherical(v, r, theta, phi);
  phi+=azimuthAngleDeflection;
  // modify the polar angle to reflect the deflection angle
  theta += deflectionAngle;

  // convert modified spherical coordinates to output direction versor
  return fromSpherical(r, theta, azimuthAngleDeflection);
}

Vector generateRandomDirection(std::mt19937& rng, const Vector& incomingDirection) {
  // Generate the scattering angle using the inverse CDF of the provided probability function
  double deflectionAngle = inverse_transform_sampling(henyey_greenstein_F, rng);//è il coseno
  std::uniform_real_distribution<double> uniform(0, 1);  // uniform distribution from 0 to 1
  double u = uniform(rng);
  // Compute the scattered direction using the incoming direction and the scattering angle
  /*
  double cosTheta = cos(theta);
  double sinTheta = sqrt(1 - cosTheta * cosTheta);
  double phi =  2 * PI * u;//atan2(incomingDirection.y, incomingDirection.x);//TODO errore qui

  Vector v;
  v.x = sinTheta * cos(phi);
  v.y = sinTheta * sin(phi);
  v.z = cosTheta;
  */
  return computeOutputVersor(incomingDirection, deflectionAngle, 2 * PI * u);
}


int main() {//eulero stocastico?
  std::ofstream out_file("TPSF.txt");
  if (!out_file) {
    std::cerr << "Error opening file" << std::endl;
    return 1;
  }
  //std::array<double> deflectionAngleList(1e4);
  // Create a random number generator
  std::mt19937 rng(12345);
  double radius_detector = 0.1;

  // Set the parameters for the simulation
  double mu_s = 10;
  double step_size;
  int numScatteringEvents = 10;
  std::vector<int> photon_count(1e5);//TODO cambia
  int tot = 0;
  // Initialize the photon
  
  for( int j = 0; j<1e5; ++j){
    Photon photon;
    photon.position.x = 0;
    photon.position.y = 0;
    photon.position.z = 0;
    photon.direction.x = 0;
    photon.direction.y = 0;
    photon.direction.z = -1;
    double length = 0;
    // Propagate the photon through the medium
    for (int i = 0; i < numScatteringEvents; ++i) {
      // Scatter the photon with probability exp(-mu_s * stepSize)
      double pos_z{photon.position.z};// pos_y{photon.position.y}, pos_x{photon.position.x};
      step_size = -log(1-rand01(rng))/mu_s;
      length+=step_size;
      //std::cout<<step_size<<std::endl;
      photon.direction = generateRandomDirection(rng, photon.direction);
        //std::cout<<photon.direction.x<<" "<<photon.direction.y<<" "<<photon.direction.z<<std::endl;
      

      // Propagate the photon in its current direction
      photon.position.x += step_size * photon.direction.x;
      photon.position.y += step_size * photon.direction.y;
      photon.position.z += step_size * photon.direction.z;
      //std::cout<<photon.position.x<<" "<<photon.position.y<<" "<<photon.position.z<<std::endl;
      //TODO detector interaction
      
      if (pos_z*photon.position.z<0){
        //std::cout<<"in"<<std::endl;
        double t =  -photon.position.z/ photon.direction.z;
        double x_intersection = photon.direction.x*t+ photon.position.x;
        double y_intersection = photon.direction.y*t+ photon.position.y;
        if (x_intersection*x_intersection+y_intersection*y_intersection < radius_detector*radius_detector){
            double time = length/C_LIGHT; // in ns
            //std::cout<<"length:"<<length<<" time"<<time<<std::endl;
            if(int(time*1e4)>1e5){
                std::cout<<"exeptionally late"<<std::endl;
                break;
                }
            ++photon_count[int(time*1e4)];
            ++tot;
            break;
            }
        
        }
    }
    
    std::cout<<j<<std::endl;
  }
  for (int i = 0; i < 1e5; ++i) {
    out_file << photon_count[i] << std::endl;
  }

  out_file.close();
  std::cout<<"tot:"<<tot<<std::endl;
  return 0;
}
