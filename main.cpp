#include <iostream>
#include <random>
#include <functional>
#include <cmath>
#include <fstream>
#include <cstdlib>

#define PI 3.14159265358979323846
#define C_LIGHT 29.9  //speed of light cm/ns

#define NUM_PHOTONS 1e8
#define SIZE_LIST_ANGLE 1000
#define TIME_LIMIT 2
#define PHOTON_INTEGRATION 1e4
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


double henyey_greenstein_F(double theta, double g)
{
    // Calculate the CDF for the Henyey-Greenstein phase function
    // g is the asymmetry parameter
    //Martelli's book 9.12 old version
    return (1-g*g)/(2*g)*(1/(1-g) -1/sqrt(1+g*g-2*g*cos(theta)));
}

std::array<double, SIZE_LIST_ANGLE> inverse_transform_sampling(std::function<double(double, double)> cdf, double g)    // CDF for the distribution
{
    std::array<double, SIZE_LIST_ANGLE> deflectionAngleArray;
    for(std::size_t i = 0; i<deflectionAngleArray.size(); ++i){
        double x = 0;    // initial value for x
        double cdf_val = cdf(x, g);    // initial value for CDF(x)

        // Find the value x that satisfies F(x) = u
        while (cdf_val < double(i)/SIZE_LIST_ANGLE) {
            x += 1e-4;    // increment x
            cdf_val = cdf(x, g);    // update CDF(x)
        }
       deflectionAngleArray[i] = x;
    }
    return deflectionAngleArray;
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
  return fromSpherical(r, theta, phi);
}

Vector generateRandomDirection(std::mt19937& rng, const Vector& incomingDirection, const std::array<double, SIZE_LIST_ANGLE> deflectionAngleArray) {
  // Generate the scattering angle using the inverse CDF of the provided probability function
  std::uniform_int_distribution<> uniform_int(0, SIZE_LIST_ANGLE); 
  int idx_deflection = uniform_int(rng);
  double deflectionAngle = deflectionAngleArray[idx_deflection];
  //std::cout<<idx_deflection<<" "<<deflectionAngle<<std::endl;
  return computeOutputVersor(incomingDirection, deflectionAngle, 2 * PI * rand01(rng));
}


int main(int argc, char **argv) {//eulero stocastico?
  std::ofstream out_file("TPSF_3.txt");
  if (!out_file) {
    std::cerr << "Error opening file" << std::endl;
    return 1;
  }
  std::array<double, SIZE_LIST_ANGLE> deflectionAngleArray;
  // Create a random number generator
  std::mt19937 rng(12345);
  deflectionAngleArray = inverse_transform_sampling(henyey_greenstein_F, atof(argv[1]));
  /*
  for(std::size_t i = 0; i<deflectionAngleArray.size();++i)
    std::cout<<deflectionAngleArray[i]<<std::endl;
  */
  double radius_detector = 0.1;

  // Set the parameters for the simulation
  double mu_s = atof(argv[2]);
  std::cout<<"g = "<<atof(argv[1])<<" mu_s = "<< atof(argv[2])<<std::endl;
  double step_size;
  //int numScatteringEvents = 3;
  std::vector<int> photon_count(1e4);//TODO cambia
  int tot = 0;
  // Initialize the photon
  double x_det = 1;
  double y_det = 0;
  for( int j = 0; j<NUM_PHOTONS; ++j){
    Photon photon;
    photon.position.x = 0;
    photon.position.y = 0;
    photon.position.z = 0;
    photon.direction.x = 0;
    photon.direction.y = 0;
    photon.direction.z = -1;
    double length = 0;
    double time = 0;
    // Propagate the photon through the medium
    while (time < TIME_LIMIT) {
      // Scatter the photon with probability exp(-mu_s * stepSize)
      double pos_z_previous{photon.position.z};// pos_y{photon.position.y}, pos_x{photon.position.x};
      step_size = -log(1-rand01(rng))/mu_s;
      length+=step_size;
      //std::cout<<step_size<<std::endl;
      photon.direction = generateRandomDirection(rng, photon.direction, deflectionAngleArray);
      //std::cout<<photon.direction.x<<" "<<photon.direction.y<<" "<<photon.direction.z<<std::endl;
      

      // Propagate the photon in its current direction
      photon.position.x += step_size * photon.direction.x;
      photon.position.y += step_size * photon.direction.y;
      photon.position.z += step_size * photon.direction.z;
      //std::cout<<photon.position.x<<" "<<photon.position.y<<" "<<photon.position.z<<std::endl;
      //TODO detector interaction
      time = length/C_LIGHT; // in ns
      if (pos_z_previous < 0 && photon.position.z>0){
        //std::cout<<"in"<<std::endl;
        double t =  -photon.position.z/photon.direction.z;
        double x_intersection = photon.direction.x*t+ photon.position.x;
        double y_intersection = photon.direction.y*t+ photon.position.y;
        double x = x_intersection - x_det;
        double y = y_intersection - y_det;
        //std::cout<<x_intersection << " "<<x_det<<" "<<y_intersection << " "<<y_det<<std::endl;
        if (x*x+y*y< radius_detector*radius_detector){
            std::cout<<tot/PHOTON_INTEGRATION<<std::endl;
            
            //std::cout<<"length:"<<length<<" time"<<time<<std::endl;
            if(int(time*1e4)>1e4){
                break;
                }
            ++photon_count[int(time*1e4)];
            ++tot;
            if(tot > PHOTON_INTEGRATION)
                goto exit;  
            break;
            }
        
        }
    }
    std::cout<<photon.position.x<<" "<<photon.position.y<<" "<<photon.position.z<<std::endl;
    //std::cout<<time<<std::endl;
    
    //std::cout<<j<<std::endl;
  }
  exit:
      for (int i = 0; i < 1e4; ++i) {
        out_file << photon_count[i] << std::endl;
      }

  out_file.close();
  std::cout<<"tot:"<<tot<<std::endl;
  return 0;
}
