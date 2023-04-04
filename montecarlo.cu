#include "montecarlo.hpp"

//TODO deflectionAngleArray accesso Random blocca
//TODO float to float o float16

__constant__  float deflectionAngleArray_const[SIZE_LIST_ANGLE];

__host__ Results::Results()
    {
    for(std::size_t i = 0; i<this->tcspc.size(); ++i)
        this->tcspc[i] = 0;
    
    for(std::size_t i = 0; i<this->cos_angle.size(); ++i)
        this->cos_angle[i] = 0;
    }
__host__ float henyey_greenstein_F(const float &theta, const float &g)
    {
    if(g == 0)
        return (1-cos(theta))/2;
    // Calculate the CDF for the Henyey-Greenstein phase function - Martelli's book 9.12 old version
    return (1-g*g)/(2*g)*(1/(1-g) -1/sqrt(1+g*g-2*g*cos(theta)));
    }
__host__ std::array<float, SIZE_LIST_ANGLE> inverse_transform_sampling(std::function<float( const float &, const float &)> cdf, const float &g)
    {
    std::array<float, SIZE_LIST_ANGLE> deflectionAngleArray;
    for(std::size_t i = 0; i<SIZE_LIST_ANGLE; ++i)
        {
        float x = 0;    // initial value for x
        float cdf_val = cdf(x, g);    // initial value for CDF(x)

        // Find the value x that satisfies F(x) = u
        while (cdf_val < float(i)/SIZE_LIST_ANGLE)
            {
            x += 1e-4;    // increment x
            cdf_val = cdf(x, g);    // update CDF(x)
            }
       deflectionAngleArray[i] = x;
        }
    return deflectionAngleArray;
    }



__device__ int sign(const float &x)
    {
    if(x>0)
        return 1;
    return -1;
    }
__host__ __device__ Vector::Vector(const float &x_,const float &y_, const float &z_):
    x{x_},y{y_},z{z_}
    {}
__device__ Vector::Vector(const Vector &v):
    x{v.x},y{v.y},z{v.z}
    {}
__device__ Vector Vector::operator *(const float &a)
    {
    Vector v(a*this->x, a*this->y, a*this->z);
    return v;
    }
__device__ float Vector::operator *(const Vector &a)
    {
    return a.x*this->x + a.y*this->y + a.z*this->z;
    }
__device__ Vector Vector::operator +(const Vector &a)
    {
    Vector v(a.x+this->x, a.y+this->y, a.z+this->z);
    return v;
    }
__device__ Vector& Vector::operator +=(const Vector &a)
    {
    this->x += a.x;
    this->y += a.y;
    this->z += a.z;
    return *this;
    }


__device__ Photon::Photon(const Vector &position_,const Vector &direction_, const float &mu_s_):
    position(position_),
    direction(direction_),
    mu_s(mu_s_)
    {
    
    length = 0;
    time = 0;
    }


__host__ __device__ Detector::Detector(const Vector &position_, const float radius_):
    position(position_),
    radius{radius_}{}
__device__ bool Detector::is_recorded(const Photon &photon, const Vector &previous_position)
    {
    if (previous_position.z < 0 && photon.position.z>0){
        float t =  -photon.position.z/photon.direction.z;
        float x_intersection = photon.direction.x*t+ photon.position.x;
        float y_intersection = photon.direction.y*t+ photon.position.y;
        float x = x_intersection - this->position.x;
        float y = y_intersection - this->position.y;
        //std::cout<<x_intersection << " "<<x_det<<" "<<y_intersection << " "<<y_det<<std::endl;
        if (x*x+y*y< this->radius*this->radius)
            return true;
        }
    return false;
    }




__device__ void find_v1(Vector &v1, const Vector &v0, float theta, float phi)
    {
    if(fabs(v0.z)>0.99)
        {
        v1.x= sin(theta)*cos(phi);
        v1.y = sin(theta)*sin(phi);
        v1.z = sign(v0.z)*cos(theta);
        return;
        }
    float c0 = sin(theta)/sqrt(1.-v0.z*v0.z);    //preso da mcx_core.cu riga 931
    v1.x = c0*(v0.x*v0.z*cos(phi)-v0.y*sin(phi))+v0.x*cos(theta);
    v1.y = c0*(v0.y*v0.z*cos(phi)+v0.x*sin(phi))+v0.y*cos(theta);
    v1.z =  -c0*(1.-v0.z*v0.z)*cos(phi)+v0.z*cos(theta);
    //std::cout<< v1*v0<<" "<<cos(theta)<<std::endl;
    }
__device__ void Photon::generateRandomDirection(curandState_t &state)
    {// generate a random direction of scattering
    // Generate the scattering angle using the inverse CDF of the provided probability function
    float theta = deflectionAngleArray_const[int(curand_uniform(&state)*(SIZE_LIST_ANGLE-1))];
    //std::cout<<cos(theta)<<" ";
    float phi = 2 * PI * curand_uniform(&state);
    Vector u(this->direction);
    find_v1(this->direction, u, theta, phi);

    }
__device__ void Photon::propagatePhoton(curandState_t &state)
    {
    float dl = -log(1-curand_uniform(&state))/this->mu_s;//TODO check log
    //float dt = dl/C_LIGHT;
    this->length+=dl;
    generateRandomDirection(state);
    this->position += this->direction*dl;    // Propagate the photon in new direction
    this->time = length/C_LIGHT; // in ns
    }

__global__ void propagation(float mu_s, float g, int *tcspc, Detector *detector, int NUM_PHOTONS, int PHOTON_INTEGRATION)
    {//TODO come inizilizzare una classe Detector
    int tot = 0;
    int repetition_per_thread = NUM_PHOTONS/(blockDim.x*gridDim.x);
    bool finish = false;
    curandState_t state;
    curand_init(blockIdx.x * blockDim.x + threadIdx.x, 0, 0, &state);
    for( int i = 0; i<repetition_per_thread && !finish; ++i) //in kernel
        {
        Vector position_start(0,0,0);
        Vector direction_start(0,0,-1);
        Photon photon(position_start, direction_start, mu_s);
        // Propagate the photon through the medium
        while (photon.time < TIME_LIMIT) 
            {
            Vector position_previous(photon.position);// Store the current position in a temporary variable
            Vector direction_previous(photon.direction);
            photon.propagatePhoton(state);
            Vector direction_new(photon.direction);
            if (detector->is_recorded(photon, position_previous))
                {
                if(int(photon.time*CH_PER_UNIT)>TIME_LIMIT*CH_PER_UNIT) break;
                atomicAdd(&tcspc[int(photon.time*CH_PER_UNIT)],1);
                ++tot;
                if(tot > PHOTON_INTEGRATION/(blockDim.x*gridDim.x)) finish = true;//return tcspc;  
                break;
                }  
            }
        }
    
    }

Results simulate(const float &g, const float &mu_s, Detector &detector, int NUM_PHOTONS, int PHOTON_INTEGRATION)
    {
    std::array<float, SIZE_LIST_ANGLE> deflectionAngleArray;
    // Create a random number generator
    deflectionAngleArray = inverse_transform_sampling(henyey_greenstein_F, g);

    Results res;
    int *tcspc;
    
    cudaEvent_t start, stop;
    check(cudaEventCreate(&start));
    check(cudaEventCreate(&stop));

    Detector *detector_dvc;
    check(cudaMalloc(&detector_dvc, 1*sizeof(Detector)));
    check(cudaMemcpy(detector_dvc,&detector, 1*sizeof(Detector), cudaMemcpyHostToDevice));
    check(cudaMalloc(&tcspc, res.tcspc.size()*sizeof(int)));
    check(cudaMemset(tcspc,0, res.tcspc.size()*sizeof(int)));
    check(cudaMemcpyToSymbol(deflectionAngleArray_const,deflectionAngleArray.data(), SIZE_LIST_ANGLE*sizeof(float)));
    dim3 dimBlock(1024);
    dim3 dimThreads(32);//TODO come settare i numeri??
    check(cudaEventRecord(start, 0));
    propagation<<<dimBlock,dimThreads>>>(mu_s, g, tcspc, detector_dvc, NUM_PHOTONS, PHOTON_INTEGRATION);
    
    cudaDeviceSynchronize();
    CHECK_LAST_CUDA_ERROR()
    check(cudaEventRecord(stop, 0));
    float time;
    check(cudaEventElapsedTime(&time, start, stop));
    printf("elapsed Time %f:", time);
    check(cudaEventDestroy(start));
    check(cudaEventDestroy(stop));
    
    check(cudaMemcpy(res.tcspc.data(),tcspc, res.tcspc.size()*sizeof(int), cudaMemcpyDeviceToHost));
    check(cudaFree(tcspc));
    check(cudaFree(detector_dvc));

    
    return res;
    }
int main(int argc, char* argv[])
    {
    if (argc < 3) 
        {
        std::cerr << "Usage: " << argv[0] << " <NUM_PHOTONS> <PHOTON_INTEGRATION>" << std::endl;
        return 1;
        }

    const int NUM_PHOTONS = std::atoi(argv[1]);
    const int PHOTON_INTEGRATION = std::atoi(argv[2]);
    Vector pos(1,0,0);//è una funzione __device__
    Detector d(pos,0.1);
    simulate(0, 10, d, NUM_PHOTONS, PHOTON_INTEGRATION);
    return 0;
    }
