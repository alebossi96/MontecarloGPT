#include "montecarlo.hpp"

int main(int argc, char **argv) 
    {//eulero stocastico?
    if(argc !=3)
        {
        std::cout<<"put the extra arguments. current argc = "<<argc<<std::endl;
        return 1;
        }
    std::ofstream out_file("TPSF_3.txt");
    if (!out_file) 
        {
        std::cerr << "Error opening file" << std::endl;
        return 1;
        }
    std::array<double, SIZE_LIST_ANGLE> deflectionAngleArray;
    // Create a random number generator
    std::mt19937 rng(12345);
    deflectionAngleArray = inverse_transform_sampling(henyey_greenstein_F, atof(argv[1]));
    Detector detector(Vector(1,0,0),0.2);
    // Set the parameters for the simulation
    std::cout<<"g = "<<atof(argv[1])<<" mu_s = "<< atof(argv[2])<<std::endl;
    //int numScatteringEvents = 3;
    std::vector<int> tcspc(1e4);//TODO cambia
    int tot = 0;
    for( int j = 0; j<NUM_PHOTONS; ++j)
        {
        Vector position_start(0,0,0);
        Vector direction_start(0,0,-1);
        Photon photon(position_start, direction_start, atof(argv[2]));
        // Propagate the photon through the medium
        while (photon.time < TIME_LIMIT) 
            {
            Vector position_previous(photon.position);// Store the current position in a temporary variable
            photon.propagatePhoton(rng, deflectionAngleArray);
            if (detector.is_recorded(photon, position_previous))
                {
                std::cout<<tot/PHOTON_INTEGRATION<<std::endl;
                if(int(photon.time*1e4)<1e4) break;
                ++tcspc[int(photon.time*1e4)];
                ++tot;
                if(tot > PHOTON_INTEGRATION) goto exit;    
                break;
                }  
            }
        //std::cout<<"position = "<<photon.position<< "time = "<<photon.time<<std::endl;
        }
    exit:
        for (int i = 0; i < 1e4; ++i)  out_file << tcspc[i] << std::endl;
    out_file.close();
    std::cout<<"tot:"<<tot<<std::endl;
    return 0;
    }
