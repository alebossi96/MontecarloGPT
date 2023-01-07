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
    std::vector<int> tcspc;
    std::cout<<"g = "<<atof(argv[1])<<" mu_s = "<< atof(argv[2])<<std::endl;
    Detector detector(Vector(1,0,0),0.1);
    tcspc = simulate(atof(argv[1]), atof(argv[2]), detector);
    for (std::size_t i = 0; i < tcspc.size(); ++i)  out_file << tcspc[i] << std::endl;
    out_file.close();
    return 0;
    }
