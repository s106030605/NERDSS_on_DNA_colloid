#include "io/io.hpp"
#include "tracing.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>

void write_complex_component(long long int simItr, unsigned frameNum, const Parameters& params, const std::vector<Molecule>& moleculeList
    , const std::vector<Complex>& complexList, const Membrane& membraneObject)
{  
    //std::cout<<"IN writing the complex component"<<std::endl; 
    if (complexList.size() > 0){
        std::cout<<"output the complex component!"<<std::endl;

        std::ofstream outfile { std::to_string(frameNum) + "components.txt" };
        //}
        auto printTime = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        outfile << std::left << std::setw(6) << "TITLE" << ' ' << std::left << std::setw(70) << "Components_ID_timestep " << simItr
            << " CREATED " << std::ctime(&printTime);
      for (int i = 0; i < complexList.size(); i++) {

        outfile << "complex_" << i <<std::endl;

        for (auto& memMol : complexList[i].memberList) {

            outfile << memMol << std::endl;

        }
            //outfile  <<std::endl;

      }
      outfile.close();
    }

}
