#include "physim.hpp"


using namespace physics::measurements; 
using namespace math::tools::descriptive_statistics; 


int main() {

    std::vector<uncertain_measurement> CW = read_uncertain_measurements("../data/CW.txt", m_s, 20);
    std::vector<uncertain_measurement> CCW = read_uncertain_measurements("../data/CCW.txt", m_s, 15);
    std::vector<uncertain_measurement> CW_CCW = read_uncertain_measurements("../data/CW_CCW.txt", m_s, 3);

    uncertain_measurement CW_mean = weighted_mean(CW);
    uncertain_measurement CCW_mean = weighted_mean(CCW);
    uncertain_measurement CW_CCW_mean = weighted_mean(CW_CCW);

    std::cout << "c_mean in CW = " << CW_mean.value() << " +- " << CW_mean.uncertainty() << " " << CW_mean.units() << "\n"; 
    std::cout << "c_mean in CCW = " << CCW_mean.value() << " +- " << CCW_mean.uncertainty() << " " << CCW_mean.units() << "\n"; 
    std::cout << "c_mean in CW_CCW = " << CW_CCW_mean.value() << " +- " << CW_CCW_mean.uncertainty() << " " << CW_CCW_mean.units() << "\n"; 


}