
/// @brief  Formulas for the calculation of the parameters of the prismatic spectrometer
/// @author Lorenzo Liuzzo (lorenzoliuzzo@outlook.com)
/// @date   2022/11/29

#pragma once 
#include "physim.hpp"


using namespace physics::measurements;
using namespace math::tools::descriptive_statistics;


constexpr double radians(double degrees) 
{

    return degrees * math::constants::pi / 180.0;

}
    
    
constexpr uncertain_measurement calculate_n(const uncertain_measurement& delta,
                                            const uncertain_measurement& alpha) 
{

    return sin((delta + alpha) / 2.) / sin(alpha / 2.);

}


constexpr uncertain_measurement calculate_A(const uncertain_measurement& lambda1,
                                            const uncertain_measurement& lambda2,
                                            const uncertain_measurement& n1,
                                            const uncertain_measurement& n2) 
{

    return (square(lambda2 * n2) - square(lambda1 * n1)) / (square(lambda2) + square(lambda1));

}


constexpr uncertain_measurement calculate_B(const uncertain_measurement& lambda1,
                                            const uncertain_measurement& lambda2,
                                            const uncertain_measurement& n1,
                                            const uncertain_measurement& n2) 
{

    return (square(n2) - square(n1)) / (square(lambda2) + square(lambda1));

}


