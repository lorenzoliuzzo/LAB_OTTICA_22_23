#include "physim.hpp"


using namespace physics::measurements; 
using namespace math::tools::descriptive_statistics;


constexpr double radians(double degrees) 
{

    return degrees * math::constants::pi / 180.0;

}


uncertain_measurement calculate_lambda(const uncertain_measurement& d, 
                                       const uncertain_measurement& delta_theta,
                                       const int32_t& k) 
{

    assert(d.units().base() == base::metre);
    assert(delta_theta.units() == rad);
    assert(k > 0);

    return d * sin(delta_theta) / k;

}


int main() {

    // passo del reticolo
    uncertain_measurement d(3.368531E-06, 2.406966E-09, m);

    // posizione angolare del primo minimo
    uncertain_measurement theta0(radians(56 + 5. / 60.), radians(1. / 60.), rad); 

    // lettura dati 
    std::ifstream data_flow;
    std::vector<uncertain_measurement> theta_viola_int; 
    std::vector<int32_t> k_viola_int;
    int32_t gradi, primi, k;

    data_flow.open("../data/mercurio/viola_int.dat");
    if (!data_flow.is_open()) throw std::runtime_error("Cannot open file!, line: " + std::to_string(__LINE__));
    while (!data_flow.eof()) {

        data_flow >> gradi >> primi >> k;

        theta_viola_int.emplace_back(uncertain_measurement(radians(gradi + primi / 60.), radians(1. / 60.), rad));
        k_viola_int.emplace_back(k);

    }
    std::cout << "Viola int data: " << theta_viola_int.size() << "\n";

    // calcolo lambda
    uncertain_measurement delta;
    std::vector<uncertain_measurement> lambda_viola_int; 
    lambda_viola_int.reserve(theta_viola_int.size());
    for (size_t i = 0; i < theta_viola_int.size(); i++) {
        delta = theta_viola_int[i] - theta0;
        if (delta.value() < 0) delta.value() *= -1; 
        lambda_viola_int.emplace_back(calculate_lambda(d, delta, k_viola_int[i]));
        std::cout << lambda_viola_int[i] << "\n";
    }

    // calcolo media
    uncertain_measurement lambda_mean = wmean(lambda_viola_int);
    std::cout << "lambda viola (I/II) = " << lambda_mean << "\n";


    return 0;
}