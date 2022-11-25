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

    assert(d.units().base() == base::metre && d.value() > 0.0);
    assert(delta_theta.units() == rad);
    assert(k > 0);

    return d * sin(delta_theta) / k; 

}


uncertain_measurement get_lambda_from_data(const std::string& data_file, 
                                           const uncertain_measurement& passo_reticolo, 
                                           const uncertain_measurement& theta0) 
{

    std::vector<uncertain_measurement> lambda_data; 
    uncertain_measurement delta; 
    int32_t gradi, primi, k;
    std::ifstream data_flow(data_file);

    if (!data_flow.is_open()) 
        throw std::runtime_error("Cannot open file: '" + data_file + "'"); 
    
    while (!data_flow.eof()) {
        
        data_flow >> gradi >> primi >> k;
        delta = uncertain_measurement(radians(gradi + primi / 60.0), radians(1. / 60.), rad) - theta0;
        if (delta.value() < 0.0) delta.value() *= -1.0; 
        lambda_data.emplace_back(calculate_lambda(passo_reticolo, delta, k));

    }

    return wmean(lambda_data);

}


int main() {

    // passo del reticolo
    uncertain_measurement d(3.368531E-06, 2.406966E-09, m);

    // posizione angolare del primo minimo
    uncertain_measurement theta0(radians(56 + 5. / 60.), radians(1. / 60.), rad); 

    // spettro del mercurio
    std::vector<std::string> colors = { "viola_int", 
                                        "viola_est", 
                                        "indaco", 
                                        "verde_acqua_int", 
                                        "verde_acqua_est", 
                                        "verde", 
                                        "giallo_int", 
                                        "giallo_est", 
                                        "rosso_int", 
                                        "rosso_est" };

    std::vector<uncertain_measurement> lambda_data;
    lambda_data.reserve(colors.size());

    for (const auto& color : colors) {

        lambda_data.emplace_back(get_lambda_from_data("../data/mercurio/" + color + ".dat", d, theta0));

    }

    for (size_t i = 0; i < colors.size(); i++) {

        std::cout << "\t" << colors[i] << "\t\t" << lambda_data[i] << "\n";

    }


    return 0;
}