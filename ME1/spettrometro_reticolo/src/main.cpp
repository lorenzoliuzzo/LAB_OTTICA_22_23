#include "physim.hpp"


using namespace physics::measurements; 
using namespace math::tools::descriptive_statistics;


constexpr double radians(double degrees) 
{

    return degrees * math::constants::pi / 180.0;

}
    

uncertain_measurement calculate_theta(int32_t gradi,
                                      int32_t primi, 
                                      const measurement& beta) 
{

    double theta = radians(gradi + primi / 60.);
    double theta_uncertainty = radians(1. / 60.);
    double sigma_ortho = beta.value() * (1 - std::cos(theta)) / std::cos(theta);

    uncertain_measurement theta_meas(theta, theta_uncertainty, rad);
    theta_meas.add_uncertainty(std::fabs(sigma_ortho));
    return theta_meas;
    
}


uncertain_measurement calculate_d(const measurement& lambda,
                                  const uncertain_measurement& delta_theta,
                                  const int32_t& k) 
{

    assert(lambda >= 0.0 * m);
    assert(delta_theta.units() == rad);
    assert(k > 0);

    return k * lambda / sin(delta_theta); 

}


uncertain_measurement calculate_N(const uncertain_measurement& L, 
                                  const uncertain_measurement& d) {

    assert(L >= 0.0 * m);
    assert(d >= 0.0 * m);

    return L / d.convert_to(L.units());

}


uncertain_measurement calculate_D(const uncertain_measurement& d,
                                  const uncertain_measurement& delta_theta,
                                  const int32_t& k) 
{

    assert(d >= 0.0 * m);
    assert(delta_theta.units() == rad);
    assert(k > 0);

    return k / (d * cos(delta_theta));

}


uncertain_measurement calculate_R(const uncertain_measurement& N,
                                  const int32_t& k) 
{

    assert(N > 0.0);
    assert(k > 0);

    return k * N; 

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


uncertain_measurement get_theta0_from_data(const std::string& data_file, 
                                           const measurement& beta) 
{

    std::ifstream data_flow(data_file);
    int32_t gradi, primi;
    std::vector<measurement> theta0_data;

    if (!data_flow.is_open()) 
        throw std::runtime_error("Cannot open file: '" + data_file + "'"); 

    while (!data_flow.eof()) {
        
        data_flow >> gradi >> primi;
        theta0_data.emplace_back(radians(gradi + primi / 60.), rad);
    
    }

    return uncertain_measurement(mean(theta0_data).value(), radians(1. / 60.), rad);    

}


uncertain_measurement get_lambda_from_data(const std::string& data_file, 
                                           const uncertain_measurement& passo_reticolo, 
                                           const uncertain_measurement& theta0,
                                           const measurement& beta) 
{

    std::vector<uncertain_measurement> lambda_data; 
    uncertain_measurement delta; 
    int32_t gradi, primi, k;
    std::ifstream data_flow(data_file);

    if (!data_flow.is_open()) 
        throw std::runtime_error("Cannot open file: '" + data_file + "'"); 
    
    while (!data_flow.eof()) {
        
        data_flow >> gradi >> primi >> k;
        delta = calculate_theta(gradi, primi, beta) - theta0;
        lambda_data.emplace_back(calculate_lambda(passo_reticolo, abs(delta), k));

    }

    return wmean(lambda_data);

}


int main() {
    

    measurement lambda1_sodio(589.0e-9, m); // prima lunghezza d'onda sodio
    uncertain_measurement lunghezza_reticolo(0.025, 0.001, m); // lunghezza reticolo

    uncertain_measurement theta0; // posizione angolare del primo ordine
    measurement beta(radians(0.22), rad); // fattore di correzione inclinazione del reticolo
    uncertain_measurement d; // passo reticolo
    uncertain_measurement D; // potere dispersivo
    uncertain_measurement N; // numero di fenditure
    uncertain_measurement R; // potere risolutivo

    
    std::vector<uncertain_measurement> d_data; // dati per il calcolo di d
    std::vector<uncertain_measurement> D_data; // dati per il calcolo di D
    std::vector<uncertain_measurement> N_data; // dati per il calcolo di N
    std::vector<uncertain_measurement> R_data; // dati per il calcolo di R
    std::vector<uncertain_measurement> lambda_data; // dati per la media di lambda

    std::vector<std::string> colors = { "viola_int", 
                                        "viola_est", 
                                        "indaco", 
                                        "verde_acqua_int", 
                                        "verde_acqua_est", 
                                        "verde", 
                                        "giallo_int", 
                                        "giallo_est", 
                                        "rosso_int", 
                                        "rosso_est" }; // spettro del mercurio


    // misura posizione angolare del primo minimo
    theta0 = get_theta0_from_data("../data/sodio/theta0.dat", beta); 


    // misura passo del reticolo, numero di fenditure, potere dispersivo e potere risolutivo
    std::ifstream data_flow("../data/sodio/theta.dat");
    int32_t gradi, primi, k;

    if (!data_flow.is_open()) 
        throw std::runtime_error("Cannot open file: '../data/sodio/theta.dat'"); 

    uncertain_measurement delta_theta;
    while (!data_flow.eof()) {
        
        data_flow >> gradi >> primi >> k;
        delta_theta = abs(calculate_theta(gradi, primi, beta)) - theta0;
        d_data.emplace_back(calculate_d(lambda1_sodio, abs(delta_theta), k));
        N_data.emplace_back(calculate_N(lunghezza_reticolo, d_data.back()));
        D_data.emplace_back(calculate_D(d_data.back(), abs(delta_theta), k));
        R_data.emplace_back(calculate_R(N_data.back(), k));

    }
    data_flow.close();

    // d = wmean(d_data);
    d = uncertain_measurement(3.37E-06,	2.4E-09, m);
    N = wmean(N_data);
    D = wmean(D_data);
    R = wmean(R_data);


    // calcolo delle lunghezze d'onda
    lambda_data.reserve(colors.size());
    for (const auto& color : colors) {

        lambda_data.emplace_back(get_lambda_from_data("../data/mercurio/" + color + ".dat", d, theta0, beta));

    }


    // stampa risultati

    std::cout << "Lunghezza reticolo: " << lunghezza_reticolo << "\n";
    std::cout << "Posizione angolare primo minimo: " << theta0 << " rad\n";
    std::cout << "Passo del reticolo: " << d << "\n";
    std::cout << "Numero di fenditure: " << N << "\n";
    std::cout << "Potere dispersivo: " << D.convert_to(rad / m) << "\n";
    std::cout << "Potere risolutivo: " << R << "\n\n";
    std::cout << "Lunghezze d'onda: " << "\n";
    for (size_t i = 0; i < colors.size(); i++) {

        std::cout << colors[i] << "\t\t" << lambda_data[i].convert_to(nm) << "\n";

    }


    return 0;


}