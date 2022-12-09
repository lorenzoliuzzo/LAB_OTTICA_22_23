#include "formulas.hpp"


std::vector<uncertain_measurement> get_angle_data(const std::string& data_file) {

    /// Gradi e primi di arco
    int32_t gradi{}; 
    float primi{}; 

    std::ifstream data_flow(data_file);
    std::vector<uncertain_measurement> theta_data;

    if (!data_flow.is_open()) 
        throw std::runtime_error("Cannot open file: '" + data_file + "'"); 

    while (!data_flow.eof()) {
        
        data_flow >> gradi >> primi;
        theta_data.emplace_back(radians(gradi + primi / 60.), radians(1. / 60.), rad);
    
    }

    std::cout << "Data loaded from file: '" << data_file << "'\n";
    return theta_data; 
    
}


uncertain_measurement get_theta_from_data(const std::string& data_file) 
{
    
    ///< Gradi e primi di arco
    int32_t gradi{}; 
    float primi{}; 

    std::ifstream data_flow(data_file);
    std::vector<measurement> theta_data;

    if (!data_flow.is_open()) 
        throw std::runtime_error("Cannot open file: '" + data_file + "'"); 

    while (!data_flow.eof()) {
        
        data_flow >> gradi >> primi;
        theta_data.emplace_back(radians(gradi + primi / 60.), rad);
    
    }

    return uncertain_measurement(mean(theta_data).value(), radians(1. / 60.), rad);    

}


uncertain_measurement get_sigma_from_data(const std::string& data_file) 
{

    ///< Gradi e primi di arco
    int32_t gradi{}; 
    float primi{}; 

    std::ifstream data_flow(data_file);
    std::vector<uncertain_measurement> sigma_data;

    if (!data_flow.is_open()) 
        throw std::runtime_error("Cannot open file: '" + data_file + "'"); 

    while (!data_flow.eof()) {
        
        data_flow >> gradi >> primi;
        sigma_data.emplace_back(radians(gradi + primi / 60.) - 2 * math::constants::pi, radians(1. / 60.), rad);
    
    }

    return mean(sigma_data);    

}


int main() {

    std::vector<std::string> colors{"giallo",
                                    "giallo2",
                                    "verde",
                                    "indaco",
                                    "viola",
                                    "rosa",
                                    "rosa2"}; ///< Spettro osservato


    std::vector<uncertain_measurement> theta1_data; ///< Primo angolo misurato
    std::vector<uncertain_measurement> theta2_data; ///< Secondo angolo misurato
    std::vector<uncertain_measurement> delta_theta; ///< Differenza theta1 theta2

    std::vector<uncertain_measurement> n_data; ///< Indice di rifrazione del primo strato
    std::vector<uncertain_measurement> A_data; ///< Parametro A
    std::vector<uncertain_measurement> B_data; ///< Parametro B
    std::vector<uncertain_measurement> lambda_isq_data; ///< Lunghezza d'onda
    std::vector<uncertain_measurement> n_expected; ///< Previsioni 

    uncertain_measurement theta0; ///< Posizione angolare del primo minimo di interferenza
    uncertain_measurement theta1; ///< Primo angolo misurato
    uncertain_measurement theta2; ///< Secondo angolo misurato
    uncertain_measurement alpha; ///< Angolo di apertura del prisma
    uncertain_measurement sigma; ///< Angolo di deviazione minima theta
    uncertain_measurement n1; ///< Indice di rifrazione del primo strato
    uncertain_measurement n2; ///< Indice di rifrazione del secondo strato
    uncertain_measurement A; ///< Parametro A
    uncertain_measurement B; ///< Parametro B


    // valori di lunghezza d'onda tabulati 
    lambda_isq_data = read_uncertain_measurements("../data/lambda.dat", nm); 
    for (auto lambda : lambda_isq_data) 
        lambda = 1 / square(lambda); 

    std::cout << "\nDeterminazione angolo di apertura prisma...\n"; 
    theta1_data = get_angle_data("../data/theta1.dat");
    theta2_data = get_angle_data("../data/theta2.dat");
    
    delta_theta.reserve(theta1_data.size());

    for (size_t i{}; i < theta1_data.size(); ++i) 
        delta_theta.emplace_back(theta2_data[i] - theta1_data[i]);

    alpha = abs(math::constants::pi * rad - mean(delta_theta));
    std::cout << "alpha = " << alpha << "\n"; 

    std::cout << "\nDeterminazione posizione angolare iniziale...\n";
    theta0 = get_theta_from_data("../data/theta0.dat");
    std::cout << "theta0 = " << theta0 << "\n";


    std::cout << "\nDeterminazione angolo di deviaizone minima e indici di rifrazione...\n";
    for (auto& color : colors) {

        theta1 = get_sigma_from_data("../data/" + color + ".dat"); 
        sigma = abs(theta1 - theta0); 

        n_data.emplace_back(-calculate_n(sigma, alpha)); 
        std::cout << "n(" << color << ") = " << n_data.back() << "\n";
    
    }

    linear_regression reg_lin; 

    reg_lin.train(lambda_isq_data, n_data); 
    A = reg_lin.slope();
    B = reg_lin.intercept(); 

    for (auto lambda : lambda_isq_data)
        n_expected.emplace_back(reg_lin.predict(lambda)); 

    std::cout << "A = " << A << "\n"; 
    std::cout << "B = " << B << "\n"; 


    return 0; 

}