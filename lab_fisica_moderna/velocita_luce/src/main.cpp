/**
 * @file main.cpp
 * @author Lorenzo Liuzzo (lorenzoliuzzo@outlook.com)
 * @brief A program to calculate the speed of light from experiment of Focault
 * @date 2022-12-09
 * 
 * @copyright Copyright (c) 2022
 */


#include "measurements.hpp"


using namespace measurements; 


uncertain_measurement calculate_c_speed(const measurement& f, 
                                        const uncertain_measurement& a, 
                                        const uncertain_measurement& D, 
                                        const uncertain_measurement& delta_omega, 
                                        const uncertain_measurement& delta_sigma){

    return (4 * square(D) * f * delta_omega) / (delta_sigma * (D + a - f)); 

}


int main() {


    uncertain_measurement a; ///< distanza tra lente e specchio rotante
    measurement f; ///< lunghezza focale della lente
    uncertain_measurement D; ///< somma delle distanze tra specchi
    std::vector<std::string> methods{"CW", "CCW", "CW_CCW"}; ///< metodi di misura 

    std::vector<uncertain_measurement> mu0; ///< frequenza iniziale
    std::vector<uncertain_measurement> mu1; ///< frequenza finale
    std::vector<uncertain_measurement> sigma0; ///< posizione iniziale
    std::vector<uncertain_measurement> sigma1; ///< posizione finale
    std::vector<uncertain_measurement> c_speed; ///< velocità della luce


    std::ofstream file_out; 
    std::ifstream file_in; 

    file_in.open("../data/distanza_lente.dat");
    if (!file_in.is_open()) {

        std::cerr << "Error: unable to open file '../data/distanza_lente.dat'\n";
        exit(EXIT_FAILURE);

    } else file_in >> a;
    file_in.close();

    file_in.open("../data/lunghezza_focale.dat");
    if (!file_in.is_open()) {

        std::cerr << "Error: unable to open file '../data/lunghezza_focale.dat'\n";
        exit(EXIT_FAILURE);

    } else file_in >> f;
    file_in.close();
    
    std::vector<uncertain_measurement> distanze_specchi = tools::read_uncertain_measurements("../data/distanze_specchi.dat");
    D = distanze_specchi[0] + distanze_specchi[1] + distanze_specchi[2];
    
    std::cout << "Condizioni iniziali\n"; 
    std::cout << "a = " << a << '\n';
    std::cout << "f = " << f << '\n';
    std::cout << "D = " << D << '\n';


    for (std::string method : methods) {

        mu0 = tools::read_uncertain_measurements("../data/" + method + "/mu0.dat", 5 * Hz); 
        mu1 = tools::read_uncertain_measurements("../data/" + method + "/mu1.dat", 5 * Hz);
        sigma0 = tools::read_uncertain_measurements("../data/" + method + "/sigma0.dat", 0.05 * mm);
        sigma1 = tools::read_uncertain_measurements("../data/" + method + "/sigma1.dat", 0.05 * mm);

        std::size_t N_samples = mu0.size();

        std::vector<uncertain_measurement> delta_omega(N_samples);
        std::vector<uncertain_measurement> delta_sigma(N_samples);
        std::vector<uncertain_measurement> c_speed_data(N_samples);

        file_out.open("../data/" + method + "/c.dat", std::ios::app);
        for (size_t i{}; i < N_samples; i++) {

            delta_omega[i] = 2 * constants::math::pi * abs(mu1[i] - mu0[i]);
            delta_sigma[i] = abs(sigma0[i] - sigma1[i]);
            c_speed_data[i] = calculate_c_speed(f, a, D, delta_omega[i], delta_sigma[i]).convert_to(m / s);
            file_out << c_speed_data[i] << '\n';

        }
        file_out.close();

        c_speed.emplace_back(statistics::weighted_mean(c_speed_data));
        std::cout << "c " << method << " = " << c_speed.back() << '\n';

    }

    std::cout << "c mean = " << statistics::weighted_mean(c_speed) << '\n';


    return 0;

}