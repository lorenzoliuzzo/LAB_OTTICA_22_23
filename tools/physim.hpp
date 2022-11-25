
/**
 * @brief   PHYSIM is a WIP c++ header file that constains the basic tools for computational physics
 * 
 * @author: Lorenzo Liuzzo
 * 
 * @email:  lorenzoliuzzo@outlook.com
 * 
 * @updated: 24/11/2022
 */


#pragma once 

#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip> 
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

#pragma pack(1)


namespace math {
            

    /// @brief Namespace contains some usefull tools for working with double precision floating point numbers
    namespace tools {


        /**
         * @brief check if the two value differ from each other less than the precision given
         * 
         * @param calculated 
         * @param expected 
         * @param epsilon 
         * 
         * @return bool
         */
        template <typename T> 
        constexpr bool are_close(const T& calculated, 
                                const T& expected, 
                                const T& epsilon) {
                                
            return (calculated > expected) ? (calculated - expected) <= epsilon : (expected - calculated) <= epsilon; 

        }


        /**
         * @brief round a value to the expected level of precision of a double
         * 
         * @param value 
         * 
         * @return double 
         */
        double cround(double value) {

            std::uint64_t bits;
            std::memcpy(&bits, &value, sizeof(bits));
            bits += 0x800ULL;
            bits &= 0xFFFFFFFFFFFFF000ULL;
            std::memcpy(&value, &bits, sizeof(bits));
            return value;

        }


        /**
         * @brief rounding compare for equality on double
         * 
         * @param val1 
         * @param val2 
         * 
         * @return bool
         */
        bool compare_round_equals(const double& val1, const double& val2) {

            static constexpr double half_precise_precision{5e-13};
            double v1 = val1 - val2;
            if (v1 == 0.0 || std::fpclassify(v1) == FP_SUBNORMAL) { return true; }
            double c1 = cround(val1);
            double c2 = cround(val2);
            return (c1 == c2) ||
                (cround(val2 * (1.0 + half_precise_precision)) == c1) ||
                (cround(val2 * (1.0 - half_precise_precision)) == c1) ||
                (cround(val1 * (1.0 + half_precise_precision)) == c2) ||
                (cround(val1 * (1.0 - half_precise_precision)) == c2);

        }


        /**
         * @brief Chech for equality between two values
         * 
         * @param val1: double
         * @param val2: double
         * 
         * @return bool
         */
        constexpr bool value_equality_check(const double& val1, const double& val2) { 
            
            return (val1 == val2) ? true : compare_round_equals(val1, val2); 
            
        } 


    } // namespace tools


    /// @brief Namespace contains some math constants
    namespace constants {

        
        constexpr double infinity = std::numeric_limits<double>::infinity();
        
        constexpr double invalid_conversion = std::numeric_limits<double>::signaling_NaN();

        constexpr double pi = 3.14159265358979323846;
        
        constexpr double e = 2.7182818284590452353603;


    } // namespace constants
    

} // namespace math


namespace physics {


    /// @brief Namespace contains some usefull tools for working with physical measurements
    namespace measurements {


        /// @brief Namespace contains some tools for working with units of measurement
        namespace units {


            /// @brief Namespace contains informations for encoding unit base exponents 
            namespace bitwidth {

                
                constexpr uint32_t base_size{sizeof(uint32_t) == 8 ? 8 : 4};
                
                constexpr uint32_t metre{(base_size == 8) ? 8 : 4};
                
                constexpr uint32_t second{(base_size == 8) ? 8 : 4};
                
                constexpr uint32_t kilogram{(base_size == 8) ? 6 : 3};
                
                constexpr uint32_t ampere{(base_size == 8) ? 6 : 3};
                
                constexpr uint32_t candela{(base_size == 8) ? 4 : 2};
                
                constexpr uint32_t kelvin{(base_size == 8) ? 6 : 3};
                
                constexpr uint32_t mole{(base_size == 8) ? 4 : 2};


            } // namespace bitwith


            /// @brief Class representing an unit base using powers of the seven SI unit bases 
            class unit_base {


                public:
                
                // =============================================
                // constructor and destructor
                // =============================================  

                    /// @brief Construct a new default unit base object (default values = 0)
                    explicit constexpr unit_base() noexcept :

                        metre_{0}, second_{0}, kilogram_{0}, 
                        ampere_{0}, kelvin_{0}, mole_{0}, candela_{0} {}


                    /**
                     * @brief Construct a new unit base object with powers of the seven SI unit bases
                     * 
                     * @param metres: power of metre
                     * @param seconds: power of second
                     * @param kilograms: power of kilogram
                     * @param amperes: power of ampere
                     * @param kelvins: power of kelvin
                     * @param moles: power of mole
                     * @param candelas: power of candela
                     */
                    explicit constexpr unit_base(const int& metres, 
                                                 const int& seconds, 
                                                 const int& kilograms, 
                                                 const int& amperes, 
                                                 const int& kelvins,  
                                                 const int& moles, 
                                                 const int& candelas) noexcept : 
                            
                        metre_(metres), second_(seconds), kilogram_(kilograms), 
                        ampere_(amperes), kelvin_(kelvins), mole_(moles), candela_(candelas) {}


                    /**
                     * @brief Construct a new unit base object from a string
                     * 
                     * @param str: string representing the unit base
                     */
                    constexpr unit_base(const std::string& unit_string) noexcept {
                        
                        if (!unit_string.empty()) {

                            size_t finder = unit_string.find("m");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') metre_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') metre_ = - std::stoi(unit_string.substr(finder + 3));
                                    else metre_ = std::stoi(unit_string.substr(finder + 2));
                                }
                                
                            }

                            finder = unit_string.find("s");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') second_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') second_ = - std::stoi(unit_string.substr(finder + 3));
                                    else second_ = std::stoi(unit_string.substr(finder + 2));
                                }
                                
                            }

                            finder = unit_string.find("kg");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') kilogram_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') kilogram_ = - std::stoi(unit_string.substr(finder + 3));
                                    else kilogram_ = std::stoi(unit_string.substr(finder + 2));
                                }

                            }

                            finder = unit_string.find("A");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') ampere_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') ampere_ = - std::stoi(unit_string.substr(finder + 3));
                                    else ampere_ = std::stoi(unit_string.substr(finder + 2));
                                }                            
                            }

                            finder = unit_string.find("K");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') kelvin_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') kelvin_ = - std::stoi(unit_string.substr(finder + 3));
                                    else kelvin_ = std::stoi(unit_string.substr(finder + 2));
                                }

                            }

                            finder = unit_string.find("mol");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') mole_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') mole_ = - std::stoi(unit_string.substr(finder + 3));
                                    else mole_ = std::stoi(unit_string.substr(finder + 2));
                                }

                            }

                            finder = unit_string.find("cd");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') candela_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') candela_ = - std::stoi(unit_string.substr(finder + 3));
                                    else candela_ = std::stoi(unit_string.substr(finder + 2));
                                }
                                
                            }

                        } else {

                            metre_ = 0; second_ = 0; kilogram_ = 0; 
                            ampere_ = 0; kelvin_ = 0; mole_ = 0; candela_ = 0;

                        }

                    }


                    /**
                     * @brief Copy construct a new unit base object from an other unit base object
                     * 
                     * @param other: unit base object to copy
                     */
                    constexpr unit_base(const unit_base& other) noexcept : 

                        metre_(other.metre_), second_(other.second_), 
                        kilogram_(other.kilogram_), ampere_(other.ampere_), 
                        kelvin_(other.kelvin_), mole_(other.mole_), candela_(other.candela_) {}


                    /// @brief Default destructor
                    ~unit_base() = default; 


                // =============================================
                // operators
                // ============================================= 

                    /**
                     * @brief Copy assignment operator
                     * 
                     * @param other: unit_base to copy 
                     * 
                     * @return constexpr unit_base&
                     */
                    constexpr unit_base& operator=(const unit_base& other) noexcept {

                        metre_ = other.metre_; 
                        second_ = other.second_; 
                        kilogram_ = other.kilogram_; 
                        ampere_ = other.ampere_; 
                        kelvin_ = other.kelvin_; 
                        mole_ = other.mole_; 
                        candela_ = other.candela_; 
                        return *this; 

                    }


                    /**
                     * @brief Multiply this unit_base to an unit_base by adding the powers together
                     * 
                     * @param other: unit_base object to multiply with
                     * 
                     * @return constexpr unit_base&
                     */
                    constexpr unit_base& operator*=(const unit_base& other) noexcept {

                        metre_ += other.metre_; 
                        second_ += other.second_; 
                        kilogram_ += other.kilogram_; 
                        ampere_ += other.ampere_; 
                        kelvin_ += other.kelvin_; 
                        mole_ += other.mole_; 
                        candela_ += other.candela_; 
                        return *this; 

                    }


                    /**
                     * @brief Divide this unit_base to an unit_base by subtracting the powers together
                     * 
                     * @param other: unit_base object to divide with
                     * 
                     * @return constexpr unit_base&
                     */
                    constexpr unit_base& operator/=(const unit_base& other) noexcept {

                        metre_ -= other.metre_; 
                        second_ -= other.second_; 
                        kilogram_ -= other.kilogram_; 
                        ampere_ -= other.ampere_; 
                        kelvin_ -= other.kelvin_; 
                        mole_ -= other.mole_; 
                        candela_ -= other.candela_; 
                        return *this; 

                    }


                    /**
                     * @brief Perform a multiplication by adding the powers together
                     * 
                     * @param other: unit_base object to multiply with
                     * 
                     * @return constexpr unit_base 
                     */
                    constexpr unit_base operator*(const unit_base& other) const noexcept {

                        return unit_base(metre_ + other.metre_, 
                                         second_ + other.second_, 
                                         kilogram_ + other.kilogram_, 
                                         ampere_ + other.ampere_, 
                                         kelvin_ + other.kelvin_, 
                                         mole_ + other.mole_, 
                                         candela_ + other.candela_); 

                    }


                    /**
                     * @brief Perform a division by subtracting the powers
                     * 
                     * @param other: unit_base object to divide with
                     * 
                     * @return constexpr unit_base 
                     */
                    constexpr unit_base operator/(const unit_base& other) const noexcept {
                    
                        return unit_base(metre_ - other.metre_,
                                         second_ - other.second_,
                                         kilogram_ - other.kilogram_,
                                         ampere_ - other.ampere_,
                                         kelvin_ - other.kelvin_,
                                         mole_ - other.mole_,
                                         candela_ - other.candela_);

                    }


                    /**
                     * @brief Equality operator
                     * 
                     * @param other: unit_base object to compare with
                     * 
                     * @return bool
                     */
                    constexpr bool operator==(const unit_base& other) const noexcept {
                    
                        return metre_ == other.metre_ && second_ == other.second_ &&
                            kilogram_ == other.kilogram_ && ampere_ == other.ampere_ &&
                            candela_ == other.candela_ && kelvin_ == other.kelvin_ && mole_ == other.mole_;   

                    }


                    /**
                     * @brief Inequality operator
                     * 
                     * @param other: unit_base object to !compare with
                     * 
                     * @return bool
                     */                    
                    constexpr bool operator!=(const unit_base& other) const noexcept { 
                        
                        return !(*this == other); 
                        
                    }


                    /**
                     * @brief Printing on video the unit_base litterals
                     * 
                     * @param os: std::ostream& 
                     * @param base: unit_base as l-value const reference
                     * 
                     * @return std::ostream&
                     */
                    friend std::ostream& operator<<(std::ostream& os, const unit_base& base) noexcept {

                        os << base.to_string(); 
                        return os; 

                    }


                    /**
                     * @brief Printing on file the unit_base litterals
                     * 
                     * @param file: std::ofstream& 
                     * @param base: unit_base as l-value const reference
                     * 
                     * @return std::ofstream&
                     */
                    friend std::ofstream& operator<<(std::ofstream& file, const unit_base& base) noexcept {

                        file << base.to_string(); 
                        return file; 

                    }


                // =============================================
                // operations
                // ============================================= 

                    /**
                     * @brief Invert the unit base
                     * 
                     * @return constexpr unit_base
                     */
                    constexpr unit_base inv() const noexcept {
                    
                        return unit_base(-metre_, 
                                         -second_, 
                                         -kilogram_, 
                                         -ampere_, 
                                         -kelvin_, 
                                         -mole_, 
                                         -candela_);

                    }


                    /**
                     * @brief Take the power of the unit base
                     * 
                     * @param power 
                     * 
                     * @return constexpr unit_base 
                     */
                    constexpr unit_base pow(const int& power) const noexcept { 

                        return unit_base(metre_ * power, 
                                         second_ * power, 
                                         kilogram_ * power, 
                                         ampere_ * power, 
                                         kelvin_ * power, 
                                         mole_ * power, 
                                         candela_ * power);

                    }


                    /**
                     * @brief Take the square unit base
                     * 
                     * @return constexpr unit_base 
                     */
                    constexpr unit_base square() const noexcept { 
                        
                        return unit_base(metre_ * 2, 
                                         second_ * 2, 
                                         kilogram_ * 2, 
                                         ampere_ * 2, 
                                         kelvin_ * 2, 
                                         mole_ * 2, 
                                         candela_ * 2);

                    }


                    /**
                     * @brief  Take the cube unit base
                     * 
                     * @return constexpr unit_base 
                     */
                    constexpr unit_base cube() const noexcept { 
                        
                        return unit_base(metre_ * 3, 
                                         second_ * 3, 
                                         kilogram_ * 3, 
                                         ampere_ * 3, 
                                         kelvin_ * 3, 
                                         mole_ * 3, 
                                         candela_ * 3);

                    }


                    /**
                     * @brief Take the root of the unit base
                     * 
                     * @param power 
                     * 
                     * @return constexpr unit_base 
                     */
                    constexpr unit_base root(const int& power) const {

                        if (!has_valid_root(power)) throw std::invalid_argument("Invalid root power");
                        else return unit_base(metre_ / power,
                                              second_ / power,
                                              kilogram_ / power,
                                              ampere_ / power,
                                              kelvin_ / power,
                                              mole_ / power,
                                              candela_ / power); 

                    }


                    /**
                     * @brief Take the square root of the unit base
                     * 
                     * @return constexpr unit_base 
                     */
                    constexpr unit_base sqrt() const { 
                        
                        return root(2); 
                    
                    }


                    /**
                     * @brief Take the cubic root of the unit base
                     * 
                     * @return constexpr unit_base 
                     */
                    constexpr unit_base cbrt() const { 
                        
                        return root(3);
                    
                    }

                                        
                // =============================================
                // get method
                // =============================================

                    /**
                     * @brief check if tha unit base has a valid root
                     * 
                     * @param power 
                     * 
                     * @return bool 
                     */
                    constexpr bool has_valid_root(const int& power) const {
                    
                        return (metre_ % power == 0) &&
                               (second_ % power == 0) &&
                               (kilogram_ % power == 0) &&
                               (ampere_ % power == 0) &&
                               (candela_ % power == 0) &&
                               (kelvin_ % power == 0) &&
                               (mole_ % power == 0);

                    }      


                    /**
                     * @brief Units litterals to string
                     * 
                     * @return std::string 
                     */
                    std::string to_string() const noexcept {
                        
                        std::string unit_base_string("");   
                        
                        if (metre_ == 1) unit_base_string.append("m");
                        else if (metre_ != 0) unit_base_string.append("m^" + std::to_string(metre_)); 
                        if (second_ == 1) unit_base_string.append("s"); 
                        else if (second_ != 0) unit_base_string.append("s^" + std::to_string(second_)); 
                        if (kilogram_ == 1) unit_base_string.append("kg"); 
                        else if (kilogram_ != 0) unit_base_string.append("kg^" + std::to_string(kilogram_)); 
                        if (ampere_ == 1) unit_base_string.append("A"); 
                        else if (ampere_ != 0) unit_base_string.append("A^" + std::to_string(ampere_)); 
                        if (kelvin_ == 1) unit_base_string.append("K");
                        else if (kelvin_ != 0) unit_base_string.append("K^" + std::to_string(kelvin_)); 
                        if (mole_ == 1) unit_base_string.append("mol"); 
                        else if (mole_ != 0) unit_base_string.append("mol^" + std::to_string(mole_)); 
                        if (candela_ == 1) unit_base_string.append("cd"); 
                        else if (candela_ != 0) unit_base_string.append("cd^" + std::to_string(candela_)); 

                        return unit_base_string; 
                    
                    }


                    /// @brief Print the unit base to the standard output
                    inline void print() const noexcept {
                        /// @brief 
                        std::cout << to_string(); 

                    }


                private: 

                // =============================================
                // class members
                // =============================================


                    /// @brief Metre exponent
                    signed int metre_ : bitwidth::metre;
                    
                    /// @brief Second exponent
                    signed int second_ : bitwidth::second;  
                    
                    /// @brief Kilogram exponent
                    signed int kilogram_ : bitwidth::kilogram;
                    
                    /// @brief Ampere exponent
                    signed int ampere_ : bitwidth::ampere;
                    
                    /// @brief Kelvin exponent
                    signed int kelvin_ : bitwidth::kelvin;
                    
                    /// @brief Mole exponent
                    signed int mole_ : bitwidth::mole;
                    
                    /// @brief Candela exponent
                    signed int candela_ : bitwidth::candela;  
        

                    /// @brief Array with the SI exponents 
                    static constexpr uint32_t bits[7] = { 

                        bitwidth::metre, bitwidth::second, 
                        bitwidth::kilogram, bitwidth::ampere, 
                        bitwidth::kelvin, bitwidth::mole, bitwidth::candela

                    };  


            }; // class unit_base
            static_assert(sizeof(unit_base) == 3, "unit_base is not the correct size");


            /// @brief Class representing an unit prefix using a multiplier (double) and a symbol (char)
            class unit_prefix {


                public: 

                // =============================================
                // constructor and destructor
                // ============================================= 

                    
                    /// @brief Construct a new unit prefix object as default
                    constexpr unit_prefix() noexcept : 

                        multiplier_{1.}, 
                        symbol_{'\0'} {}


                    /**
                     * @brief Create the unit_prefix object from a multiplier and a symbol
                     * 
                     * @param mult: double multiplier for scaling the measurement 
                     * @param symbol: char symbol for the string reppresentation
                     * 
                     * @note mult must be positive (> 0)
                     * 
                     */
                    constexpr unit_prefix(const double& mult, 
                                          const char& symbol) : 

                        multiplier_{mult}, 
                        symbol_{symbol} {

                            if (mult <= 0) throw std::invalid_argument("Unit prefix multiplier must be positive");

                        }
                

                    /**
                     * @brief Destroy the unit prefix object
                     * 
                     */
                    ~unit_prefix() = default; 


                // =============================================
                // operators
                // ============================================= 

                    /**
                     * @brief Copy assignment operator
                     * 
                     * @param other: unit_prefix to copy 
                     * 
                     * @return constexpr unit_prefix&
                     * 
                     */
                    constexpr unit_prefix& operator=(const unit_prefix& other) noexcept {

                        multiplier_ = other.multiplier_; 
                        symbol_ = other.symbol_; 
                        return *this; 

                    }


                    /**
                     * @brief Multiply this unit_prefix to an unit_prefix by multiplying the multipliers together
                     * 
                     * @param other: unit_prefix object to multiply with
                     * 
                     * @return constexpr unit_prefix&
                     * 
                     */
                    constexpr unit_prefix& operator*=(const unit_prefix& other) noexcept {

                        multiplier_ *= other.multiplier_; 
                        return *this; 

                    }


                    /**
                     * @brief Divide this unit_prefix to an unit_prefix by dividing the multipliers together
                     * 
                     * @param other: unit_prefix object to divide with
                     * 
                     * @return constexpr unit_prefix&
                     * 
                     */
                    constexpr unit_prefix& operator/=(const unit_prefix& other) noexcept {

                        multiplier_ /= other.multiplier_; 
                        return *this; 

                    }


                    /**
                     * @brief Perform a multiplication between unit_prefixes 
                     * 
                     * @param other: unit_prefix object to multiply with
                     * 
                     * @return constexpr unit_prefix 
                     * 
                     */
                    constexpr unit_prefix operator*(const unit_prefix& other) const noexcept {

                        return unit_prefix(multiplier_ * other.multiplier_, symbol_); 

                    }


                    /**
                     * @brief Perform a division between unit_prefixes 
                     * 
                     * @param other: unit_prefix object to muldividetiply with
                     * 
                     * @return constexpr unit_prefix 
                     * 
                     */
                    constexpr unit_prefix operator/(const unit_prefix& other) const noexcept {

                        return unit_prefix(multiplier_ / other.multiplier_, symbol_); 

                    }


                    /**
                     * @brief Equality operator
                     * 
                     * @param other: unit_prefix object to compare with
                     * 
                     * @return bool
                     * 
                     */
                    constexpr bool operator==(const unit_prefix& other) const noexcept {
                    
                        return (multiplier_ == other.multiplier_ && symbol_ == other.symbol_);   

                    }


                    /**
                     * @brief Inequality operator
                     * 
                     * @param other: unit_prefix object to !compare with
                     * 
                     * @return bool
                     * 
                     */                    
                    constexpr bool operator!=(const unit_prefix& other) const noexcept { 
                        
                        return !(*this == other); 
                        
                    }


                    /**
                     * @brief Printing on video the unit_prefix litterals
                     * 
                     * @param os: std::ostream& 
                     * @param prefix: unit_prefix as l-value const reference
                     * 
                     */
                    friend std::ostream& operator<<(std::ostream& os, const unit_prefix& prefix) noexcept {

                        os << prefix.symbol_; 
                        return os; 

                    }


                    /**
                     * @brief Printing on file the unit_prefix litterals
                     * 
                     * @param file: std::ofstream& 
                     * @param prefix: unit_prefix as l-value const reference
                     * 
                     */
                    friend std::ofstream& operator<<(std::ofstream& file, const unit_prefix& prefix) noexcept {

                        file << prefix.symbol_; 
                        return file; 

                    }


                // =============================================
                // operations
                // ============================================= 

                    /**
                     * @brief Invert the unit_prefix
                     * 
                     * @return constexpr unit_prefix 
                     * 
                     */
                    constexpr unit_prefix inv() const noexcept { 
                        
                        return unit_prefix(1. / multiplier_, symbol_); 
                        
                    }


                    /**
                     * @brief Take the power of the unit_prefix
                     * 
                     * @param power 
                     * 
                     * @return constexpr unit_prefix 
                     * 
                     */
                    constexpr unit_prefix pow(const int& power) const noexcept { 
                        
                        return unit_prefix(std::pow(multiplier_, power), symbol_); 
                        
                    }


                    /**
                     * @brief Take the square of the unit_prefix
                     * 
                     * @return constexpr unit_prefix 
                     * 
                     */
                    constexpr unit_prefix square() const noexcept { 
                        
                        return unit_prefix(std::pow(multiplier_, 2), symbol_); 
                        
                    }


                    /**
                     * @brief Take the cube of the unit_prefix
                     * 
                     * @return constexpr unit_prefix 
                     * 
                     */
                    constexpr unit_prefix cube() const noexcept { 
                        
                        return unit_prefix(std::pow(multiplier_, 3), symbol_); 
                        
                    }


                    /**
                     * @brief Take the root of the unit_prefix
                     * 
                     * @param power 
                     * 
                     * @return constexpr unit_prefix 
                     * 
                     */
                    constexpr unit_prefix root(const int& power) const noexcept { 
                        
                        return unit_prefix(std::pow(multiplier_, power), symbol_); 
                        
                    }


                    /**
                     * @brief Take the square root of the unit_prefix
                     * 
                     * @return constexpr unit_prefix 
                     * 
                     */
                    constexpr unit_prefix sqrt() const noexcept { 
                        
                        return unit_prefix(std::pow(multiplier_, 1. / 2.), symbol_); 
                        
                    }


                    /**
                     * @brief Take the cube root of the unit_prefix
                     * 
                     * @return constexpr unit_prefix 
                     * 
                     */
                    constexpr unit_prefix cbrt() const noexcept { 
                        
                        return unit_prefix(std::pow(multiplier_, 1. / 3.), symbol_); 
                        
                    }


                // =============================================
                // get methods
                // ============================================= 

                    /**
                     * @brief Get the multiplier of the unit_prefix
                     * 
                     * @return constexpr double 
                     * 
                     */
                    constexpr double multiplier() const noexcept {

                        return multiplier_; 

                    }

                    
                    /**
                     * @brief Get the symbol of the unit_prefix
                     * 
                     * @return constexpr char
                     * 
                     */
                    constexpr char symbol() const noexcept {

                        return symbol_; 

                    }


                private:    

                    /// @brief multiplier of the unit_prefix
                    double multiplier_; 

                    /// @brief symbol of the unit_prefix
                    char symbol_; 


            }; // class unit_prefix
            static_assert(sizeof(unit_prefix) == 9, "unit_prefix is not the correct size");


            /// @brief Class representing an unit of measurement as an unit_base and an unit_prefix
            class unit {


                public:

                // =============================================
                // constructors & destructor
                // ============================================= 

                    /// @brief Construct a new default unit object
                    explicit constexpr unit() noexcept : 
                        
                        base_{}, 
                        prefix_{} {}


                    /**
                     * @brief Construct a new unit object from a prefix and an unit_base 
                     * 
                     * @param prefix: double as l-value const reference
                     * @param base: unit_base as l-value const reference
                     * 
                     */
                    explicit constexpr unit(const unit_prefix& prefix, 
                                            const unit_base& base) noexcept : 
                        
                        base_{base}, 
                        prefix_{prefix} {}


                    /**
                     * @brief Construct a new unit object from a unit_prefix and an unit_base
                     * 
                     * @param prefix: unit_prefix as r-value
                     * @param base: unit_base as r-value
                     * 
                     */
                    explicit constexpr unit(unit_prefix&& prefix, 
                                            unit_base&& base) noexcept : 
                        
                        base_{std::move(base)}, 
                        prefix_{std::move(prefix)} {}
                        

                    /**
                     * @brief Construct a new unit object from an unit_base and an unit_prefix
                     * 
                     * @param base: unit_base object as l-value const reference
                     * @param prefix: unit_prefix as l-value const reference
                     * 
                     */
                    explicit constexpr unit(const unit_base& base, 
                                            const unit_prefix& prefix = unit_prefix()) noexcept : 
                        
                        base_{base}, 
                        prefix_{prefix} {}


                    /**
                     * @brief Construct a new unit object from a unit_base and unit_prefix
                     * 
                     * @param base: unit_base as r-value
                     * @param prefix: unit_prefix as r-value
                     * 
                     */
                    explicit constexpr unit(unit_base&& base, 
                                            unit_prefix&& prefix) noexcept : 
                        
                        base_{std::move(base)}, 
                        prefix_{std::move(prefix)} {}


                    /**
                     * @brief Construct a new unit object from an unit_prefix and an unit
                     * 
                     * @param prefix: unit_prefix as l-value const reference
                     * @param unit: unit as l-value const reference
                     * 
                     */
                    explicit constexpr unit(const unit_prefix& prefix, 
                                            const unit& unit) noexcept : 
                        
                        base_{unit.base_}, 
                        prefix_{prefix * unit.prefix_} {}


                    /**
                     * @brief Copy construct a new unit object from another unit
                     * 
                     * @param other: unit object to copy
                     * 
                     */
                    constexpr unit(const unit& other) noexcept : 

                        base_{other.base_}, 
                        prefix_{other.prefix_} {}


                    /**
                     * @brief Move construct a new unit object from another unit
                     * 
                     * @param other: unit object to move
                     * 
                     */
                    constexpr unit(unit&& other) noexcept : 

                        base_{std::move(other.base_)}, 
                        prefix_{std::move(other.prefix_)} {}


                    /// @brief Destroy the unit object
                    ~unit() = default; 


                // =============================================
                // operators
                // ============================================= 

                    /**
                     * @brief Copy assignment operator
                     * 
                     * @param other: unit object to copy
                     * 
                     * @return constexpr unit&
                     * 
                     */
                    constexpr unit& operator=(const unit& other) noexcept {

                        base_ = other.base_;
                        prefix_ = other.prefix_;
                        return *this;

                    }


                    /**
                     * @brief Move assignment operator
                     * 
                     * @param other: unit object to move
                     * 
                     * @return constexpr unit& 
                     * 
                     */
                    constexpr unit& operator=(unit&& other) noexcept {

                        base_ = std::move(other.base_);
                        prefix_ = std::move(other.prefix_);
                        return *this;

                    }


                    /**
                     * @brief Multiply this unit with an unit
                     * 
                     * @param other: unit object to multiply with as l-value
                     *                      
                     * @return constexpr unit&
                     * 
                     */
                    constexpr unit& operator*=(const unit& other) noexcept { 
                        
                        base_ *= other.base_; 
                        prefix_ *= other.prefix_; 
                        return *this; 
                        
                    }


                    /**
                     * @brief Multiply this unit with an unit
                     * 
                     * @param other: unit object to multiply with as r-value
                     *                      
                     * @return constexpr unit&
                     * 
                     */
                    constexpr unit& operator*=(unit&& other) noexcept { 
                        
                        base_ *= std::move(other.base_); 
                        prefix_ *= std::move(other.prefix_); 
                        return *this; 
                        
                    }


                    /**
                     * @brief Divide this unit with an unit
                     * 
                     * @param other: unit object to divide with as l-value
                     *                      
                     * @return constexpr unit&
                     * 
                     */
                    constexpr unit& operator/=(const unit& other) noexcept { 
                        
                        base_ /= other.base_; 
                        prefix_ /= other.prefix_; 
                        return *this; 
                        
                    }                 


                    /**
                     * @brief Divide this unit with an unit
                     * 
                     * @param other: unit object to divide with as r-value
                     *                      
                     * @return constexpr unit&
                     * 
                     */
                    constexpr unit& operator/=(unit&& other) noexcept { 
                        
                        base_ /= std::move(other.base_); 
                        prefix_ /= std::move(other.prefix_); 
                        return *this; 
                        
                    }  


                    /**
                     * @brief Perform a multiplication by multiply the bases and the prefixes
                     * 
                     * @param other: unit object to multiply with
                     *                      
                     * @return constexpr unit
                     * 
                     */
                    constexpr unit operator*(const unit& other) const noexcept { 
                        
                        return unit(base_ * other.base_, prefix_ * other.prefix_); 
                        
                    }


                    /**
                     * @brief Perform a division by divide the bases and the prefixes
                     * 
                     * @param other: unit object to divide with
                     *                      
                     * @return constexpr unit
                     * 
                     */
                    constexpr unit operator/(const unit& other) const noexcept { 
                        
                        return unit(base_ / other.base_, prefix_ / other.prefix_);
                        
                    }                    


                    /**
                     * @brief Equality operator
                     * 
                     * @param other: unit object to compare
                     * 
                     * @return constexpr bool 
                     * 
                     */
                    constexpr bool operator==(const unit& other) const noexcept {

                        if (base_ != other.base_) return false; 
                        return prefix_ == other.prefix_; 

                    }


                    /**
                     * @brief Inequality operator
                     * 
                     * @param other: unit object to compare
                     * 
                     * @return constexpr bool 
                     * 
                     */
                    constexpr bool operator!=(const unit& other) const noexcept { 
                        
                        return !operator==(other); 
                        
                    }


                    /**
                     * @brief Output operator for an unit of measurement
                     * 
                     * @param os: std::ostream&
                     * @param units: unit of measurement as l-value const reference
                     * 
                     * @return std::ostream&
                     *  
                     */
                    friend std::ostream& operator<<(std::ostream& os, const unit& units) noexcept {

                        os << units.to_string(); 
                        return os; 

                    }


                    /**
                     * @brief Output operator for an unit of measurement
                     * 
                     * @param file: std::ofstream&
                     * @param units: unit of measurement as l-value const reference
                     * 
                     * @return std::ofstream&
                     *  
                     */
                    friend std::ofstream& operator<<(std::ofstream& file, const unit& units) noexcept {

                        file << units.to_string(); 
                        return file; 

                    }


                // =============================================
                // operations
                // ============================================= 

                    /**
                     * @brief Invert the unit
                     * 
                     * @return constexpr unit
                     * 
                     */
                    constexpr unit inv() const noexcept { 
                        
                        return unit(base_.inv(), prefix_.inv()); 
                        
                    }


                    /**
                     * @brief Take the power of the unit
                     * 
                     * @param power 
                     * 
                     * @return constexpr unit 
                     * 
                     */
                    constexpr unit pow(const int& power) const noexcept { 
                        
                        return unit(base_.pow(power), prefix_.pow(power)); 
                        
                    }


                    /**
                     * @brief Take the square of the unit
                     * 
                     * @return constexpr unit 
                     * 
                     */
                    constexpr unit square() const noexcept { 
                        
                        return unit(base_.square(), prefix_.square()); 
                        
                    }


                    /**
                     * @brief Take the cube of the unit
                     * 
                     * @return constexpr unit 
                     * 
                     */
                    constexpr unit cube() const noexcept { 
                        
                        return unit(base_.cube(), prefix_.cube()); 
                        
                    }


                    /**
                     * @brief Take the root of the unit
                     * 
                     * @param power 
                     * 
                     * @return constexpr unit 
                     * 
                     */
                    constexpr unit root(const int& power) const { 
                        
                        return unit(base_.root(power), prefix_.root(power)); 
                        
                    }


                    /**
                     * @brief Take the square root of the unit
                     * 
                     * @return constexpr unit 
                     * 
                     */
                    constexpr unit sqrt() const { 
                        
                        return unit(base_.sqrt(), prefix_.sqrt()); 
                        
                    }


                    /**
                     * @brief Take the cube root of the unit
                     * 
                     * @return constexpr unit 
                     * 
                     */
                    constexpr unit cbrt() const { 
                        
                        return unit(base_.cbrt(), prefix_.cbrt()); 
                        
                    }


                // =============================================
                // get methods
                // ============================================= 

                    /**
                     * @brief Get the unit base object
                     * 
                     * @return constexpr unit_base 
                     * 
                     */
                    constexpr unit_base base() const noexcept { 
                        
                        return base_; 
                        
                    }


                    /**
                     * @brief Get the unit base object
                     * 
                     * @return constexpr unit_base& 
                     * 
                     */
                    constexpr unit_base& base() noexcept { 
                        
                        return base_; 
                        
                    }


                    /**
                     * @brief Get the unit prefix object
                     * 
                     * @return constexpr unit_prefix 
                     * 
                     */
                    constexpr unit_prefix prefix() const noexcept { 
                        
                        return prefix_; 
                        
                    }


                    /**
                     * @brief Get the unit prefix object
                     * 
                     * @return constexpr unit_prefix&
                     * 
                     */
                    constexpr unit_prefix& prefix() noexcept { 
                        
                        return prefix_; 
                        
                    }


                    /**
                     * @brief Get the unit object
                     * 
                     * @return constexpr unit
                     * 
                     */
                    constexpr unit units() const noexcept { 
                        
                        return *this; 
                        
                    }


                    /**
                     * @brief Get the unit object
                     * 
                     * @return constexpr unit& 
                     */
                    constexpr unit& units() noexcept { 
                        
                        return *this; 
                        
                    }


                    /**
                     * @brief Generate a conversion factor between this unit and another unit 
                     * 
                     * @param result: desired unit
                     * 
                     * @return constexpr double 
                     * @note the units will only convert if they have the same base unit
                     * 
                     */
                    constexpr double convertion_factor(const unit& other) const noexcept { 
                        
                        return (base_ == other.base_) ? prefix_.multiplier() / other.prefix_.multiplier() : math::constants::invalid_conversion; 
                        
                    }


                    /**
                     * @brief Convert a value in the current unit to the target unit 
                     * 
                     * @param value: value to convert 
                     * @param other: desired unit unit
                     * 
                     * @return constexpr double 
                     * @note the units will only convert if they have the same base unit
                     * 
                     */
                    constexpr double convert(const double& value, const unit& other) const noexcept { 
                        
                        return (base_ == other.base_) ? value * prefix_.multiplier() / other.prefix_.multiplier() : math::constants::invalid_conversion; 
                        
                    }


                    /**
                     * @brief Get the unit string
                     * 
                     * @return std::string 
                     */
                    inline std::string to_string() const noexcept {

                        return prefix_.symbol() + base_.to_string(); 

                    }


                    /**
                     * @brief Print the unit
                     * 
                     */
                    inline void print() const noexcept {

                        std::cout << *this;

                    }


                private:
            
                // =============================================
                // class members
                // ============================================= 

                    unit_base base_;

                    unit_prefix prefix_;


            }; // class unit     
            static_assert(sizeof(unit) == 12, "unit is not the correct size");


            #define CREATE_UNIT_BASE(NAME, METRE, SECOND, KILOGRAM, AMPERE, KELVIN, MOLE, CANDELA) \
                constexpr unit_base NAME(METRE, SECOND, KILOGRAM, AMPERE, KELVIN, MOLE, CANDELA); 
            

            #define CREATE_UNIT_PREFIX(NAME, MULTIPLIER, SYMBOL) constexpr unit_prefix NAME(MULTIPLIER, SYMBOL); 


            #define CREATE_UNIT(NAME, PREFIX, BASE) constexpr unit NAME(PREFIX, BASE); 


            /// @brief Namespace containing the units of measurement
            namespace SI {


                /// @brief Namespace with the SI unit bases
                namespace base {


                    CREATE_UNIT_BASE(default_type, 0, 0, 0, 0, 0, 0, 0)

                    CREATE_UNIT_BASE(metre, 1, 0, 0, 0, 0, 0, 0)
                    
                    CREATE_UNIT_BASE(second, 0, 1, 0, 0, 0, 0, 0)
                    
                    CREATE_UNIT_BASE(kilogram, 0, 0, 1, 0, 0, 0, 0)
                    
                    CREATE_UNIT_BASE(Ampere, 0, 0, 0, 1, 0, 0, 0)
                    
                    CREATE_UNIT_BASE(Kelvin, 0, 0, 0, 0, 1, 0, 0)
                    
                    CREATE_UNIT_BASE(mole, 0, 0, 0, 0, 0, 1, 0)
                    
                    CREATE_UNIT_BASE(candela, 0, 0, 0, 0, 0, 0, 1)
                

                } // namespace SI


                /// @brief Namespace containg the SI prefixes
                namespace prefix {


                    CREATE_UNIT_PREFIX(default_type, 1, '\0')

                    CREATE_UNIT_PREFIX(yocto, 1e-24, 'y')
                    
                    CREATE_UNIT_PREFIX(zepto, 1e-21, 'z')
                    
                    CREATE_UNIT_PREFIX(atto, 1e-18, 'a')
                    
                    CREATE_UNIT_PREFIX(femto, 1e-15, 'f')
                    
                    CREATE_UNIT_PREFIX(pico, 1e-12, 'p')
                    
                    CREATE_UNIT_PREFIX(nano, 1e-9, 'n')
                    
                    CREATE_UNIT_PREFIX(micro, 1e-6, 'u')
                    
                    CREATE_UNIT_PREFIX(milli, 1e-3, 'm')
                    
                    CREATE_UNIT_PREFIX(centi, 1e-2, 'c')
                    
                    CREATE_UNIT_PREFIX(deci, 1e-1, 'd')
                                    
                    CREATE_UNIT_PREFIX(hecto, 1e2, 'h')
                    
                    CREATE_UNIT_PREFIX(kilo, 1e3, 'k')
                    
                    CREATE_UNIT_PREFIX(mega, 1e6, 'M')
                    
                    CREATE_UNIT_PREFIX(giga, 1e9, 'G')
                    
                    CREATE_UNIT_PREFIX(tera, 1e12, 'T')
                    
                    CREATE_UNIT_PREFIX(peta, 1e15, 'P')
                    
                    CREATE_UNIT_PREFIX(exa, 1e18, 'E')
                    
                    CREATE_UNIT_PREFIX(zetta, 1e21, 'Z')
                    
                    CREATE_UNIT_PREFIX(yotta, 1e24, 'Y')


                } // namespace prefix


                // unitless 
                constexpr unit unitless(prefix::default_type, base::default_type);

                // SI units
                constexpr unit m(prefix::default_type, base::metre);
                constexpr unit s(prefix::default_type, base::second);
                constexpr unit kg(prefix::default_type, base::kilogram);
                constexpr unit k(prefix::default_type, base::Kelvin);
                constexpr unit A(prefix::default_type, base::Ampere);
                constexpr unit mol(prefix::default_type, base::mole);
                constexpr unit cd(prefix::default_type, base::candela);

                // length units
                constexpr unit ym(prefix::yocto, base::metre);
                constexpr unit zm(prefix::zepto, base::metre);
                constexpr unit am(prefix::atto, base::metre);
                constexpr unit fm(prefix::femto, base::metre);
                constexpr unit pm(prefix::pico, base::metre);
                constexpr unit nm(prefix::nano, base::metre);
                constexpr unit um(prefix::micro, base::metre);
                constexpr unit mm(prefix::milli, base::metre);
                constexpr unit cm(prefix::centi, base::metre);
                constexpr unit dm(prefix::deci, base::metre);
                constexpr unit hm(prefix::hecto, base::metre);
                constexpr unit km(prefix::kilo, base::metre);
                constexpr unit Mm(prefix::mega, base::metre);
                constexpr unit Gm(prefix::giga, base::metre);
                constexpr unit Tm(prefix::tera, base::metre);
                constexpr unit Pm(prefix::peta, base::metre);
                constexpr unit Em(prefix::exa, base::metre);

                // time units
                constexpr unit ys(prefix::yocto, base::second);
                constexpr unit zs(prefix::zepto, base::second);
                constexpr unit as(prefix::atto, base::second);
                constexpr unit fs(prefix::femto, base::second);
                constexpr unit ps(prefix::pico, base::second);
                constexpr unit ns(prefix::nano, base::second);
                constexpr unit us(prefix::micro, base::second);
                constexpr unit ms(prefix::milli, base::second);
                constexpr unit cs(prefix::centi, base::second);
                constexpr unit ds(prefix::deci, base::second);
                constexpr unit hs(prefix::hecto, base::second);
                constexpr unit ks(prefix::kilo, base::second);
                constexpr unit Ms(prefix::mega, base::second);
                constexpr unit Gs(prefix::giga, base::second);
                constexpr unit Ts(prefix::tera, base::second);
                constexpr unit Ps(prefix::peta, base::second);
                constexpr unit Es(prefix::exa, base::second);

                constexpr unit rad(prefix::default_type, base::default_type);
                constexpr unit m_s(prefix::default_type, base::metre / base::second); 
                constexpr unit km_s(prefix::kilo, base::metre / base::second); 

                // composed units
                constexpr unit hertz(s.inv());
                constexpr unit Hz = hertz;

                constexpr unit volt(unit_base(2, -3, 1, -1, 0, 0, 0));
                constexpr unit V = volt;

                constexpr unit newton(kg * m / s);
                constexpr unit N = newton;

                constexpr unit Pa(unit_base(-1, -2, 1, 0, 0, 0, 0));
                constexpr unit pascal = Pa;

                constexpr unit joule(unit_base(2, -2, 1, 0, 0, 0, 0));
                constexpr unit J = joule;

                constexpr unit watt(unit_base(2, -3, 1, 0, 0, 0, 0));
                constexpr unit W = watt;

                constexpr unit coulomb(unit_base(0, 1, 0, 1, 0, 0, 0));
                constexpr unit C = coulomb;

                constexpr unit farad(unit_base(-2, 4, -1, 2, 0, 0, 0));
                constexpr unit F = farad;

                constexpr unit weber(unit_base(2, -2, 1, -1, 0, 0, 0));
                constexpr unit Wb = weber;

                constexpr unit tesla(unit_base(0, -2, 1, -1, 0, 0, 0));
                constexpr unit T = tesla;

                constexpr unit henry(unit_base(2, -2, 1, -2, 0, 0, 0));                    
                constexpr unit H = henry;


            } // namespace SI


        } // namespace units


        using namespace units;
        using namespace units::SI;


        /// @brief A class for representing a physical quantity with a numerical value and an unit
        class measurement {


            public:

            // =============================================                                                                                         
            // constructors & destructor 
            // =============================================  

                /// @brief Construct a new default measurement object with value 0.0 and unitless
                explicit constexpr measurement() noexcept : 

                    value_{0.0}, 
                    units_{} {}


                /**
                 * @brief Construct a new measurement object from a value and unit
                 * 
                 * @param value: double value of the measurement
                 * @param unit: unit of measurement
                 * 
                 */
                explicit constexpr measurement(const double& value, 
                                               const unit& units = unit()) noexcept : 
                                               
                    value_{value}, 
                    units_{units} {}


                /**
                 * @brief Copy constuct a new measurement object from another measurement
                 * 
                 * @param other: the measurement to copy
                 * 
                 */
                constexpr measurement(const measurement& other) noexcept : 
                
                    value_{other.value_}, 
                    units_{other.units_} {}
                    
                
                /**
                 * @brief Move construct a new measurement object from another measurement
                 * 
                 * @param other: the measurement to move from
                 * 
                 */
                constexpr measurement(measurement&& other) noexcept :
                    
                    value_{std::move(other.value_)}, 
                    units_{std::move(other.units_)} {}


                /// @brief Default destructor of the measurement object
                ~measurement() = default;


            // =============================================                                                                                         
            // operators
            // =============================================  
                
                /**
                 * @brief Copy assign a measurement from another measurement
                 * 
                 * @param other: the measurement to copy
                 * 
                 * @return measurement& 
                 * 
                 */
                constexpr measurement& operator=(const measurement& other) noexcept {

                    value_ = other.value_;
                    units_ = other.units_;
                    return *this;

                }     


                /**
                 * @brief Move assign a measurement from another measurement
                 * 
                 * @param other: the measurement to move from
                 * 
                 * @return measurement& 
                 * 
                 */
                constexpr measurement& operator=(measurement&& other) noexcept {

                    value_ = std::move(other.value_);
                    units_ = std::move(other.units_);
                    return *this;

                }


                /**
                 * @brief Add a measurement to this measurement
                 * 
                 * @param other: the measurement to add as l-value
                 * 
                 * @return measurement& 
                 * 
                 */
                constexpr measurement& operator+=(const measurement& other) { 
                    
                    if (units_.base() != other.units_.base()) throw std::invalid_argument("Cannot add measurements with different unit bases");
                    if (units_ != unitless) value_ += other.value_as(units_); 
                    else {
                        value_ += other.value_; 
                        units_ = other.units_; 
                    }                 
                    return *this; 
                
                }


                /**
                 * @brief Add a measurement to this measurement
                 * 
                 * @param other: the measurement to add as r-value
                 * 
                 * @return measurement& 
                 * 
                 */
                constexpr measurement& operator+=(measurement&& other) { 

                    if (units_.base() != other.units_.base()) throw std::invalid_argument("Cannot add measurements with different unit bases");
                    if (units_ != unitless) value_ += std::move(other.value_as(units_));   
                    else {
                        value_ += std::move(other.value_);  
                        units_ = std::move(other.units_); 
                    }                
                    return *this; 
                
                }


                /**
                 * @brief Subtract a measurement to this measurement
                 * 
                 * @param other: the measurement to subtract as l-value
                 * 
                 * @return measurement& 
                 * 
                 */
                constexpr measurement& operator-=(const measurement& other) { 

                    if (units_.base() != other.units_.base()) throw std::invalid_argument("Cannot subtract measurements with different unit bases");
                    if (units_ != unitless) value_ -= other.value_as(units_);   
                    else {
                        value_ -= other.value_;    
                        units_ = other.units_; 
                    }              
                    return *this; 
                
                }


                /**
                 * @brief Subtract a measurement to this measurement
                 * 
                 * @param other: the measurement to subtract as r-value
                 * 
                 * @return measurement& 
                 * 
                 */
                constexpr measurement& operator-=(measurement&& other) { 

                    if (units_.base() != other.units_.base()) throw std::invalid_argument("Cannot subtract measurements with different unit bases");
                    if (units_ != unitless) value_ -= std::move(other.value_as(units_));   
                    else {
                        value_ -= std::move(other.value_);    
                        units_ = std::move(other.units_); 
                    }              
                    return *this;   

                }


                /**
                 * @brief Multiply this measurement and a measurement
                 * 
                 * @param other: the measurement to multiply as l-value
                 * 
                 * @return measurement& 
                 * 
                 */
                constexpr measurement& operator*=(const measurement& meas) noexcept { 

                    value_ *= meas.value_;
                    units_ *= meas.units_;                                    
                    return *this; 
                
                }


                /**
                 * @brief Multiply this measurement and a measurement
                 * 
                 * @param other: the measurement to multiply as r-value
                 * 
                 * @return measurement& 
                 * 
                 */
                constexpr measurement& operator*=(measurement&& meas) noexcept { 

                    value_ *= std::move(meas.value_);
                    units_ *= std::move(meas.units_);                                    
                    return *this; 
                
                }


                /**
                 * @brief Divide this measurement and a measurement
                 * 
                 * @param other: the measurement to divide as l-value
                 * 
                 * @return measurement& 
                 * 
                 */
                constexpr measurement& operator/=(const measurement& meas) { 

                    if (meas.value_ == 0.0) throw std::runtime_error("Cannot divide measurement by 0");
                    value_ /= meas.value_;
                    units_ /= meas.units_;                                    
                    return *this; 
                
                }


                /**
                 * @brief Divide this measurement and a measurement
                 * 
                 * @param other: the measurement to divide as r-value
                 * 
                 * @return measurement& 
                 * 
                 */
                constexpr measurement& operator/=(measurement&& meas) { 
                    
                    if (meas.value_ == 0.0) throw std::runtime_error("Cannot divide measurement by 0");
                    value_ /= std::move(meas.value_);
                    units_ /= std::move(meas.units_);                                    
                    return *this; 
                
                }


                /**
                 * @brief Multiply this measurement and a scalar
                 * 
                 * @param value 
                 * 
                 * @return constexpr measurement& 
                 * 
                 */
                constexpr measurement& operator*=(const double& value) noexcept { 

                    value_ *= value;    
                    return *this; 
                
                }


                /**
                 * @brief Divide this measurement and a scalar
                 * 
                 * @param value 
                 * 
                 * @return constexpr measurement& 
                 * 
                 */
                constexpr measurement& operator/=(const double& value) { 
                    
                    if (value == 0.0) throw std::runtime_error("Cannot divide measurement by 0");
                    value_ /= value;                    
                    return *this; 
                
                }


                /**
                 * @brief Sum a measurement to this measurement
                 * 
                 * @param other: measurement to add as l-value const reference
                 * 
                 * @return measurement 
                 * 
                 */
                constexpr measurement operator+(const measurement& other) const { 
                    
                    if (units_.base() != other.units_.base()) throw std::invalid_argument("Cannot sum measurements with different unit bases");
                    return measurement(value_ + other.value_as(units_), units_); 
                
                }

                
                /**
                 * @brief Return the opposite of this measurement
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement operator-() const noexcept { 
                    
                    return measurement(-value_, units_); 
                
                }


                /**
                 * @brief Subtract a measurement to this measurement
                 * 
                 * @param other: measurement to subtract as l-value const reference
                 *  
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement operator-(const measurement& other) const { 

                    if (units_.base() != other.units_.base()) throw std::invalid_argument("Cannot subtract measurements with different unit bases");
                    return measurement(value_ - other.value_as(units_), units_); 
                
                }

                
                /**
                 * @brief Multiply this measurement and a measurement
                 * 
                 * @param other: measurement to multiply as l-value const reference
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement operator*(const measurement& other) const noexcept { 
                    
                    return measurement(value_ * other.value_, units_ * other.units_); 
                
                }


                /**
                 * @brief Divide this measurement and a measurement
                 * 
                 * @param other: measurement to divide as l-value const reference
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement operator/(const measurement& other) const { 
                    
                    if (other.value_ == 0.0) throw std::runtime_error("Cannot divide measurement by 0");
                    return measurement(value_ / other.value_, units_ / other.units_); 
                
                }

                
                /**
                 * @brief Multiply this measurement and a scalar
                 * 
                 * @param value 
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement operator*(const double& value) const noexcept { 
                    
                    return measurement(value_ * value, units_); 
                    
                }

                
                /**
                 * @brief Divide this measurement and a scalar
                 * 
                 * @param value
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement operator/(const double& value) const { 
                    
                    if (value_ == 0.0) throw std::runtime_error("Cannot divide measurement by 0");
                    return measurement(value_ / value, units_); 
                    
                } 


                /**
                 * @brief Perform a product between a double and a measurement
                 * 
                 * @param val: double
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                friend constexpr measurement operator*(const double& val, 
                                                       const measurement& meas) noexcept { 
                                                    
                    return meas * val; 
                    
                }
                

                /**
                 * @brief Perform a division between a double and a measurement
                 * 
                 * @param val: double
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                friend constexpr measurement operator/(const double& val, 
                                                       const measurement& meas) { 

                    if (meas.value_ == 0.0) throw std::runtime_error("Cannot divide by zero");          
                    return measurement(val / meas.value_, meas.units_.inv()); 
                    
                }


                /**
                 * @brief Equality operator
                 * 
                 * @param other: measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator==(const measurement& other) const { 
                    
                    return math::tools::value_equality_check(value_, (units_ == other.units_) ? other.value_ : other.value_as(units_)); 
                    
                }

                
                /**
                 * @brief Inequality operator
                 * 
                 * @param other: measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator!=(const measurement& other) const { 
                    
                    return !math::tools::value_equality_check(value_, (units_ == other.units_) ? other.value_ : other.value_as(units_)); 
                    
                }
                
                
                /**
                 * @brief More than operator
                 * 
                 * @param other: measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator>(const measurement& other) const { 
                    
                    return value_ > other.value_as(units_); 
                    
                }

                
                /**
                 * @brief Less than operator
                 * 
                 * @param other: measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator<(const measurement& other) const { 
                    
                    return value_ < other.value_as(units_); 
                    
                }

                
                /**
                 * @brief More than or equal operator
                 * 
                 * @param other: measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator>=(const measurement& other) const { 
                    
                    return (value_ > other.value_as(units_)) ? true : math::tools::value_equality_check(value_, other.value_as(units_)); 
                    
                }

                
                /**
                 * @brief Less than or equal operator
                 * 
                 * @param other: measurement to compare as l-value const reference 
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator<=(const measurement& other) const { 
                    
                    return (value_ < other.value_as(units_)) ? true : math::tools::value_equality_check(value_, other.value_as(units_)); 
                    
                }       


                /**
                 * @brief Equality operator of value
                 * 
                 * @param val
                 * 
                 * @return bool
                 *  
                 */
                constexpr bool operator==(const double& val) const { 
                    
                    return (value_ == val) ? true : math::tools::compare_round_equals(value_, val); 
                    
                }

                
                /**
                 * @brief Inequality operator of value
                 * 
                 * @param val
                 * 
                 * @return bool
                 *  
                 */
                constexpr bool operator!=(const double& val) const { 
                    
                    return !operator==(val); 
                    
                }
                
               
                /**
                 * @brief More than with value
                 * 
                 * @param val 
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator>(const double& val) const { 
                    
                    return value_ > val; 
                    
                }
               
               
                /**
                 * @brief Less than with value
                 * 
                 * @param val 
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator<(const double& val) const { 
                    
                    return value_ < val; 
                    
                }
                
                
                /**
                 * @brief More than or equal with value
                 * 
                 * @param val 
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator>=(const double& val) const { 
                    
                    return (value_ >= val) ? true : operator==(val); 
                    
                }
                
                
                /**
                 * @brief Less than or equal with value
                 * 
                 * @param val 
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator<=(const double& val) const { 
                    
                    return value_ <= val ? true : operator==(val); 
                    
                }


                /**
                 * @brief Output operator for a measurement
                 * 
                 * @param os: std::ostream&
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return std::ostream&
                 *  
                 */
                friend std::ostream& operator<<(std::ostream& os, const measurement& meas) { 
                    
                    os << meas.value_ << " " << meas.units_; 
                    return os; 
                    
                }


                /**
                 * @brief Output operator for a measurement
                 * 
                 * @param file: std::ofstream&
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return std::ofstream&
                 *  
                 */
                friend std::ofstream& operator<<(std::ofstream& file, const measurement& meas) { 
                    
                    file << meas.value_ << " " << meas.units_; 
                    return file; 
                    
                }


            // =============================================
            // operations
            // ============================================= 

                /**
                 * @brief Get the absolute measurement object
                 * 
                 * @return constexpr measurement
                 * 
                 */
                friend constexpr measurement abs(const measurement& meas) noexcept { 
                    
                    return (meas.value_ < 0.0) ? -meas : meas; 
                
                }


                /**
                 * @brief Invert the measurement
                 * 
                 * @return constexpr measurement 
                 * 
                 * @note Cannot invert a measurement with a zero value
                 * 
                 */
                constexpr measurement inv() const { 
                    
                    if (value_ == 0) throw std::runtime_error("Cannot invert a measurement with a zero value");
                    return measurement(1 / value_, units_.inv()); 
                
                }

                
                /**
                 * @brief Take the power of the measurement
                 * 
                 * @param power 
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement pow(const int& power) const noexcept { 
                    
                    return measurement(std::pow(value_, power), units_.pow(power)); 
                
                }

                
                /**
                 * @brief Take the square of the measurement
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement square() const noexcept { 
                    
                    return measurement(std::pow(value_, 2), units_.square()); 
                
                }

                
                /**
                 * @brief Take the cube of the measurement
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement cube() const noexcept { 
                    
                    return measurement(std::pow(value_, 3), units_.cube()); 
                
                }


                /**
                 * @brief Take the root power of the measurement
                 * 
                 * @param power 
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement root(const int& power) const { 
                    
                    return measurement(std::pow(value_, 1.0 / power), units_.root(power)); 
                
                }

                
                /**
                 * @brief Take the square root of the measurement
                 * 
                 * @return constexpr measurement
                 *  
                 */
                constexpr measurement sqrt() const { 
                    
                    return measurement(std::sqrt(value_), units_.sqrt()); 
                
                }

                
                /**
                 * @brief Take the cubic root of the measurement
                 * 
                 * @return constexpr measurement
                 *  
                 */
                constexpr measurement cbrt() const { 
                    
                    return measurement(std::cbrt(value_), units_.cbrt()); 
                
                }


                /**
                 * @brief Take the sine of a measurement
                 * 
                 * @param meas: measurement 
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement sin(const measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the sine of a measurement that is not in radians"); 
                    else return measurement(std::sin(meas.value_), unitless); 
                
                }


                /**
                 * @brief Take the cosine of a measurement
                 * 
                 * @param meas: measurement 
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement cos(const measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the cosine of a measurement that is not in radians"); 
                    else return measurement(std::cos(meas.value_), unitless); 
                
                }


                /**
                 * @brief Take the tangent of a measurement
                 * 
                 * @param meas: measurement 
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement tan(const measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the tangent of a measurement that is not in radians"); 
                    else return measurement(std::tan(meas.value_), unitless);

                }


                /**
                 * @brief Take the arcsine of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement asin(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the arcsine of a measurement that is not unitless"); 
                    else return measurement(std::asin(meas.value_), rad);

                }


                /**
                 * @brief Take the arccosine of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement acos(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the arccosine of a measurement that is not unitless"); 
                    else return measurement(std::acos(meas.value_), rad);

                }


                /**
                 * @brief Take the arctangent of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement atan(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the arctangent of a measurement that is not unitless"); 
                    else return measurement(std::atan(meas.value_), rad);

                }


                /**
                 * @brief Take the hyperbolic sine of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement sinh(const measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the hyperbolic sine of a measurement that is not in radians"); 
                    else return measurement(std::sinh(meas.value_), unitless);

                }


                /**
                 * @brief Take the hyperbolic cosine of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement cosh(const measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the hyperbolic cosine of a measurement that is not in radians"); 
                    else return measurement(std::cosh(meas.value_), unitless);
                
                }


                /**
                 * @brief Take the hyperbolic tangent of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement tanh(const measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the hyperbolic tangent of a measurement that is not in radians"); 
                    else return measurement(std::tanh(meas.value_), unitless);

                }


                /**
                 * @brief Take the hyperbolic arcsine of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement asinh(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the hyperbolic arcsine of a measurement that is not unitless"); 
                    else return measurement(std::asinh(meas.value_), rad);

                }


                /**
                 * @brief Take the hyperbolic arccosine of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement acosh(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the hyperbolic arccosine of a measurement that is not unitless"); 
                    else return measurement(std::acosh(meas.value_), rad);
                
                }


                /**
                 * @brief Take the hyperbolic arctangent of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement atanh(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the hyperbolic arctangent of a measurement that is not unitless"); 
                    else return measurement(std::atanh(meas.value_), rad);

                }


            // =============================================                                                                                         
            // set and get methods
            // =============================================  

                /**
                 * @brief Get the value of the measurement
                 * 
                 * @return constexpr const double
                 * 
                 */
                constexpr double value() const noexcept { 
                    
                    return value_; 
                
                }                        
                
                
                /**
                 * @brief Get the value of the measurement
                 * 
                 * @return constexpr double& 
                 * 
                 */
                constexpr double& value() noexcept { 
                    
                    return value_; 
                
                }  


                /**
                 * @brief Get the value of the measurement expressed in another units of the measurement
                 * 
                 * @param desired_units 
                 * 
                 * @return constexpr double 
                 * 
                 */
                constexpr double value_as(const unit& desired_units) const { 
                    
                    return (units_ == desired_units) ? value_ : units_.convert(value_, desired_units); 
                
                }


                /**
                 * @brief Get the units of the measurement
                 * 
                 * @return constexpr unit 
                 * 
                 */
                constexpr unit units() const noexcept { 
                    
                    return units_; 
                
                }


                /**
                 * @brief Get the units of the measurement
                 * 
                 * @return constexpr unit&
                 * 
                 */
                constexpr unit& units() noexcept { 
                    
                    return units_; 
                
                }


                /**
                 * @brief Get the measurement object
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement as_measurement() const noexcept {

                    return *this; 

                }


                /**
                 * @brief Get the measurement object
                 * 
                 * @return constexpr measurement&
                 *  
                 */
                constexpr measurement& as_measurement() noexcept {

                    return *this; 

                }
            
                
                /**
                 * @brief Convert the measurement to another units
                 * 
                 * @param desired_units: desired unit of measurement 
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement convert_to(const unit& desired_units) const { 
                    
                    return measurement(units_.convert(value_, desired_units), desired_units); 
                
                }


                /**
                 * @brief Print the measurement
                 * 
                 * @param newline: if set to true it prints a newline character at the end of the measurement
                 * 
                 */
                inline void print(const bool& newline = true) const noexcept { 

                    std::cout << value_ << " " << units_;
                    if (newline) std::cout << "\n"; 

                }

            
            private:

            // =============================================                                                                                         
            // class members
            // =============================================  

                /// @brief The numerical value of the measurement
                double value_;

                /// @brief The units of the measurement
                unit units_;


                /// @brief Friend class to allow access to private members
                friend class uncertain_measurement;
                

        }; // class measurement
        static_assert(sizeof(measurement) == 20, "measurement is not the correct size");


        /// @brief A class for representing a physical quantity with a numerical value, an uncertanty and an unit
        class uncertain_measurement {     
            

            public: 

            // =============================================                                                                                         
            // constructors & destructor 
            // =============================================  

                /// @brief default constructor
                constexpr uncertain_measurement() noexcept :

                    value_{0.0}, 
                    uncertainty_{0.0}, 
                    units_{} {}


                /**
                 * @brief Construct a new uncertain measurement object
                 * 
                 * @param value: the numerical value of the measurement
                 * @param uncertainty: the uncertainty of the measurement
                 * @param units: the units of the measurement
                 * 
                 */
                explicit constexpr uncertain_measurement(const double& val, 
                                                         const double& uncertainty_val, 
                                                         const unit& unit) :

                    value_{val}, 
                    uncertainty_{uncertainty_val}, 
                    units_{unit} {

                        if (uncertainty_ < 0.0) throw std::invalid_argument("Uncertainty cannot be negative");
                        if (uncertainty_ > std::fabs(value_)) throw std::invalid_argument("Uncertainty cannot be greater than the value");

                    }


                /**
                 * @brief Construct a new uncertain measurement object with 0.0 uncertainty
                 * 
                 * @param value: the numerical value of the measurement
                 * @param units: the units of the measurement
                 * 
                 */
                explicit constexpr uncertain_measurement(const double& val, 
                                                         const unit& unit) noexcept :

                    value_{val}, 
                    uncertainty_{0.0},
                    units_(unit) {}


                /**
                 * @brief Construct a new uncertain measurement object from a measurement and an uncertainty
                 * 
                 * @param measurement: measurement (value and unit)
                 * @param uncertainty: uncertainty
                 * 
                 */
                explicit constexpr uncertain_measurement(const measurement& other, 
                                                         const double& uncertainty_val) : 

                    value_{other.value_}, 
                    uncertainty_{uncertainty_val}, 
                    units_{other.units_} {

                        if (uncertainty_ < 0.0) throw std::invalid_argument("Uncertainty cannot be negative");
                        if (uncertainty_ > std::fabs(value_)) throw std::invalid_argument("Uncertainty cannot be greater than the value");

                    }


                /**
                 * @brief Construct a new uncertain measurement object from two measurement 
                 * 
                 * @param value: measurement 
                 * @param uncertainty: measurement 
                 * 
                 * @note the units of the two measurements must have the same base_unit
                 */
                explicit constexpr uncertain_measurement(const measurement& value, 
                                                         const measurement& uncertainty) : 

                    value_{value.value_}, 
                    uncertainty_{uncertainty.value_as(value.units_)}, 
                    units_{value.units_} {

                        if (uncertainty_ < 0.0) 
                            throw std::invalid_argument("Uncertainty cannot be negative");
                        if (uncertainty_ > std::fabs(value_)) 
                            throw std::invalid_argument("Uncertainty cannot be greater than the value");
                        if (value.units_.base() != uncertainty.units_.base()) 
                            throw std::invalid_argument("The units of the two measurements must have the same base_unit");

                    }


                /**
                 * @brief Copy construct a new uncertain measurement object
                 * 
                 * @param other: uncertain measurement to copy as a l-value const reference
                 * 
                 */
                constexpr uncertain_measurement(const uncertain_measurement& other) noexcept :

                    value_{other.value_},
                    uncertainty_{other.uncertainty_},
                    units_{other.units_} {}


                /**
                 * @brief Move construct a new uncertain measurement object
                 * 
                 * @param other: uncertain measurement to move from as a r-value reference
                 * 
                 */
                constexpr uncertain_measurement(uncertain_measurement&& other) noexcept :

                    value_{std::move(other.value_)},
                    uncertainty_{std::move(other.uncertainty_)},
                    units_{std::move(other.units_)} {}


                /// @brief default destructor
                ~uncertain_measurement() = default;
                

            // =============================================                                                                                         
            // operators
            // =============================================  

                /**
                 * @brief Copy assign a uncertain_measurement from another uncertain_measurement
                 * 
                 * @param other: the uncertain_measurement to copy
                 * 
                 * @return uncertai_measurement& 
                 * 
                 */
                constexpr uncertain_measurement& operator=(const uncertain_measurement& other) noexcept {
                   
                    value_ = other.value_;
                    uncertainty_ = other.uncertainty_; 
                    units_ = other.units_;
                    return *this;

                }   


                /**
                 * @brief Move assign a uncertain_measurement from another uncertain_measurement
                 * 
                 * @param other: the uncertain_measurement to move from
                 * 
                 * @return uncertai_measurement& 
                 * 
                 */
                constexpr uncertain_measurement& operator=(uncertain_measurement&& other) noexcept {
                   
                    value_ = std::move(other.value_);
                    uncertainty_ = std::move(other.uncertainty_); 
                    units_ = std::move(other.units_);
                    return *this;
                    
                }   


                /**
                 * @brief Copy assign a uncertain_measurement from a measurement
                 * 
                 * @param other: the uncertain_measurement to copy
                 * 
                 * @return uncertai_measurement& 
                 * 
                 */
                constexpr uncertain_measurement& operator=(const measurement& other) noexcept {
                   
                    value_ = other.value_;
                    uncertainty_ = 0.0; 
                    units_ = other.units_;
                    return *this;

                }   


                /**
                 * @brief Move assign a uncertain_measurement from a measurement
                 * 
                 * @param other: the uncertain_measurement to move from
                 * 
                 * @return uncertai_measurement& 
                 * 
                 */
                constexpr uncertain_measurement& operator=(measurement&& other) noexcept {
                   
                    value_ = std::move(other.value_);
                    uncertainty_ = 0.0; 
                    units_ = std::move(other.units_);
                    return *this;
                    
                }   


                /**
                 * @brief Compute a product and calculate the new uncertainties using the root sum of squares(rss) method
                 * 
                 * @param other: uncertain_measurement to multiply with
                 *  
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement operator*(const uncertain_measurement& other) const noexcept {

                    double tval1 = uncertainty_ / value_;
                    double tval2 = other.uncertainty_ / other.value_;
                    double ntol = std::sqrt(tval1 * tval1 + tval2 * tval2);
                    double nval = value_ * other.value_;
                    return uncertain_measurement(nval, nval * ntol, units_ * other.units_);

                }
        

                /**
                 * @brief Compute a product and calculate the new uncertainties using the simple method for uncertainty propagation
                 * 
                 * @param other: uncertain_measurement to multiply with
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement simple_product(const uncertain_measurement& other) const noexcept {
                    
                    double ntol = uncertainty_ / value_ + other.uncertainty_ / other.value_;
                    double nval = value_ * other.value_;
                    return uncertain_measurement(nval, nval * ntol, units_ * other.units_);
                
                }
                

                /**
                 * @brief Compute a product with a measurement, equivalent to uncertain_measurement multiplication with 0 uncertainty
                 * 
                 * @param other: measurement to multiply with
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement operator*(const measurement& other) const noexcept {
                    
                    return uncertain_measurement(value_ * other.value_, other.value_ * uncertainty_, units_ * other.units_);
                
                }


                /**
                 * @brief Compute a product with a double and calculate the new uncertainties using the simple uncertainty multiplication method
                 * 
                 * @param val: the double to multiply with
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement operator*(const double& val) const noexcept { 
                    
                    return uncertain_measurement(value_ * val, uncertainty_ * val, units_); 
                
                }


                /**
                 * @brief Compute a division and calculate the new uncertainties using the root sum of squares(rss) method
                 * 
                 * @param other: uncertain_measurement to divide by
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement operator/(const uncertain_measurement& other) const {
                    
                    if (other.value_ == 0.0) throw std::invalid_argument("Cannot divide uncertain_measurement by 0");
                    double tval1 = uncertainty_ / value_;
                    double tval2 = other.uncertainty_ / other.value_;
                    double ntol = std::sqrt(tval1 * tval1 + tval2 * tval2);
                    double nval = value_ / other.value_;
                    return uncertain_measurement(nval, nval * ntol, units_ / other.units_);
                
                }


                /**
                 * @brief Compute a division and calculate the new uncertainties using simple method for uncertainty propagation
                 * 
                 * @param other: uncertain_measurement to divide by
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement simple_divide(const uncertain_measurement& other) const {
                    
                    if (other.value_ == 0.0) throw std::invalid_argument("Cannot divide uncertain_measurement by 0");
                    double ntol = uncertainty_ / value_ + other.uncertainty_ / other.value_;
                    double nval = value_ / other.value_;
                    return uncertain_measurement(nval, nval * ntol, units_ / other.units_);
                
                }


                /**
                 * @brief Compute a division with a measurement, equivalent to uncertain_measurement division with 0 uncertainty
                 * 
                 * @param other: measurement to divide by
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement operator/(const measurement& other) const {
                    
                    if (other.value_ == 0.0) throw std::invalid_argument("Cannot divide uncertain_measurement by 0");
                    return uncertain_measurement(value_ / other.value_, uncertainty_ / other.value_, units_ / other.units_);
                
                }


                /**
                 * @brief Divide this uncertain_measurement with a double
                 * 
                 * @param val: the double to divide with
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement operator/(const double& val) const {

                    if (val == 0.0) throw std::invalid_argument("Cannot divide uncertain_measurement by 0");
                    return uncertain_measurement(value_ / val, uncertainty_ / val, units_);
                
                }


                /**
                 * @brief Compute an addition and calculate the new uncertainties using the root sum of squares(rss) method
                 * 
                 * @param other: uncertain_measurement to sum with
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement operator+(const uncertain_measurement& other) const {
                    
                    if (units_.base() != other.units_.base()) throw std::invalid_argument("Cannot add uncertain_measurements with different unit bases");
                    double cval = other.units_.convertion_factor(units_);
                    double ntol = std::sqrt(uncertainty_ * uncertainty_ + cval * cval * other.uncertainty_ * other.uncertainty_);
                    return uncertain_measurement(value_ + cval * other.value_, ntol, units_);
               
                }


                /**
                 * @brief Compute an addition and calculate the new uncertainties using the simple uncertainty summation method
                 * 
                 * @param other: uncertain_measurement to sum with
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement simple_add(const uncertain_measurement& other) const {
                    
                    if (units_.base() != other.units_.base()) throw std::invalid_argument("Cannot add uncertain_measurements with different unit bases");
                    double cval = other.units_.convertion_factor(units_);
                    double ntol = uncertainty_ + other.uncertainty_ * cval;
                    return uncertain_measurement(value_ + cval * other.value_, ntol, units_);
                
                }


                /**
                 * @brief Compute an addition with a measurement and calculate the new uncertainties using the simple uncertainty summation method
                 * 
                 * @param other: measurement to sum with
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement operator+(const measurement& other) const {
                    
                    if (units_.base() != other.units_.base()) throw std::invalid_argument("Cannot add uncertain_measurement and measurement with different unit bases");
                    return uncertain_measurement(value_ + other.value_as(units_), uncertainty_, units_);
                
                }


                /**
                 * @brief Return the opposite of this uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement operator-() const noexcept {

                    return uncertain_measurement(-value_, uncertainty_, units_);
                
                }
                

                /**
                 * @brief Compute a subtraction and calculate the new uncertainties using the root sum of squares(rss) method
                 * 
                 * @param other: uncertain_measurement to subtract with
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement operator-(const uncertain_measurement& other) const {
                    
                    if (units_.base() != other.units_.base()) throw std::invalid_argument("Cannot subtract uncertain_measurements with different unit bases");
                    double cval = other.units_.convertion_factor(units_);
                    double ntol = std::sqrt(uncertainty_ * uncertainty_ + cval * cval * other.uncertainty_ * other.uncertainty_);
                    return uncertain_measurement(value_ - cval * other.value_, ntol, units_);
               
                }


                /**
                 * @brief Compute a subtraction and calculate the new uncertainties using the simple uncertainty summation method
                 * 
                 * @param other: uncertain_measurement to subtract with
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement simple_subtract(const uncertain_measurement& other) const {
                    
                    if (units_.base() != other.units_.base()) throw std::invalid_argument("Cannot subtract uncertain_measurements with different unit bases");
                    auto cval = other.units_.convertion_factor(units_);
                    double ntol = uncertainty_ + other.uncertainty_ * cval;
                    return uncertain_measurement(value_ - cval * other.value_, ntol, units_);
                
                }

                
                /**
                 * @brief Compute a subtraction with a measurement and calculate the new uncertainties using the simple uncertainty summation method
                 * 
                 * @param other: measurement to subtract with
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement operator-(const measurement& other) const {

                    if (units_.base() != other.units_.base()) throw std::invalid_argument("Cannot subtract uncertain_measurement and measurement with different unit bases");
                    return uncertain_measurement(value_ - other.value_as(units_), uncertainty_, units_);
                
                }
  

                /**
                 * @brief Equality operator between uncertain_measurement and measurement
                 * 
                 * @param other: measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator==(const measurement& other) const noexcept {
                    
                    if (uncertainty_ == 0.0) return (value_ == other.value_as(units_)) ? true : math::tools::compare_round_equals(value_, other.value_as(units_)); 
                    return (other.value_as(units_) >= (value_ - uncertainty_) && other.value_as(units_) <= (value_ + uncertainty_));
               
                }


                /**
                 * @brief Equality operator
                 * 
                 * @param other: uncertain_measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator==(const uncertain_measurement& other) const noexcept { 
                    
                    return (simple_subtract(other) == measurement(0.0, units_)); 
                    
                }


                /**
                 * @brief Inequality operator between uncertain_measurement and measurement
                 * 
                 * @param other: measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator!=(const measurement& other) const noexcept { 
                    
                    return !operator==(other); 
                    
                }


                /**
                 * @brief Inequality operator
                 * 
                 * @param other: uncertain_measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator!=(const uncertain_measurement& other) const noexcept { 
                    
                    return !operator==(other); 
                    
                }

                
                /**
                 * @brief More than operator
                 * 
                 * @param other: uncertain_measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator>(const uncertain_measurement& other) const noexcept { 
                    
                    return value_ > other.value_as(units_); 
                    
                }

                
                /**
                 * @brief More than operator
                 * 
                 * @param other: measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator>(const measurement& other) const noexcept { 
                    
                    return value_ > other.value_as(units_); 
                    
                }


                /** 
                 * @brief More than operator
                 * 
                 * @param val: double to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator>(const double& val) const noexcept { 
                    
                    return value_ > val; 
                    
                }

                
                /**
                 * @brief Less than operator
                 * 
                 * @param other: uncertain_measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator<(const uncertain_measurement& other) const noexcept { 
                    
                    return value_ < other.value_as(units_); 
                    
                }
               
                
                /**
                 * @brief Less than operator
                 * 
                 * @param other: measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator<(const measurement& other) const noexcept { 
                    
                    return value_ < other.value_as(units_); 
                    
                }

                
                /**
                 * @brief Less than operator
                 * 
                 * @param val: double to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator<(const double& val) const noexcept { 
                    
                    return value_ < val; 
                    
                }

                
                /**
                 * @brief More than or equal operator
                 * 
                 * @param other: uncertain_measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator>=(const uncertain_measurement& other) const noexcept {
                    
                    return (simple_subtract(other).value_ >= 0.0) ? true : (simple_subtract(other) == measurement(0.0, units_));
                
                }


                /**
                 * @brief More than or equal operator
                 * 
                 * @param other: measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator>=(const measurement& other) const noexcept {
                    
                    return (value_ >= other.value_as(units_)) ? true : operator==(measurement(other.value_as(units_), units_));
                
                }


                /**
                 * @brief More than or equal operator
                 * 
                 * @param val: double to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator>=(const double& val) const noexcept { 

                    return value_ >= val; 
                    
                }


                /**
                 * @brief Less than or equal operator
                 * 
                 * @param other: uncertain_measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator<=(const uncertain_measurement& other) const noexcept {
                    
                    return (simple_subtract(other).value_ <= 0.0) ? true : (simple_subtract(other) == measurement(0.0, units_));
               
                }


                /**
                 * @brief Less than or equal operator
                 * 
                 * @param other: measurement to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator<=(const measurement& other) const noexcept {
                    
                    return (value_ <= other.value_as(units_)) ? true : operator==(measurement(other.value_as(units_), units_));
               
                }


                /**
                 * @brief Less than or equal operator
                 * 
                 * @param val: double to compare as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator<=(const double& val) const noexcept { 

                    return value_ <= val; 
                    
                }

                
                /**
                 * @brief Perform a product between a measurement and an uncertain_measurement
                 * 
                 * @param meas: measurement
                 * @param umeas: uncertain_measurement
                 *  
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                friend constexpr uncertain_measurement operator*(const measurement& meas, 
                                                                 const uncertain_measurement& umeas) noexcept { 
                                    
                    return umeas.operator*(meas); 
                                
                }


                /**
                 * @brief Perform a product between a double and an uncertain_measurement
                 * 
                 * @param value: double
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                friend constexpr uncertain_measurement operator*(const double& value, 
                                                                 const uncertain_measurement& umeas) noexcept { 
                                    
                    return umeas.operator*(value); 
                                
                }


                /**
                 * @brief Perform a division between a measurement and an uncertain_measurement
                 * 
                 * @param meas: measurement
                 * @param umeas: uncertain_measurement
                 *  
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                friend constexpr uncertain_measurement operator/(const measurement& meas, 
                                                                 const uncertain_measurement& umeas) {
                                                            
                    if (umeas.value() == 0.0) throw std::runtime_error("Cannot divide by zero");
                    double ntol = umeas.uncertainty() / umeas.value();
                    double nval = meas.value() / umeas.value();
                    return uncertain_measurement(nval, nval * ntol, meas.units() / umeas.units());
                
                }


                /**
                 * @brief Perform a division between a double and an uncertain_measurement
                 * 
                 * @param value: double
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                friend constexpr uncertain_measurement operator/(const double& v1, 
                                                                 const uncertain_measurement& umeas) {

                    if (umeas.value() == 0.0) throw std::runtime_error("Cannot divide by zero");
                    double ntol = umeas.uncertainty() / umeas.value();
                    double nval = v1 / umeas.value();
                    return uncertain_measurement(nval, nval * ntol, umeas.units().inv());
                
                }

                
                /**
                 * @brief Perform a sum between a measurement and an uncertain_measurement
                 * 
                 * @param meas: measurement 
                 * @param umeas: uncertain_measurement
                 *  
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                friend constexpr uncertain_measurement operator+(const measurement& meas, 
                                                                 const uncertain_measurement& umeas) {
                    
                    if (meas.units().base() != umeas.units().base()) throw std::invalid_argument("Cannot sum measurement and uncertain_measurement with different unit bases");
                    double cval = umeas.units().convertion_factor(meas.units());
                    double ntol = umeas.uncertainty() * cval;
                    return uncertain_measurement(meas.value() + cval * umeas.value(), ntol, meas.units());
                
                }


                /**
                 * @brief Perform a subtraction between a measurement and an uncertain_measurement
                 * 
                 * @param meas: measurement 
                 * @param umeas: uncertain_measurement
                 *  
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                friend constexpr uncertain_measurement operator-(const measurement& meas, 
                                                                 const uncertain_measurement& umeas) {

                    if (meas.units().base() != umeas.units().base()) throw std::invalid_argument("Cannot subtract measurement and uncertain_measurement with different unit bases");
                    double cval = umeas.units().convertion_factor(meas.units());
                    double ntol = umeas.uncertainty() * cval;
                    return uncertain_measurement(meas.value() - cval * umeas.value(), ntol, meas.units());
                
                }


                /**
                 * @brief Equality operator between a measurement and an uncertain_measurement
                 * 
                 * @param meas: measurement
                 * @param umeas: uncertain_measurement
                 * 
                 * @return bool
                 * 
                 */
                friend constexpr bool operator==(const measurement& meas, 
                                                 const uncertain_measurement& umeas) noexcept { 
                    
                    return umeas == meas; 
                
                }


                /**
                 * @brief Inequality operator between a measurement and an uncertain_measurement
                 * 
                 * @param meas: measurement
                 * @param umeas: uncertain_measurement
                 * 
                 * @return bool
                 * 
                 */
                friend constexpr bool operator!=(const measurement& meas, 
                                                 const uncertain_measurement& umeas) noexcept { 
                    
                    return umeas != meas; 
                
                }
                

                /**
                 * @brief More operator between a measurement and an uncertain_measurement
                 * 
                 * @param meas: measurement
                 * @param umeas: uncertain_measurement
                 * 
                 * @return bool
                 * 
                 */
                friend constexpr bool operator>(const measurement& meas, 
                                                const uncertain_measurement& umeas) noexcept { 
                    
                    return meas.value() > umeas.value(); 
                
                }
                

                /**
                 * @brief Less operator between a measurement and an uncertain_measurement
                 * 
                 * @param meas: measurement
                 * @param umeas: uncertain_measurement
                 * 
                 * @return bool
                 *  
                 */
                friend constexpr bool operator<(const measurement& meas, 
                                                const uncertain_measurement& umeas) noexcept { 
                    
                    return meas.value() < umeas.value(); 
                
                }
                

                /**
                 * @brief More or equal operator between a measurement and an uncertain_measurement
                 * 
                 * @param meas: measurement
                 * @param umeas: uncertain_measurement
                 * 
                 * @return bool
                 * 
                 */
                friend constexpr bool operator>=(const measurement& meas, 
                                                 const uncertain_measurement& umeas) noexcept { 
                    
                    return (meas > umeas) ? true : (umeas == meas); 
                
                }
                

                /**
                * @brief Less or equal operator between a measurement and an uncertain_measurement
                * 
                * @param meas: measurement
                * @param umeas: uncertain_measurement
                * 
                * @return bool
                * 
                */
                friend constexpr bool operator<=(const measurement& meas, 
                                                 const uncertain_measurement& umeas) noexcept { 
                    
                    return (meas < umeas) ? true : (umeas == meas); 
                
                }
        

                /**
                 * @brief Output operator for a uncertain_measurement
                 * 
                 * @param os: std::ostream&
                 * @param meas: uncertain_measurement as l-value const reference
                 * 
                 * @note if the precision of the uncertain_measurement is 0, the uncertainty is not printed
                 * @note scientific notation is used if the value is greater than 1e4 or less than 1e-4
                 * 
                 * @return std::ostream&
                 */
                friend std::ostream& operator<<(std::ostream& os, const uncertain_measurement& umeas) noexcept { 

                    double abs_value = std::fabs(umeas.value_);
                    
                    // first significative digit positions
                    int32_t n_val = ((umeas.uncertainty_ >= 1) ? 
                                        std::ceil(std::log10(abs_value)) : 
                                        ((abs_value >= 1) ? 
                                            std::ceil(std::log10(abs_value)) : 
                                            std::floor(std::log10(abs_value)))); 

                    int32_t n_unc = ((umeas.uncertainty_ >= 1) ? 
                                        std::ceil(std::log10(umeas.uncertainty_)) : 
                                        std::floor(std::log10(umeas.uncertainty_))); 

                    // check if the uncertainty needs to be printed
                    if (umeas.uncertainty_ == 0.0) os << umeas.as_measurement(); 
                    else if (umeas.uncertainty_ >= 1) { 
                        
                        // check for scientific notation
                        if (abs_value > 1.e4 || umeas.uncertainty_ > 1.e4) {

                            os << std::scientific << std::setprecision(1 + n_val - n_unc) << umeas.value_; 

                        } else {

                            os << std::setprecision(2 + n_val - n_unc) << umeas.value_; 
                            os << std::fixed;

                        } 

                        os << std::setprecision(0) << " ± " << umeas.uncertainty_;
                        
                    } else {
                        
                        // check for scientific notation
                        if (abs_value < 1.e-4 || umeas.uncertainty_ < 1.e-4 || abs_value > 1.e4) {

                            os << std::scientific << std::setprecision(1 + std::fabs(n_val) + n_unc) << umeas.value_; 
                            os << std::setprecision(0); 
                            
                        } else os << std::fixed << std::setprecision(std::fabs(n_unc)) << umeas.value_; 

                        os << " ± " << umeas.uncertainty_;

                    }

                    // printing the units 
                    os << " " << umeas.units_; 

                    os << std::defaultfloat << std::setprecision(6);
                    return os; 
                    
                }


                /**
                 * @brief Output operator for a uncertain_measurement
                 * 
                 * @param file: std::ofstream&
                 * @param meas: uncertain_measurement as l-value const reference
                 * 
                 * @return std::ofstream&
                 *  
                 */
                friend std::ofstream& operator<<(std::ofstream& file, const uncertain_measurement& umeas) noexcept { 

                    double abs_value = std::fabs(umeas.value_);
                    
                    // first significative digit positions
                    int32_t n_val = ((umeas.uncertainty_ >= 1) ? 
                                        std::ceil(std::log10(abs_value)) : 
                                        ((abs_value >= 1) ? 
                                            std::ceil(std::log10(abs_value)) : 
                                            std::floor(std::log10(abs_value)))); 

                    int32_t n_unc = ((umeas.uncertainty_ >= 1) ? 
                                        std::ceil(std::log10(umeas.uncertainty_)) : 
                                        std::floor(std::log10(umeas.uncertainty_))); 

                    // check if the uncertainty needs to be printed
                    if (umeas.uncertainty_ == 0.0) file << umeas.as_measurement(); 
                    else if (umeas.uncertainty_ >= 1) { 
                        
                        // check for scientific notation
                        if (abs_value > 1.e4 || umeas.uncertainty_ > 1.e4) {

                            file << std::scientific << std::setprecision(1 + n_val - n_unc) << umeas.value_; 

                        } else {

                            file << std::setprecision(2 + n_val - n_unc) << umeas.value_; 
                            file << std::fixed;

                        } 

                        file << std::setprecision(0) << " ± " << umeas.uncertainty_;
                        
                    } else {
                        
                        // check for scientific notation
                        if (abs_value < 1.e-4 || umeas.uncertainty_ < 1.e-4 || abs_value > 1.e4) {

                            file << std::scientific << std::setprecision(1 + std::fabs(n_val) + n_unc) << umeas.value_; 
                            file << std::setprecision(0); 
                            
                        } else file << std::fixed << std::setprecision(std::fabs(n_unc)) << umeas.value_; 

                        file << " ± " << umeas.uncertainty_;

                    }

                    // printing the units 
                    file << " " << umeas.units_; 

                    file << std::defaultfloat;
                    return file; 
                    
                    
                }


            // =============================================
            // operations
            // ============================================= 

                /**                  
                 * @brief Invert the uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 * @note Cannot invert an uncertain_measurement with a zero value
                 * @note The uncertainty is not inverted
                 * 
                 */
                constexpr uncertain_measurement inv() const { 
                    
                    if (value_ == 0) throw std::runtime_error("Cannot invert an uncertain_measurement with a zero value");
                    return uncertain_measurement(1 / value_, uncertainty_ / std::pow(value_, 2), units_.inv());
                    
                } 


                /**
                 * @brief Take the power of the uncertain_measurement
                 * 
                 * @param power 
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement pow(const int& power) const noexcept { 
                    
                    return uncertain_measurement(std::pow(value_, power), power * std::pow(value_, power - 1) * uncertainty_, units_.pow(power)); 
                    
                }

                
                /**
                 * @brief Take the square of the uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement square() const noexcept { 
                    
                    return uncertain_measurement(std::pow(value_, 2), 2 * value_ * uncertainty_, units_.square()); 
                    
                }

                
                /**
                 * @brief Take the cube of the uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement cube() const noexcept { 
                    
                    return uncertain_measurement(std::pow(value_, 3), 3 * std::pow(value_, 2) * uncertainty_, units_.cube()); 
                    
                }

                
                /**
                 * @brief Take the root power of the uncertain_measurement
                 * 
                 * @param power 
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement root(const int& power) const { 
                    
                    return uncertain_measurement(std::pow(value_, 1.0 / power), std::pow(value_, 1.0 / power - 1) * uncertainty_ / power, units_.root(power)); 
                    
                }

                
                /**
                 * @brief Take the square root of the uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement
                 *  
                 */
                constexpr uncertain_measurement sqrt() const { 
                    
                    return uncertain_measurement(std::sqrt(value_), uncertainty_ / (2 * std::sqrt(value_)), units_.sqrt()); 
                    
                }


                /**
                 * @brief Take the cubic root of the uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement
                 *  
                 */                
                constexpr uncertain_measurement cbrt() const { 
                    
                    return uncertain_measurement(std::cbrt(value_), std::pow(value_, - 2. / 3.) * uncertainty_ / 3, units_.cbrt());
                    
                }


                /**
                 * @brief Take the sine of a measurement
                 * 
                 * @param meas: measurement 
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement sin(const uncertain_measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the sine of a measurement that is not in radians"); 
                    else return uncertain_measurement(std::sin(meas.value_), std::fabs(std::cos(meas.value_)) * meas.uncertainty_, unitless); 
                
                }


                /**
                 * @brief Take the cosine of a measurement
                 * 
                 * @param meas: measurement 
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement cos(const uncertain_measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the cosine of a measurement that is not in radians"); 
                    else return uncertain_measurement(std::cos(meas.value_), -std::sin(meas.uncertainty_), unitless); 
                
                }


                /**
                 * @brief Take the tangent of a measurement
                 * 
                 * @param meas: measurement 
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement tan(const uncertain_measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the tangent of a measurement that is not in radians"); 
                    else return uncertain_measurement(std::tan(meas.value_), 1 + std::pow(meas.uncertainty_, 2), unitless);

                }


                /**
                 * @brief Take the arcsine of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement asin(const uncertain_measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the arcsine of a measurement that is not unitless"); 
                    else return uncertain_measurement(std::asin(meas.value_), 1. / std::sqrt(1 - std::pow(meas.uncertainty_, 2)), rad);

                }


                /**
                 * @brief Take the arccosine of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement acos(const uncertain_measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the arccosine of a measurement that is not unitless"); 
                    else return uncertain_measurement(std::acos(meas.value_), - 1. / std::sqrt(1 - std::pow(meas.uncertainty_, 2)), rad);

                }


                /**
                 * @brief Take the arctangent of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement atan(const uncertain_measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the arctangent of a measurement that is not unitless"); 
                    else return uncertain_measurement(std::atan(meas.value_), 1. / (1 + std::pow(meas.uncertainty_, 2)), rad);

                }


                /**
                 * @brief Take the hyperbolic sine of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement sinh(const uncertain_measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the hyperbolic sine of a measurement that is not in radians"); 
                    else return uncertain_measurement(std::sinh(meas.value_), std::cosh(meas.uncertainty_), unitless);

                }


                /**
                 * @brief Take the hyperbolic cosine of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement cosh(const uncertain_measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the hyperbolic cosine of a measurement that is not in radians"); 
                    else return uncertain_measurement(std::cosh(meas.value_), std::sinh(meas.uncertainty_), unitless);
                
                }


                /**
                 * @brief Take the hyperbolic tangent of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement tanh(const uncertain_measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the hyperbolic tangent of a measurement that is not in radians"); 
                    else return uncertain_measurement(std::tanh(meas.value_), 1 - std::pow(meas.uncertainty_, 2), unitless);

                }


                /**
                 * @brief Take the hyperbolic arcsine of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement asinh(const uncertain_measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the hyperbolic arcsine of a measurement that is not unitless"); 
                    else return uncertain_measurement(std::asinh(meas.value_), 1. / std::sqrt(std::pow(meas.uncertainty_, 2) + 1), rad);

                }


                /**
                 * @brief Take the hyperbolic arccosine of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement acosh(const uncertain_measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the hyperbolic arccosine of a measurement that is not unitless"); 
                    else return uncertain_measurement(std::acosh(meas.value_), 1. / std::sqrt(std::pow(meas.uncertainty_, 2) - 1), rad);
                
                }


                /**
                 * @brief Take the hyperbolic arctangent of a measurement
                 * 
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement atanh(const uncertain_measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the hyperbolic arctangent of a measurement that is not unitless"); 
                    else return uncertain_measurement(std::atanh(meas.value_), 1. / std::sqrt(1 - std::pow(meas.uncertainty_, 2)), rad);

                }


            // =============================================                                                                                         
            // set & get methods
            // =============================================  
                
                /**
                 * @brief Get the value of the measurement
                 * 
                 * @return constexpr const double
                 * 
                 */
                constexpr double value() const noexcept { 
                    
                    return value_; 
                
                }                        
                
                
                /**
                 * @brief Get the value of the measurement
                 * 
                 * @return constexpr double& 
                 * 
                 */
                constexpr double& value() noexcept { 
                    
                    return value_; 
                
                }  


                /**
                 * @brief Get the value of the measurement expressed in another units of the measurement
                 * 
                 * @param desired_units 
                 * @return constexpr double 
                 * 
                 */
                constexpr double value_as(const unit& desired_units) const { 
                    
                    return (units_ == desired_units) ? value_ : units_.convert(value_, desired_units); 
                
                }


                /**
                 * @brief Cast to a measurement
                 * 
                 * @return measurement 
                 * 
                 */
                constexpr measurement as_measurement() const noexcept { 
                    
                    return measurement(value_, units_); 
                    
                }


                /**
                 * @brief Get the uncertainty of the measurement
                 * 
                 * @return constexpr const double
                 * 
                 */
                constexpr double uncertainty() const noexcept { 
                    
                    return uncertainty_; 
                
                }                        
                
                
                /**
                 * @brief Get the uncertainty of the measurement
                 * 
                 * @return constexpr double& 
                 * 
                 */
                constexpr double& uncertainty() noexcept { 
                    
                    return uncertainty_; 
                
                }  


                /**
                 * @brief Get the uncertainty of the measurement expressed in another units of the measurement
                 * 
                 * @param desired_units 
                 * @return constexpr double 
                 * 
                 */
                constexpr double uncertainty_as(const unit& desired_units) const { 
                    
                    return (units_ == desired_units) ? uncertainty_ : units_.convert(uncertainty_, desired_units); 
                
                }


                /**
                 * @brief Get the relative uncertainty of the measurement
                 * 
                 * @return constexpr double
                 */
                constexpr double relative_uncertainty() const noexcept { 
                    
                    return uncertainty_ / value_;

                }


                /**
                 * @brief Get the weight of the measurement
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement weight() const {
                        
                    return this->uncertainty_as_measurement().inv().square();
    
                }


                /**
                 * @brief Get the uncertainty as a separate measurement
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement uncertainty_as_measurement() const noexcept { 
                    
                    return measurement(uncertainty_, units_); 
                
                }


                /**
                 * @brief Get the units of the measurement
                 * 
                 * @return constexpr unit 
                 * 
                 */
                constexpr unit units() const noexcept { 
                    
                    return units_; 
                
                }


                /**
                 * @brief Get the units of the measurement
                 * 
                 * @return constexpr unit&
                 * 
                 */
                constexpr unit& units() noexcept { 
                    
                    return units_; 
                
                }


                /**
                 * @brief Get the uncertain_measurement object
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement as_uncertain_measurement() const noexcept {

                    return *this; 

                }


                /**
                 * @brief Get the uncertain_measurement object
                 * 
                 * @return constexpr uncertain_measurement&
                 *  
                 */
                constexpr uncertain_measurement& as_uncertain_measurement() noexcept {

                    return *this; 

                }
                

                /**
                 * @brief Convert the uncertain_measurement to another units
                 * 
                 * @param desired_units: desired unit of measurement 
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 */
                constexpr uncertain_measurement convert_to(const unit& newUnits) const noexcept {

                    auto cval = units_.convertion_factor(newUnits);
                    return uncertain_measurement(cval * value_, uncertainty_ * cval, newUnits);

                }


                /**
                 * @brief Print the measurement to a ostream
                 * 
                 * @param newLine: bool
                 * 
                 * @return void
                 * 
                 */
                void print(const bool& newline = true) const noexcept { 

                    std::cout << *this; 
                    if (newline) std::cout << "\n";

                }


            private:

            // =============================================                                                                                         
            // class members
            // =============================================  

                /// @brief the numerical value of the measurement
                double value_; 
                
                /// @brief the uncertainty of the measurement
                double uncertainty_;

                /// @brief the units of the measurement
                unit units_;


        }; // class uncertain_measurement
        static_assert(sizeof(uncertain_measurement) == 28, "uncertain_measurement is not the correct size");


        // =============================================
        // some more operators
        // =============================================

        /**
         * @brief Create a measurement by multiplying a double with an unit
         * 
         * @param val: double
         * @param units: unit
         * 
         * @return constexpr measurement 
         * 
         */
        constexpr measurement operator*(const double& val, 
                                        const unit& units) noexcept { 
                                            
            return measurement(val, units); 
            
        }


        /**
         * @brief Create a measurement by dividing a double by an unit
         * 
         * @param val: double
         * @param units: unit
         * 
         * @return constexpr measurement 
         * 
         */
        constexpr measurement operator/(const double& val, 
                                        const unit& units) noexcept { 
                                            
            return measurement(val, units.inv()); 
            
        }


        /**
         * @brief Read measurements from a ifstream
         * 
         * @param file_path_name: std::string
         * @param units: unit
         * @param size: uint32_t (default = 0)
         * 
         * @return std::vector<measurement> 
         */
        std::vector<measurement> read_measurements(const std::string& file_path_name, const unit& units, const uint32_t& size = 0) {

            double value{}; 
            std::vector<measurement> data;
            if (size > 0) data.reserve(size);
            std::ifstream file(file_path_name);
            if (file.is_open()) {

                std::string line;
                while (std::getline(file, line)) {

                    std::stringstream ss(line);
                    ss >> value;
                    data.emplace_back(measurement(value, units));

                }

            } else throw std::runtime_error("Could not open file: " + file_path_name);

            return data;

        }


        /**
         * @brief Read uncertain_measurements from a ifstream
         * 
         * @param file_path_name: std::string
         * @param units: unit
         * @param size: uint32_t (default = 0)
         * 
         * @return std::vector<uncertain_measurements> 
         */
        std::vector<uncertain_measurement> read_uncertain_measurements(const std::string& file_path_name, const unit& units, const uint32_t& size = 0) {

            double value{}, uncertainty{};
            std::vector<uncertain_measurement> data;
            if (size > 0) data.reserve(size);
            std::ifstream file(file_path_name);
            if (file.is_open()) {

                std::string line;
                while (std::getline(file, line)) {

                    std::stringstream ss(line);
                    ss >> value >> uncertainty;
                    data.emplace_back(uncertain_measurement(value, uncertainty, units));

                }

            } else throw std::runtime_error("Could not open file: " + file_path_name);

            return data;

        }
        
        
    } // namespace measurements


    /// @brief Namespace contains some usefull class and containers for representing physical quantities
    namespace tools {
  

        using namespace physics::measurements; 


        /**
         * @brief Class expressing a generic vector of measurements in a n-dimentional system
         * 
         * @tparam DIM: the number of dimensions
         * 
         */
        template <size_t DIM> 
        class vector {


            using scalar = double; 


            public: 

            // =============================================
            // constructors and destructor
            // =============================================

                /**
                 * @brief Construct a new vector object
                 * 
                 * @param unit: unit of measurement (default unit is unitless)
                 * 
                 */
                explicit constexpr vector(const unit_base& base = base::default_type) noexcept {

                    data_.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) data_.emplace_back(0.0 * unit(base));

                }


                /**
                 * @brief Construct a new vector object from an std::vector of measurements
                 * 
                 * @param data: std::vector<measurement> as l-value const reference
                 * 
                 */
                constexpr vector(const std::vector<measurement>& data) {

                    if (data.size() != DIM) throw std::invalid_argument("Cannot construct a vector with an std::vector of measurements with different size");
                    data_.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) data_.emplace_back(data[i]); 

                }


                /**
                 * @brief Construct a new vector object from an std::vector of measurements
                 * 
                 * @param data: std::vector<measurement> as r-value reference
                 * 
                 */
                constexpr vector(std::vector<measurement>&& data) {

                    if (data.size() != DIM) throw std::invalid_argument("Cannot construct a vector with an std::vector of measurements with different size");
                    data_.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) data_.emplace_back(std::move(data[i]));

                }


                /**
                 * @brief Copy construct a new vector object from a vector
                 * 
                 * @param other: vector as l-value const reference
                 * 
                 */
                constexpr vector(const vector& other) noexcept {
                
                    data_.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) data_.emplace_back(other.data_[i]); 

                }


                /**
                 * @brief Move construct a new vector object from a vector
                 * 
                 * @param other: vector as r-value reference
                 * 
                 */
                constexpr vector(vector&& other) noexcept {

                    data_.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) data_.emplace_back(std::move(other.data_[i])); 

                }


                /// @brief Default destructor
                ~vector() = default;


            // =============================================
            // operators
            // =============================================

                /**
                 * @brief Copy assignment operator
                 * 
                 * @param other: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 * 
                 */
                constexpr vector& operator=(const vector& other) noexcept {

                    data_.clear(); 
                    for (size_t i{}; i < DIM; ++i) data_[i] = other.data_[i]; 
                    return *this; 

                }


                /**
                 * @brief Move assignment operator
                 * 
                 * @param other: vector as r-value reference
                 * 
                 * @return constexpr vector& 
                 * 
                 */
                constexpr vector& operator=(vector&& other) noexcept {
                    
                    data_.clear();
                    for (size_t i{}; i < DIM; ++i) data_[i] = std::move(other.data_[i]); 
                    return *this; 

                }


                /**
                 * @brief Add a vector to the current vector
                 * 
                 * @param other: vector to add as l-value const reference
                 * 
                 * @return constexpr vector& 
                 * 
                 * @note: the two vectors must have the same unit of measurement and the same size
                 * 
                 */
                constexpr vector& operator+=(const vector& other) noexcept {

                    for (size_t i{}; i < DIM; ++i) data_[i] += other.data_[i];
                    return *this; 

                }


                /**
                 * @brief Add a vector to the current vector
                 * 
                 * @param other: vector to add as r-value reference
                 * 
                 * @return constexpr vector& 
                 * 
                 * @note: the two vectors must have the same unit of measurement and the same size
                 * 
                 */
                constexpr vector& operator+=(vector&& other) noexcept {

                    for (size_t i{}; i < DIM; ++i) data_[i] += std::move(other.data_[i]);
                    return *this; 

                }

                
                /**
                 * @brief Subtract a vector to the current vector
                 * 
                 * @param other: vector to subtract as l-value const reference
                 * 
                 * @return constexpr vector& 
                 * 
                 * @note: the two vectors must have the same unit of measurement and the same size
                 * 
                 */
                constexpr vector& operator-=(const vector& other) noexcept {

                    for (size_t i{}; i < DIM; ++i) data_[i] -= other.data_[i];
                    return *this; 

                }


                /**
                 * @brief Subtract a vector to the current vector
                 * 
                 * @param other: vector to subtract as r-value reference
                 * 
                 * @return constexpr vector& 
                 * 
                 * @note: the two vectors must have the same unit of measurement and the same size
                 * 
                 */
                constexpr vector& operator-=(vector&& other) noexcept {

                    for (size_t i{}; i < DIM; ++i) data_[i] -= std::move(other.data_[i]);
                    return *this; 

                }


                /**
                 * @brief Multiply the current vector by a measurement
                 * 
                 * @param meas: measurement to multiply with as l-value const reference
                 * 
                 * @return constexpr vector& 
                 * 
                 */
                constexpr vector& operator*=(const measurement& meas) noexcept {

                    for (size_t i{}; i < DIM; ++i) data_[i] *= meas;
                    return *this; 

                }


                /**
                 * @brief Multiply the current vector by a measurement
                 * 
                 * @param meas: measurement to multiply with as r-value reference
                 * 
                 * @return constexpr vector&
                 * 
                 */
                constexpr vector& operator*=(measurement&& meas) noexcept {

                    for (size_t i{}; i < DIM; ++i) data_[i] *= std::move(meas);
                    return *this; 

                }


                /**
                 * @brief Divide the current vector by a measurement
                 * 
                 * @param meas: measurement to divide by as l-value const reference
                 *  
                 * @return constexpr vector& 
                 * 
                 */
                constexpr vector& operator/=(const measurement& meas) {
                    
                    if (meas.value() == 0.0) throw std::runtime_error("Cannot divide by zero");
                    else {
                        for (size_t i{}; i < DIM; ++i) data_[i] /= meas;
                        return *this; 
                    }

                }


                /**
                 * @brief Divide the current vector by a measurement
                 * 
                 * @param meas: measurement to divide by as r-value reference
                 * 
                 * @return constexpr vector&
                 * 
                 */
                constexpr vector& operator/=(measurement&& meas) noexcept {

                    for (size_t i{}; i < DIM; ++i) {
                        if (meas.value() == 0.0) throw std::runtime_error("Cannot divide by zero");
                        else data_[i] /= meas;
                    }
                    return *this; 

                }


                /**
                 * @brief Multiply the current vector by a scalar
                 * 
                 * @param scalar: double as l-value const reference
                 * 
                 * @return constexpr vector& 
                 * 
                 */
                constexpr vector& operator*=(const scalar& scalar) noexcept {

                    for (size_t i{}; i < DIM; ++i) data_[i] *= scalar;
                    return *this; 

                }


                /**
                 * @brief Divide the current vector by a scalar
                 * 
                 * @param scalar: double as l-value const reference
                 *  
                 * @return constexpr vector& 
                 * 
                 */
                constexpr vector& operator/=(const scalar& scalar) {

                    if (scalar == 0.0) throw std::runtime_error("Cannot divide by zero");
                    else {
                        for (size_t i{}; i < DIM; ++i) data_[i] /= scalar;
                        return *this; 
                    }
                    
                }


                /**
                 * @brief Sum the current vector and another vector
                 * 
                 * @param other: vector to add as l-value const reference
                 * 
                 * @return constexpr vector 
                 * 
                 * @note: the two vectors must have the same unit of measurement and the same size
                 * 
                 */
                constexpr vector operator+(const vector& other) const noexcept {

                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i] + other.data_[i]); 
                    return result;

                }


                /**
                 * @brief Subtract the current vector and another vector
                 * 
                 * @param other: vector to subtract as l-value const reference
                 * 
                 * @return constexpr vector 
                 * 
                 * @note: the two vectors must have the same unit of measurement and the same size
                 * 
                 */
                constexpr vector operator-(const vector& other) const noexcept {

                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i] - other.data_[i]); 
                    return result;

                }


                /**
                 * @brief Return the opposite of the current vector
                 * 
                 * @return constexpr vector 
                 * 
                 */
                constexpr vector operator-() const noexcept {

                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(-data_[i]); 
                    return result;

                }


                /**
                 * @brief Multiply the current vector by a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr vector 
                 * 
                 */
                constexpr vector operator*(const measurement& meas) const noexcept {

                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i] * meas); 
                    return result;

                }


                /**
                 * @brief Divide the current vector by a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr vector 
                 * 
                 */
                constexpr vector operator/(const measurement& meas) const noexcept {

                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i] / meas); 
                    return result;

                }


                /**
                 * @brief Multiply the a measurement and a vector
                 * 
                 * @param meas: measurement as l-value const reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 * 
                 */
                friend constexpr vector operator*(const measurement& meas, const vector& vec) noexcept {

                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(meas * vec.data_[i]); 
                    return result;

                }


                /**
                 * @brief Divide the a measurement and a vector
                 * 
                 * @param meas: measurement as l-value const reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 * 
                 */
                friend constexpr vector operator/(const measurement& meas, const vector& other) noexcept {

                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(meas / other.data_[i]); 
                    return result;

                }


                /**
                 * @brief Multiply the current vector by a scalar
                 * 
                 * @param scalar: scalar as l-value const reference
                 * 
                 * @return constexpr vector
                 *  
                 */
                constexpr vector operator*(const scalar& scalar) const noexcept {

                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i] * scalar); 
                    return result;

                }


                /**
                 * @brief Divide the current vector by a scalar
                 * 
                 * @param scalar: scalar as l-value const reference
                 * 
                 * @return constexpr vector
                 *  
                 */
                constexpr vector operator/(const scalar& scalar) const noexcept {

                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i] / scalar); 
                    return result;

                }


                /**
                 * @brief Multiply the a scalar and a vector
                 * 
                 * @param scalar: scalar as l-value const reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 * 
                 */
                friend constexpr vector operator*(const scalar& scalar, const vector& vec) noexcept {

                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(scalar * vec.data_[i] ); 
                    return result;

                }


                /**
                 * @brief Divide the a scalar and a vector
                 * 
                 * @param scalar: scalar as l-value const reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 * 
                 */                
                friend constexpr vector operator/(const scalar& scalar, const vector& vec) noexcept {

                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(scalar / vec.data_[i] ); 
                    return result;

                }


                constexpr vector operator*(const std::vector<scalar>& scalar_vec) const noexcept {

                    assert(scalar_vec.size() == DIM);
                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i] * scalar_vec[i]); 
                    return result;

                }


                constexpr vector operator/(const std::vector<scalar>& scalar_vec) const noexcept {

                    assert(scalar_vec.size() == DIM);
                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i] / scalar_vec[i]); 
                    return result;

                }


                friend constexpr vector operator*(const std::vector<scalar>& scalar_vec, const vector& other) noexcept {

                    assert(scalar_vec.size() == DIM);
                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(scalar_vec[i] * other.data_[i]); 
                    return result;

                }


                friend constexpr vector operator/(const std::vector<scalar>& scalar_vec, const vector& other) noexcept {

                    assert(scalar_vec.size() == DIM);
                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(scalar_vec[i] / other.data_[i]); 
                    return result;

                }


                /**
                 * @brief Equality operator
                 * 
                 * @param other: vector as l-value const reference
                 * 
                 * @return bool
                 * 
                 */
                constexpr bool operator==(const vector& other) const noexcept {

                    for (size_t i{}; i < DIM; ++i) if (data_[i] != other.data_[i]) return false;
                    return true;

                }


                /**
                 * @brief Inequality operator
                 * 
                 * @param other: vector as l-value const reference
                 * 
                 * @return bools
                 *  
                 */
                constexpr bool operator!=(const vector& other) const noexcept {

                    for (size_t i{}; i < DIM; ++i) if (data_[i] != other.data_[i]) return true;
                    return false;

                }
            

                /**
                 * @brief Access the i-th element of the vector
                 * 
                 * @param index: size_t
                 * 
                 * @return constexpr measurement 
                 * 
                 * @note: index must be in the range [0, DIM)
                 */
                constexpr const measurement& operator[](const size_t& index) const { 
                    
                    if (index >= DIM) throw std::out_of_range("Cannot access a vector with an index out of range");
                    return data_[index]; 
                    
                }


                /**
                 * @brief Access the i-th element of the vector
                 * 
                 * @param index: size_t 
                 * 
                 * @return constexpr measurement& 
                 * 
                 * @note: index must be in the range [0, DIM)
                 * 
                 */
                constexpr measurement& operator[](const size_t& index) { 
                    
                    if (index >= DIM) throw std::out_of_range("Cannot access a vector with an index out of range");
                    return data_[index]; 
                    
                }


                /**
                 * @brief Print the vector to the standard output
                 * 
                 * @param os: ostream as l-value reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return std::ostream&
                 * 
                 */
                friend std::ostream& operator<<(std::ostream& os, const vector& vec) noexcept {

                    os << "{";
                    for (size_t i{}; i < DIM; ++i) os << vec.data_[i] << (i != DIM - 1 ? ", " : "}");
                    return os;

                }


                /**
                 * @brief Print the vector to file
                 * 
                 * @param file: ofstream as l-value reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return std::ofstream&
                 * 
                 */
                friend std::ofstream& operator<<(std::ofstream& file, const vector& vec) noexcept {

                    file << "{";
                    for (size_t i{}; i < DIM; ++i) file << vec.data_[i] << (i != DIM - 1 ? ", " : "}");
                    return file;

                }


            // =============================================
            // operations
            // =============================================

                /**
                 * @brief Compute the cross product between two vectors
                 * 
                 * @param v1: vector as l-value const reference
                 * @param v2: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 * 
                 */
                constexpr friend vector cross(const vector& v1, const vector& v2) {

                    std::vector<measurement> cross_vec;
                    cross_vec.reserve(DIM); 
                    for (size_t i{}; i < DIM; ++i) cross_vec.emplace_back(v1[(i + 1) % v1.size()] * v2[(i + 2) % v1.size()] - v1[(i + 2) % v1.size()] * v2[(i + 1) % v1.size()]); 
                    return cross_vec;
                
                }

                
                /**
                 * @brief Compute the dot product between two vectors
                 * 
                 * @param v1: vector as l-value const reference
                 * @param v2: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 * 
                 */
                constexpr friend measurement dot(const vector& v1, const vector& v2) noexcept {

                    measurement result(v1[0].units() * v2[0].units());
                    for (size_t i{}; i < v1.size(); ++i) result += v1[i] * v2[i]; 
                    return result;
                
                }


                /**
                 * @brief Invert the vector
                 * 
                 * @return constexpr vector 
                 * 
                 */
                constexpr vector inv() const {

                    std::vector<measurement> result; 
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i].inv()); 
                    return result;

                }


                /**
                 * @brief Take the power of the vector
                 * 
                 * @param power
                 *  
                 * @return constexpr vector 
                 * 
                 */
                constexpr vector pow(const int& power) const noexcept {

                    std::vector<measurement> result;
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i].pow(power));
                    return result;

                }


                /**
                 * @brief Take the square of the vector
                 * 
                 * @return constexpr vector 
                 * 
                 */
                constexpr vector square() const noexcept {

                    std::vector<measurement> result;
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i].square());
                    return result;

                }


                /**
                 * @brief Take the cube of the vector
                 * 
                 * @return constexpr vector 
                 * 
                 */
                constexpr vector cube() const noexcept {

                    std::vector<measurement> result;
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i].cube());
                    return result;

                }


                /**
                 * @brief Take the root of the vector
                 * 
                 * @param power 
                 * @return constexpr vector 
                 */
                constexpr vector root(const int& power) const {

                    std::vector<measurement> result;
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i].root(power));
                    return result;

                }


                /**
                 * @brief Take the square root of the vector
                 * 
                 * @return constexpr vector 
                 * 
                 */
                constexpr vector sqrt() const {

                    std::vector<measurement> result;
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i].sqrt());
                    return result;

                }


                /**
                 * @brief Take the cube root of the vector
                 * 
                 * @return constexpr vector 
                 * 
                 */
                constexpr vector cbrt() const {

                    std::vector<measurement> result;
                    result.reserve(DIM);
                    for (size_t i{}; i < DIM; ++i) result.emplace_back(data_[i].cbrt());
                    return result;

                }


            // =============================================
            // set & get methods
            // =============================================

                /**
                 * @brief Get the size of the vector
                 * 
                 * @return constexpr size_t 
                 * 
                 */
                constexpr size_t size() const noexcept { 
                    
                    return DIM; 
                    
                }


                /**
                 * @brief Get the first element of the vector
                 * 
                 * @return constexpr measurement
                 * 
                 * @note the vector must have at least one element
                 * 
                 */
                constexpr measurement x() const noexcept requires (DIM >= 1) { 
                    
                    return data_[0]; 
                
                }
                                

                /**
                 * @brief Get the second element of the vector
                 * 
                 * @return constexpr measurement
                 * 
                 * @note the vector must have at least two elements
                 * 
                 */
                constexpr measurement y() const noexcept requires (DIM >= 2) { 
                    
                    return data_[1]; 
                
                }


                /**
                 * @brief Get the third element of the vector
                 * 
                 * @return constexpr measurement
                 * 
                 * @note the vector must have at least three elements
                 * 
                 */
                constexpr measurement z() const noexcept requires (DIM >= 3) { 
                    
                    return data_[2]; 
                
                }


                /**
                 * @brief Get the first element of the vector
                 * 
                 * @return constexpr measurement&
                 * 
                 * @note the vector must have at least one element
                 * 
                 */
                constexpr measurement& x() noexcept requires (DIM >= 1) { 
                    
                    return data_[0]; 
                
                }
                                

                /**
                 * @brief Get the second element of the vector
                 * 
                 * @return constexpr measurement&
                 * 
                 * @note the vector must have at least two elements
                 * 
                 */
                constexpr measurement& y() noexcept requires (DIM >= 2) { 
                    
                    return data_[1]; 
                
                }


                /**
                 * @brief Get the third element of the vector
                 * 
                 * @return constexpr measurement&
                 * 
                 * @note the vector must have at least three elements
                 * 
                 */
                constexpr measurement& z() noexcept requires (DIM >= 3) { 
                    
                    return data_[2]; 
                
                }


                /**
                 * @brief Get the data of the vector
                 * 
                 * @return constexpr std::vector<measurement> 
                 * 
                 */
                constexpr std::vector<measurement> data() const noexcept { 
                    
                    return data_; 
                
                }


                /**
                 * @brief Get the vector object
                 * 
                 * @return constexpr vector<DIM> 
                 * 
                 */
                constexpr vector<DIM> as_vector() const noexcept { 
                    
                    return *this; 
                
                }


                /**
                 * @brief Get the norm of the vector
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement norm() const noexcept { 

                    if constexpr (DIM == 1) return data_[0];
                    measurement norm; 
                    for (size_t i{}; i < DIM; ++i) norm += data_[i].square();
                    return norm.sqrt(); 

                }


                /**
                 * @brief Get the squared norm of the vector
                 * 
                 * @return constexpr measurement 
                 * 
                 */
                constexpr measurement norm2() const noexcept { 

                    if constexpr (DIM == 1) return data_[0]; 
                    measurement norm(0.0 * data_.front().units().square()); 
                    for (size_t i{}; i < DIM; ++i) norm += data_[i].square();
                    return norm; 
                    
                }


                /**
                 * @brief Get the normalization of the vector
                 * 
                 * @return constexpr vector 
                 * 
                 */
                constexpr vector versor() const {
                    
                    vector result;
                    measurement norm = this->norm();
                    for (size_t i{}; i < DIM; ++i) result[i] = data_[i] / norm;
                    return result; 

                } 


                /**
                 * @brief Get the polar angle
                 * 
                 * @return constexpr measurement
                 * 
                 * @note the vector must have at least two elements
                 * 
                 */
                constexpr measurement phi() const noexcept requires (DIM >= 2) { 
                    
                    return std::atan2(data_[1].value(), data_[0].value()) * rad; 
                
                }     
                

                /**
                 * @brief Get the azimuthal angle
                 * 
                 * @return constexpr measurement
                 * 
                 * @note the vector must have at least three elements
                 * 
                 */
                constexpr measurement theta() const requires (DIM >= 3) { 

                    if (data_[2] == 0.0 * m) return 0.0 * rad;
                    return std::acos((data_[2] / norm()).value()) * rad;

                }
                

                /// @brief Print the vector to the standard output
                constexpr void print() const noexcept {

                    std::cout << "{\n";
                    for (size_t i{}; i < DIM; ++i) std::cout << "\t" << data_[i] << "\n";
                    std::cout << "}\n";

                }   


                /**
                 * @brief Save the vector to a file
                 * 
                 * @param file_name: the name of the file
                 * @param units: desired units for the output
                 */
                void save(const std::string& file_name, const unit& units) const {

                    std::ofstream file_out(file_name, std::ios::app);
                    if (file_out.is_open()) {
                        for (size_t i{}; i < DIM; ++i) {
                            file_out << data_[i].value_as(units); 
                            if (i < DIM - 1) file_out << ' '; 
                        }
                    } else std::cerr << "PROBLEM: Unable to open '" << file_name << "'\n";
                    file_out << '\n';
                    file_out.close();

                }


            protected:

            // =============================================
            // class members
            // =============================================

                /// @brief Vector data
                std::vector<measurement> data_;


        }; // class vector


        /// @brief Class expressing a generic matrix of vectors of measurements in a n-dimentional system
        template <size_t rows, size_t cols = rows> 
        class matrix {


            public:

                // =============================================
                // class members
                // =============================================

                size_t rows_; 

                size_t cols_;

                using mat = std::vector<vector<rows>>;

                mat data_;


                // =============================================
                // constructors & destructor
                // =============================================

                explicit constexpr matrix() noexcept {

                    assert(rows != 0 && cols != 0);
                    rows_ = rows;
                    cols_ = cols;
                    for (size_t i{}; i < cols_; i++) data_.emplace_back(vector<rows>());    

                }

                
                constexpr matrix(const mat& data) noexcept {

                    assert(data.size() == cols);
                    rows_ = rows;
                    cols_ = cols;
                    for (auto& i : data) {
                        assert(i.size() == rows);
                        data_.emplace_back(i);
                    }

                }


                // =============================================
                // operators
                // =============================================

                constexpr matrix operator=(const matrix& other) noexcept {

                    assert(other.rows_ != 0 || other.cols_ != 0);
                    rows_ = other.rows_;
                    cols_ = other.cols_;
                    for (size_t i{}; i < other.rows_; i++) data_[i] = other.data_[i]; 
                    return *this;

                }
                

                constexpr matrix operator+=(const matrix& other) noexcept {

                    for(size_t i{}; i < cols_; i++) data_[i] += other.data_[i];
                    return *this; 

                }



                constexpr matrix operator-=(const matrix& other) noexcept {

                    for(size_t i{}; i < cols_; i++) data_[i] -= other.data_[i];
                    return *this; 

                }


                constexpr matrix operator+(const matrix& other) const noexcept {

                    assert(rows_ == other.rows_ && cols_ == other.cols_);
                    matrix result;
                    for(size_t i{}; i < cols_; i++) result.data_[i] = data_[i] + other.data_[i];
                    return result;

                }


                constexpr matrix operator-(const matrix& other) const noexcept {

                    assert(rows_ == other.rows_ && cols_ == other.cols_);
                    matrix result;
                    for(size_t i{}; i < cols_; i++) result[i] = data_[i] - other.data_[i];
                    return result;

                }


                constexpr matrix operator-() const noexcept {

                    matrix result;
                    for(size_t i{}; i < cols_; i++) result[i] = -data_[i];
                    return result;

                }


                constexpr matrix operator*(const std::vector<measurement>& vec) const noexcept {

                    assert(vec.size() == cols_);
                    matrix result;
                    for (size_t j{}; j < cols_; j++) result.data_[j] = data_[j] * vec[j];
                    return result;

                }

                constexpr matrix operator*(const measurement& meas) const noexcept {
                    
                    assert(meas.value() != 0.0); 
                    matrix result;
                    for (size_t j{}; j < cols_; j++) result.data_[j] = data_[j] * meas;
                    return result;

                }


                friend constexpr matrix operator*(const measurement& meas, const matrix& mat) noexcept {
                    
                    assert(meas.value() != 0.0); 
                    matrix result;
                    for (size_t j{}; j < mat.cols_; j++) result.data_[j] = meas * mat.data_[j];
                    return result;

                }


                constexpr matrix operator/(const measurement& meas) const noexcept {
                    
                    assert(meas.value() != 0.0);
                    matrix result;
                    for (size_t i{}; i < rows_; i++) 
                        for (size_t j{}; j < cols_; j++) 
                            result.data_[i][j] = data_[i][j] / meas;
                    
                    return result;

                }


                friend constexpr matrix operator/(const measurement& meas, const matrix& mat) noexcept {
                    
                    assert(meas.value() != 0.0); 
                    matrix result;
                    for (size_t j{}; j < mat.cols_; j++) result.data_[j] = meas / mat.data_[j];
                    return result;

                }


                constexpr matrix operator*(const double& scalar) const noexcept {
                    
                    assert(scalar != 0.0); 
                    matrix result;
                    for (size_t i{}; i < cols_; i++) result.data_[i] = data_[i] * scalar;
                    return result;

                }


                constexpr matrix operator/(const double& scalar) const noexcept {
                    
                    assert(scalar != 0.0);
                    matrix result;
                    for (size_t i{}; i < cols_; i++) result.data_[i] = data_[i] * scalar;
                    return result;

                }


                // =============================================
                // matrix methods
                // =============================================

                constexpr size_t rows_size() const noexcept { 
                    
                    return rows_; 
                    
                }


                constexpr size_t cols_size() const noexcept { 
                    
                    return cols_; 
                    
                }


                constexpr vector<rows> operator[](const size_t& i) const noexcept { 
                    
                    return data_[i]; 
                    
                }


                constexpr mat data() const noexcept { 
                    
                    return data_; 
                    
                }


                constexpr vector<rows> column(const size_t& n_col) const noexcept {
                    
                    return data_[n_col]; 
                    
                }


                constexpr vector<cols> row(const size_t& m_row) const noexcept {

                    vector<cols> appo; 
                    for (size_t i{}; i < cols_; i++) appo[i] = data_[i][m_row];
                    return appo; 

                }


                constexpr matrix<cols, rows> transpose() const noexcept {

                    matrix<cols, rows> result;
                    for (size_t i{}; i < cols_; i++) result.data_[i] = row(i);
                    return result;

                }


                constexpr matrix<rows, cols> as_matrix() const noexcept {

                    return *this;

                }


                constexpr matrix<rows, cols> units() const noexcept {

                    matrix<rows, cols> result;
                    for (size_t i{}; i < cols_; i++) {
                        for (size_t j{}; j < rows_; j++) {
                            result.data_[i][j] = measurement(0.0, data_[i][j].units());
                        }
                    }
                    return result;

                }


                constexpr void print() const noexcept {

                    std::cout << "matrix = {\n";
                    for (size_t i{}; i < rows_; i++) {
                        for (size_t j{}; j < cols_; j++) { 
                            std::cout << "\t"; 
                            data_[j][i].print(false); 
                        
                        }
                        std::cout << "\n"; 
                    }
                    std::cout << "}\n";

                }


        }; // class matrix


        /// @brief Class expressing a time measurement
        class time {


            public:

            // =============================================
            // constructors and destructor
            // =============================================   

                explicit time() noexcept : 

                    time_(0.0, s) {}


                explicit constexpr time(const double& time, const unit& unit = s) {

                    if (time < 0.0) throw std::invalid_argument("Time cannot be negative");
                    if (unit.base() != base::second) throw std::invalid_argument("Wrong time unit, the unit_base must be seconds");
                    time_ = measurement(time, unit);
                    
                }


                constexpr time(const measurement& time) {

                    if (time.value() < 0.0) throw std::invalid_argument("Time cannot be negative");
                    if (time.units().base() != base::second) throw std::invalid_argument("Wrong time unit, the unit_base must be seconds");
                    time_ = time;

                }


                constexpr time(const time& other) noexcept :
                    
                    time_{other.time_} {}   

                
                constexpr time(time&& other) noexcept :
                    
                    time_{std::move(other.time_)} {}
                

                ~time() = default; 


                // =============================================
                // operators
                // =============================================

                constexpr time operator=(const time& other) {
                    
                    time_ = other.time_;
                    return *this;

                }


                constexpr time& operator=(time&& other) {

                    time_ = std::move(other.time_);
                    return *this;

                }


                constexpr time operator+=(const time& other) {

                    time_ += std::move(other.time_);
                    return *this;

                }


                constexpr time operator-=(const time& other) {

                    time_ -= std::move(other.time_);
                    return *this;

                }


                constexpr time operator+(const time& other) const { 
                    
                    return { time_ + other.time_ }; 
                    
                } 


                constexpr time operator-(const time& other) const { 
                    
                    return { time_ - other.time_ }; 
                    
                } 


                constexpr bool operator==(const time& other) const { 
                    
                    return { time_ == other.time_ }; 
                    
                }


                constexpr bool operator<(const time& other) const { 
                    
                    return time_ < other.time_; 
                
                } 


                constexpr bool operator>(const time& other) const { 
                    
                    return time_ > other.time_; 
                
                } 


                constexpr bool operator<=(const time& other) const { 
                    
                    return time_ <= other.time_; 
                
                } 


                constexpr bool operator>=(const time& other) const { 
                    
                    return time_ >= other.time_; 
                
                } 


                // =============================================
                // set & get & print methods
                // =============================================          
                            
                constexpr void increase(const time& other) { 
                    
                    time_ += other.time_; 
                
                } 

                
                constexpr void reset() { 
                    
                    time_.value() = 0.0; 
                
                }

    
                constexpr measurement as_measurement() const { 
                    
                    return time_; 
                
                }

                
                constexpr time as_time() const { 
                    
                    return *this; 
                
                }


                inline void print() const { 
                    
                    time_.print(); 
                
                }


                void save(const std::string& file_name, const unit& units = s) const {

                    std::ofstream file_out(file_name, std::ios::app);
                    if (file_out.is_open()) file_out << time_.value_as(units) << " "; 
                    else std::cerr << "PROBLEM: Unable to open '" << file_name << "'\n";
                    file_out.close();

                }


            private:     

            // =============================================
            // class members
            // =============================================     

                measurement time_; 


        }; // class time
 

        /**
         * @brief Class expressing the position of a generic object as a vector of measurements
         * 
         * @tparam DIM: size_t dimentions of the system
         * 
         */
        template <size_t DIM> 
        class position : public vector<DIM> {


            public: 

                // =============================================
                // constructors & destructor
                // =============================================

                explicit constexpr position() noexcept : 
                
                    vector<DIM>(base::metre) {}


                constexpr position(const std::vector<measurement>& pos) noexcept : 
                    
                    vector<DIM>(pos) {

                        for (auto& i : pos) assert(i.units().base() == base::metre);
                    
                    }


                constexpr position(const vector<DIM>& pos) noexcept : 
                    
                    vector<DIM>(pos) {

                        for (auto& i : pos.data()) assert(i.units().base() == base::metre);
                    
                    }


                constexpr position(const position& other) noexcept : 
                    
                    vector<DIM>(other) {}


                ~position() = default; 


                // =============================================
                // operators
                // =============================================

                constexpr position& operator=(const position& other) noexcept {
                    
                    *this = other; 
                    return *this;

                }


                constexpr position& operator+=(const position& other) noexcept {

                    this += other;
                    return *this;

                }


                constexpr position& operator-=(const position& other) noexcept {

                    this -= other;
                    return *this;

                }


                constexpr position operator+(const position& other) const noexcept { 
                    
                    return *this + other.as_vector();

                }


                constexpr position operator-(const position& other) const noexcept { 
                    
                    return *this - other.as_vector();

                }


                constexpr position operator-() const noexcept { 
                    
                    return -*this;

                }


                constexpr vector<DIM> operator*(const measurement& meas) noexcept {

                    return *this * meas;

                }


                constexpr vector<DIM> operator/(const measurement& meas) noexcept {

                    return *this / meas;

                }


                friend vector<DIM> operator*(const measurement& meas, const position& pos) noexcept {

                    return meas * pos.as_vector();

                }


                friend vector<DIM> operator/(const measurement& meas, const position& pos) noexcept {

                    return meas / pos.as_vector();
                    
                }


                constexpr position operator*=(const double& val) noexcept {

                    this *= val;
                    return *this;

                }


                constexpr position operator/=(const double& val) noexcept {

                    this /= val;
                    return *this;

                }


                constexpr position operator*(const double& val) const noexcept { 
                    
                    return position(this * val); 
                
                }


                constexpr position operator/(const double& val) const noexcept { 
                    
                    return position(this / val); 
                
                }


                friend constexpr position operator*(const double& val, const position& pos) noexcept { 
                    
                    return val * pos.as_vector();
                
                }


                friend constexpr vector<DIM> operator/(const double& val, const position& pos) noexcept { 
                    
                    return val / pos.as_vector();
                
                }


                // =============================================
                // position methods
                // =============================================

                constexpr measurement distance(const position& other) const noexcept {    
                    
                    return (*this - other).norm();
                    
                }       


                constexpr measurement distance2(const position& other) const noexcept {    
                    
                    return (*this - other).norm2();
                    
                }       


                constexpr measurement phi() const noexcept {

                    return vector<DIM>::phi();

                }

                    
                constexpr measurement phi(const position& other) const noexcept requires (DIM >= 2) { 

                    return measurement(std::atan2((other[1] - this[1]).value(), (other[0] - this[0]).value()), rad);
                
                }


                constexpr measurement theta() const noexcept {

                    return vector<DIM>::theta();

                }


                constexpr measurement theta(const position& other) const noexcept requires (DIM >= 3) {

                    if (other.pos_[2] == this[2]) return 0 * rad;
                    return std::acos(((other[2] - this[2]) / this.distance(other)).value()) * rad; 
                
                }


                constexpr std::vector<double> versor() const noexcept {

                    return this.versor();

                }

                
                constexpr std::vector<double> versor(const position& other) const noexcept requires (DIM >= 2) {

                    if constexpr (DIM == 2) return { std::cos(phi(other).value()), std::sin(phi(other).value()) };
                    if constexpr (DIM == 3) return { std::cos(phi(other).value()), std::sin(phi(other).value()), (other[2] - this[2]) / this.distance(other)};

                } 


                constexpr void print() const noexcept {

                    std::cout << "position = ";
                    vector<DIM>::print();

                }   


        }; // class position


        // class expressing the linear velocity of an generic object 
        // as a vector of measurements (base::metre / base::second) in a n-dimentional system
        template <size_t DIM> 
        class linear_velocity : public vector<DIM> {

            
            public: 

                // =============================================
                // constructors & destructor
                // =============================================

                explicit constexpr linear_velocity() noexcept : 
                
                    vector<DIM>(base::metre / base::second) {}


                constexpr linear_velocity(const std::vector<measurement>& vel) noexcept : 
                    
                    vector<DIM>(vel) {

                        for (auto& i : vel) assert(i.units().base() == base::metre / base::second);

                    }


                constexpr linear_velocity(const vector<DIM>& vel) noexcept : 
                    
                    vector<DIM>(vel) {

                        for (auto& i : vel.data()) assert(i.units().base() == base::metre / base::second);

                    }


                constexpr linear_velocity(const linear_velocity& other) noexcept : 
                    
                    vector<DIM>(other) {}


                ~linear_velocity() = default; 


                // =============================================
                // operators
                // =============================================

                constexpr linear_velocity operator=(const linear_velocity& other) noexcept {

                    vector<DIM>::operator=(other);
                    return *this;

                }


                constexpr linear_velocity operator+=(const linear_velocity& other) noexcept {

                    vector<DIM>::operator+=(other);
                    return *this;

                }


                constexpr linear_velocity operator-=(const linear_velocity& other) noexcept {

                    vector<DIM>::operator-=(other);
                    return *this;

                }


                constexpr linear_velocity operator+(const linear_velocity& other) const noexcept { 
                    
                    return vector<DIM>::operator+(other.as_vector());

                }


                constexpr linear_velocity operator-(const linear_velocity& other) const noexcept { 
                    
                    return vector<DIM>::operator-(other.as_vector());

                }


                constexpr linear_velocity operator-() const noexcept { 
                    
                    return vector<DIM>::operator-();

                }


                constexpr vector<DIM> operator*(const measurement& meas) noexcept {

                    return vector<DIM>::operator*(meas);

                }


                constexpr vector<DIM> operator/(const measurement& meas) noexcept {

                    return vector<DIM>::operator/(meas);

                }


                friend vector<DIM> operator*(const measurement& meas, const linear_velocity& pos) noexcept {

                    return meas * pos.as_vector();

                }


                friend vector<DIM> operator/(const measurement& meas, const linear_velocity& pos) noexcept {

                    return meas / pos.as_vector();
                    
                }


                constexpr linear_velocity operator*=(const double& val) noexcept {

                    vector<DIM>::operator*=(val);
                    return *this;

                }


                constexpr linear_velocity operator/=(const double& val) noexcept {

                    vector<DIM>::operator/=(val);
                    return *this;

                }


                constexpr linear_velocity operator*(const double& val) const noexcept { 
                    
                    return linear_velocity(vector<DIM>::operator*(val)); 
                
                }


                constexpr linear_velocity operator/(const double& val) const noexcept { 
                    
                    return linear_velocity(vector<DIM>::operator/(val)); 
                
                }


                friend constexpr linear_velocity operator*(const double& val, const linear_velocity& vel) noexcept { 
                    
                    return val * vel.as_vector();
                
                }


                friend constexpr vector<DIM> operator/(const double& val, const linear_velocity& vel) noexcept { 
                    
                    return val / vel.as_vector();
                
                }


                friend constexpr linear_velocity operator/(const position<DIM>& pos, const time& t) noexcept {

                    return pos.as_vector() / t.as_measurement();

                }


                // =============================================
                // print method
                // =============================================

                constexpr void print() const {

                    std::cout << "linear velocity = ";
                    vector<DIM>::print();

                }   


        }; // class linear_velocity
    

        // class expressing the linear acceleration of an generic object 
        // as a vector of measurements (base::metre / base::second.square()) in a n-dimentional system
        template <size_t DIM> 
        class linear_acceleration : public vector<DIM> {

            
            public: 

                // =============================================
                // constructors & destructor
                // =============================================

                explicit constexpr linear_acceleration() noexcept : 
                
                    vector<DIM>(base::metre / base::second.square()) {}


                constexpr linear_acceleration(const std::vector<measurement>& vel) noexcept : 
                    
                    vector<DIM>(vel) {

                        for (auto& i : vel) assert(i.units().base() == base::metre / base::second.square());

                    }


                constexpr linear_acceleration(const vector<DIM>& vel) noexcept : 
                    
                    vector<DIM>(vel) {

                        for (auto& i : vel.data()) assert(i.units().base() == base::metre / base::second.square());

                    }


                constexpr linear_acceleration(const linear_acceleration& other) noexcept : 
                    
                    vector<DIM>(other) {}


                ~linear_acceleration() = default; 


                // =============================================
                // operators
                // =============================================

                constexpr linear_acceleration operator=(const linear_acceleration& other) noexcept {

                    vector<DIM>::operator=(other);
                    return *this;

                }


                constexpr linear_acceleration operator+=(const linear_acceleration& other) noexcept {

                    vector<DIM>::operator+=(other);
                    return *this;

                }


                constexpr linear_acceleration operator-=(const linear_acceleration& other) noexcept {

                    vector<DIM>::operator-=(other);
                    return *this;

                }


                constexpr linear_acceleration operator+(const linear_acceleration& other) const noexcept { 
                    
                    return vector<DIM>::operator+(other.as_vector());

                }


                constexpr linear_acceleration operator-(const linear_acceleration& other) const noexcept { 
                    
                    return vector<DIM>::operator-(other.as_vector());

                }


                constexpr linear_acceleration operator-() const noexcept { 
                    
                    return vector<DIM>::operator-();

                }


                constexpr vector<DIM> operator*(const measurement& meas) noexcept {

                    return vector<DIM>::operator*(meas);

                }


                constexpr vector<DIM> operator/(const measurement& meas) noexcept {

                    return vector<DIM>::operator/(meas);

                }


                friend vector<DIM> operator*(const measurement& meas, const linear_acceleration& pos) noexcept {

                    return meas * pos.as_vector();

                }


                friend vector<DIM> operator/(const measurement& meas, const linear_acceleration& pos) noexcept {

                    return meas / pos.as_vector();
                    
                }


                constexpr linear_acceleration operator*=(const double& val) noexcept {

                    vector<DIM>::operator*=(val);
                    return *this;

                }


                constexpr linear_acceleration operator/=(const double& val) noexcept {

                    vector<DIM>::operator/=(val);
                    return *this;

                }


                constexpr linear_acceleration operator*(const double& val) const noexcept { 
                    
                    return linear_acceleration(vector<DIM>::operator*(val)); 
                
                }


                constexpr linear_acceleration operator/(const double& val) const noexcept { 
                    
                    return linear_acceleration(vector<DIM>::operator/(val)); 
                
                }


                friend constexpr linear_acceleration operator*(const double& val, const linear_acceleration& acc) noexcept { 
                    
                    return val * acc.as_vector();
                
                }


                friend constexpr vector<DIM> operator/(const double& val, const linear_acceleration& acc) noexcept { 
                    
                    return val / acc.as_vector();
                
                }


                friend constexpr linear_acceleration<DIM> operator/(const linear_velocity<DIM>& vel, const time& t) noexcept {

                    return vel.as_vector() / t.as_measurement();

                }   


                // =============================================
                // print method
                // =============================================

                constexpr void print() const {

                    std::cout << "linear acceleration = ";
                    vector<DIM>::print();

                }   


        }; // class linear_acceleration


        // class expressing the force of an generic object 
        // as a vector of measurements (base::kilogram * base::metre / base::second.square()) in a n-dimentional system
        template <size_t DIM> 
        class force : public vector<DIM> {

            
            public: 

                // =============================================
                // constructors & destructor
                // =============================================

                explicit constexpr force() noexcept : 
                
                    vector<DIM>(base::kilogram * base::metre / base::second.square()) {}


                constexpr force(const std::vector<measurement>& vel) noexcept : 
                    
                    vector<DIM>(vel) {

                        for (auto& i : vel) assert(i.units().base() == base::kilogram * base::metre / base::second.square());

                    }


                constexpr force(const vector<DIM>& vel) noexcept : 
                    
                    vector<DIM>(vel) {

                        for (auto& i : vel.data()) assert(i.units().base() == base::kilogram * base::metre / base::second.square());

                    }


                constexpr force(const force& other) noexcept : 
                    
                    vector<DIM>(other) {}


                ~force() = default; 


                // =============================================
                // operators
                // =============================================

                constexpr force operator=(const force& other) noexcept {

                    vector<DIM>::operator=(other);
                    return *this;

                }


                constexpr force operator+=(const force& other) noexcept {

                    vector<DIM>::operator+=(other);
                    return *this;

                }


                constexpr force operator-=(const force& other) noexcept {

                    vector<DIM>::operator-=(other);
                    return *this;

                }


                constexpr force operator+(const force& other) const noexcept { 
                    
                    return vector<DIM>::operator+(other.as_vector());

                }


                constexpr force operator-(const force& other) const noexcept { 
                    
                    return vector<DIM>::operator-(other.as_vector());

                }


                constexpr force operator-() const noexcept { 
                    
                    return vector<DIM>::operator-();

                }


                constexpr vector<DIM> operator*(const measurement& meas) noexcept {

                    return vector<DIM>::operator*(meas);

                }


                constexpr vector<DIM> operator/(const measurement& meas) noexcept {

                    return vector<DIM>::operator/(meas);

                }


                friend vector<DIM> operator*(const measurement& meas, const force& pos) noexcept {

                    return meas * pos.as_vector();

                }


                friend vector<DIM> operator/(const measurement& meas, const force& pos) noexcept {

                    return meas / pos.as_vector();
                    
                }


                constexpr force operator*=(const double& val) noexcept {

                    vector<DIM>::operator*=(val);
                    return *this;

                }


                constexpr force operator/=(const double& val) noexcept {

                    vector<DIM>::operator/=(val);
                    return *this;

                }


                constexpr force operator*(const double& val) const noexcept { 
                    
                    return force(vector<DIM>::operator*(val)); 
                
                }


                constexpr force operator/(const double& val) const noexcept { 
                    
                    return force(vector<DIM>::operator/(val)); 
                
                }


                friend constexpr force operator*(const double& val, const force& acc) noexcept { 
                    
                    return val * acc.as_vector();
                
                }


                friend constexpr vector<DIM> operator/(const double& val, const force& acc) noexcept { 
                    
                    return val / acc.as_vector();
                
                }


                // =============================================
                // print method
                // =============================================

                constexpr void print() const {

                    std::cout << "force = ";
                    vector<DIM>::print();

                }   


        }; // class force


        template <size_t DIM>
        constexpr position<DIM> operator*(const linear_velocity<DIM>& vel, const time& t) noexcept {

            return vel.as_vector() * t.as_measurement();

        }


        template <size_t DIM>
        constexpr linear_velocity<DIM> operator*(const linear_acceleration<DIM>& acc, const time& t) noexcept {

            return acc.as_vector() * t.as_measurement();

        }


        /// @brief Class for timing the execution of a generic function/process
        class timer {

            
            public:

                // =============================================
                // constructor and destructor
                // =============================================   

                /**
                 * @brief Construct a new timer object
                 * 
                 */
                constexpr timer() = default;


                /**
                 * @brief Destroy the timer object
                 * 
                 */
                ~timer() = default;

                
                // =============================================
                // timer methods
                // =============================================   

                /**
                 * @brief start the timer
                 * 
                 */
                inline void start() { 
                    
                    start_ = std::chrono::high_resolution_clock::now(); 
                    
                }


                /**
                 * @brief stop the timer
                 * 
                 */
                inline void pause() { 
                    
                    pause_ = std::chrono::high_resolution_clock::now(); 
                    
                }
                

                constexpr measurement elapsed(const unit& units = s) const {

                    if (units != s) return measurement(s.convert(static_cast<std::chrono::duration<double>>(pause_ - start_).count(), units), units);
                    else return measurement(static_cast<std::chrono::duration<double>>(pause_ - start_).count(), units);

                }


                /**
                 * @brief print the elapsed time
                 * 
                 */
                inline void print(const unit& units = s) const { 

                    std::cout << "elapsed time = " << elapsed(units) << "\n"; 

                }


            protected: 

                // =============================================
                // class members
                // =============================================     
                
                std::chrono::time_point<std::chrono::high_resolution_clock> start_, pause_;


        }; // class timer


    } // namespace tools


} // namespace physics



namespace std {


    using measurement = physics::measurements::measurement; 

    
    std::vector<measurement> operator*(const std::vector<measurement>& vec, const double& value) {
        
        std::vector<measurement> result; 
        result.reserve(vec.size()); 
        for (size_t i{}; i < result.size(); ++i) result[i] = vec[i] * value; 
        return result; 

    }

    
    std::vector<measurement> operator/(const std::vector<measurement>& vec, const double& value) {
        
        std::vector<measurement> result; 
        result.reserve(vec.size()); 
        for (size_t i{}; i < result.size(); ++i) result[i] = vec[i] / value; 
        return result; 

    }


}


namespace math {


    namespace equations {
        

        using namespace physics::tools; 


        template <size_t DIM, typename... Args>
        struct ode {
            

            std::function<vector<DIM>(const vector<DIM>&, const Args&...)> eval_; 


            std::function<vector<DIM>(const vector<DIM>&, const Args&...)> diff_; 


            constexpr ode(const std::function<vector<DIM>(const vector<DIM>&, const Args&...)>& evaluate, 
                          const std::function<vector<DIM>(const vector<DIM>&, const Args&...)>& differentiate) noexcept :

                eval_(evaluate), 
                diff_(differentiate) {}

            
            constexpr vector<DIM> eval(const vector<DIM>& init, const Args&... param) const noexcept {

                return eval_(init, param...); 

            }


            constexpr vector<DIM> diff(const vector<DIM>& init, const Args&... param) const noexcept {

                return diff_(init, param...); 

            } 


        }; // struct ode


        template <size_t DIM, size_t N_EQ, typename... Args>
        struct ode_system {


            std::function<matrix<DIM, N_EQ>(const matrix<DIM, N_EQ>&, const Args&...)> diff_; 


            constexpr ode_system(const std::function<matrix<DIM, N_EQ>(const matrix<DIM, N_EQ>&, const Args&...)>& differentiate) noexcept :

                diff_(differentiate) {}


            constexpr matrix<DIM, N_EQ> diff(const matrix<DIM, N_EQ>& init, const Args&... param) const noexcept {

                return diff_(init, param...); 

            } 


        }; // struct ode_system


        template <size_t DIM, size_t N_EQ, typename... Args>
        struct ode_solver {


            template <class ODE, class Type, class Incr>
            constexpr Type euler(const ODE& ode, const Type& init, const Args&... param, const Incr& incr) const noexcept {

                return init + ode.diff(init, param...) * incr;
                
            }


            template <class ODE, class Type, class Incr>
            constexpr Type RK4(const ODE& ode, const Type& init, const Args&... param, const Incr& incr) const noexcept {

                matrix<DIM, N_EQ> k1(ode.diff(init, param...));
                matrix<DIM, N_EQ> k2(ode.diff(init + k1 * (incr / 2.), param...)); 
                matrix<DIM, N_EQ> k3(ode.diff(init + k2 * (incr / 2.), param...));
                matrix<DIM, N_EQ> k4(ode.diff(init + k3 * incr, param...)); 
                return init + (k1 + k2 * 2. + k3 * 2. + k4) * (incr / 6.); 

            }


        }; // struct ode_solver


        template <size_t DIM, typename... Args>
        class hamiltonian {

            
            private: 
            

                ode<DIM, Args...>* potential_;


                ode_system<DIM, 2, Args...> system_;


                ode_solver<DIM, 2, Args...> solver_;


            public: 
                
                
                constexpr hamiltonian(ode<DIM, Args...>* potential) noexcept :

                    potential_{potential}, 
                    system_{
                    
                        [this](const matrix<DIM, 2>& init, const Args&... param) -> matrix<DIM, 2> { 
                        
                            return matrix<DIM, 2>({init[1], - potential_->diff(init[0], param...)}); 

                        }

                    } {}

/*
                constexpr hamiltonian(ode<DIM>* potential) noexcept {

                    potential_.emplace_back(potential);

                }


                constexpr hamiltonian(const std::vector<ode<DIM>*>& potentials) noexcept : 

                    potential_(potentials) {}


                constexpr void add_potential(ode<DIM>* pot) noexcept {

                    potential_.emplace_back(pot);

                }


                constexpr void add_potential(const std::vector<ode<DIM>*>& potentials) noexcept {

                    potential_.insert(potential_.end(), potentials.begin(), potentials.end());

                }


                constexpr void remove_potential(ode<DIM>* pot) noexcept {

                    potential_.erase(std::remove(potential_.begin(), potential_.end(), pot), potential_.end());

                }


                constexpr vector<DIM> eval(const measurement& mass, const position<DIM>& pos, const linear_velocity<DIM>& vel) const noexcept {

                    vector<DIM> pot_eval(base::kilogram * base::metre.square() / base::second.square());
                    for (auto& pot : potential_) pot_eval += pot->eval(pos.as_vector());
                    return mass * vel.square() / 2. + pot_eval;

                }

*/


                constexpr matrix<DIM, 2> solve(const measurement& mass, const matrix<DIM, 2>& init, const Args&... param, const physics::tools::time& incr) const noexcept {
                    
                    return solver_.RK4(system_, init, param..., std::vector<measurement>({incr.as_measurement(), incr.as_measurement() / mass}));

                }


        }; // class hamiltonian


    } // namespace equations 


} // namespace math


namespace physics {


    /// @brief Namespace constains some physical constants
    namespace constants {


        using namespace physics::measurements; 

        constexpr measurement G = measurement(6.67430e-17, N * km.square() / kg.square());
        
        constexpr measurement electron_mass = measurement(9.109383701528e-31, kg);     
        
        constexpr measurement proton_mass = measurement(1.672621923695e-27, kg);
        
        constexpr measurement electron_charge = measurement(-1.602176634e-19, C);
        
        constexpr measurement proton_charge = measurement(1.602176634e-19, C); 


    } // namespace constants


    namespace potentials {


        using namespace tools;  
        using namespace math::equations;       

        
        template <size_t DIM>
        class gravitational_potential : public ode<DIM, measurement, measurement, position<DIM>> {


            public:


                constexpr gravitational_potential() noexcept : 
                    
                    ode<DIM, measurement, measurement, position<DIM>>{
                        
                        [this](const vector<DIM>& init, const measurement& mass, const measurement& source_mass, const position<DIM>& source_position) -> vector<DIM> { 
                            
                            return - physics::constants::G * mass * source_mass * (init - source_position.as_vector()) / (init - source_position.as_vector()).norm2();
                            
                        }, 

                        [this](const vector<DIM>& init, const measurement& mass, const measurement& source_mass, const position<DIM>& source_position) -> vector<DIM> { 

                            return physics::constants::G * mass * source_mass * (init - source_position.as_vector()) / (init - source_position.as_vector()).norm().cube();

                        }

                    } {}



        }; // class gravitational_potential


    } // namespace potential


    namespace objects {


        using namespace physics::tools;
        using namespace physics::potentials; 


        template <size_t DIM>
        class object {


            protected:
                    
                // =============================================
                // class members
                // =============================================
                
                uint32_t id_;

                std::string type_;

                std::string name_;

                // shape* shape_;

                position<DIM> position_;

                linear_velocity<DIM> linear_velocity_;

                linear_acceleration<DIM> linear_acceleration_;

                force<DIM> force_; 

                    
            public: 

                // =============================================
                // constructors and destructor 
                // =============================================

                explicit constexpr object(const std::string& type,
                                          // shape shape = geometry::shapes::circle(0.0 * m),
                                          const position<DIM>& pos, 
                                          const linear_velocity<DIM>& vel = linear_velocity<DIM>(), 
                                          const linear_acceleration<DIM>& acc = linear_acceleration<DIM>(), 
                                          const force<DIM>& F = force<DIM>(),
                                          const uint32_t& id = 0, 
                                          const std::string& name = "") noexcept : 

                    id_{id}, 
                    type_{type},
                    name_{name},
                    // shape_(std::move(&shape)), 
                    position_{pos},
                    linear_velocity_{vel}, 
                    linear_acceleration_{acc},
                    force_{F} {}


                constexpr object(const object& other) noexcept : 
                
                    id_{other.id_}, 
                    type_{other.type_},
                    name_{other.name_},
                    // shape_{other.shape_},
                    position_{other.position_}, 
                    linear_velocity_{other.linear_velocity_}, 
                    linear_acceleration_{other.linear_acceleration_},
                    force_{other.force_} {

                        std::cout << "object copy constructor called\n";

                    } 

    /*
                    constexpr object(object&& other) noexcept : 
                    
                        id_{std::move(other.id_)}, 
                        type_{std::move(other.type_)},
                        name_{std::move(other.name_)},
                        // shape_{std::move(other.shape_)},
                        position_{std::move(other.position_)}, 
                        linear_velocity_{std::move(other.linear_velocity_)}, 
                        linear_acceleration_{std::move(other.linear_acceleration_)},
                        force_{std::move(other.force_)} {

                            std::cout << "object move constructor called\n";

                        } 
    */
                
                ~object() noexcept {

                    std::cout << "object destructor called\n";
                    
                }


                // =============================================
                // operators
                // =============================================

                constexpr object& operator=(const object& other) noexcept {

                    std::cout << "object copy assignment operator called\n";

                    if (this != &other) {

                        id_ = other.id_;
                        type_ = other.type_;
                        name_ = other.name_;
                        // shape_ = other.shape_;
                        position_ = other.position_;
                        linear_velocity_ = other.linear_velocity_;
                        linear_acceleration_ = other.linear_acceleration_;
                        force_ = other.force_;

                    }

                    return *this;

                }


                constexpr object& operator=(object&& other) noexcept {

                    std::cout << "object move assignment operator called\n";

                    if (this != &other) {

                        id_ = std::move(other.id_);
                        type_ = std::move(other.type_);
                        name_ = std::move(other.name_);
                        // shape_ = std::move(other.shape_);
                        position_ = std::move(other.position_);
                        linear_velocity_ = std::move(other.linear_velocity_);
                        linear_acceleration_ = std::move(other.linear_acceleration_);
                        force_ = std::move(other.force_);

                    }

                    return *this;

                }


                constexpr bool operator==(const object& other) const noexcept {

                    return (id_ == other.id_ && 
                            type_ == other.type_ && 
                            name_ == other.name_ && 
                            // shape_ == other.shape_ && 
                            position_ == other.position_ && 
                            linear_velocity_ == other.linear_velocity_ && 
                            linear_acceleration_ == other.linear_acceleration_ && 
                            force_ == other.force_);

                }


                constexpr bool operator!=(const object& other) const noexcept {

                    return !(*this == other);

                }


                // =============================================
                // object methods
                // =============================================

                constexpr void add_force(const force<DIM>& other) {

                    force_ += other; 

                }


                constexpr void add_force(const vector<DIM>& other) {
                    
                    for (auto i : other.data()) assert(i.units() == newton);
                    force_ += other; 

                }


                constexpr void reset_force() { 
                    
                    force_ = force<DIM>(); 
                    
                }


                constexpr virtual void update_linear_acceleration() {}


/*
                template <size_t N_EQ>
                constexpr void move(const ode<DIM, N_EQ>& ode, const tools::time& dt) {
                    
                    set_state(solver_.RK4(ode, get_state(), dt.as_measurement()));
                    update_linear_acceleration(); 

                }
*/

                // =============================================
                // set methods
                // =============================================

                constexpr void set_id(const uint32_t& id) { 
                    
                    id_ = id; 
                    
                }


                constexpr void set_type(const char* type) { 
                    
                    type_ = type; 
                    
                }


                constexpr void set_name(const char* name) { 
                    
                    name_ = name; 
                    
                }

/*
                constexpr void shape(shape_2D* shape) { 
                    
                    shape_ = shape; 
                    
                }
*/


                constexpr void set_position(const std::vector<measurement>& pos) { 
                    
                    position_ = pos; 
                    
                }

                                
                constexpr void set_position(const position<DIM>& pos) { 
                    
                    position_ = pos; 
                    
                }


                constexpr void set_linear_velocity(const std::vector<measurement>& vel) { 
                    
                    linear_velocity_ = vel; 
                    
                }
                

                constexpr void set_linear_velocity(const linear_velocity<DIM>& other) { 
                    
                    linear_velocity_ = other; 
                    
                }


                inline void set_linear_acceleration(const std::vector<measurement>& lin_acc) { 
                    
                    linear_acceleration_ = lin_acc; 
                    
                }


                constexpr void set_linear_acceleration(const linear_acceleration<DIM>& lin_acc) { 
                    
                    linear_acceleration_ = lin_acc; 
                    
                }


                constexpr void set_state(const matrix<DIM, 2>& state) {

                    position_ = state[0];
                    linear_velocity_ = state[1];
                    // linear_acceleration_ = state[2]; 

                }


                // =============================================
                // get and print methods
                // =============================================

                constexpr uint32_t id() const { 
                    
                    return id_;
                    
                }

                
                constexpr std::string type() const { 
                    
                    return type_;
                    
                }


                constexpr std::string name() const { 
                    
                    return name_;
                    
                }

/*
                inline shape& shape() const { 
                    
                    return *shape_;
                
                }
*/


                constexpr position<DIM> get_position() const { 
                    
                    return position_; 
                
                }


                constexpr linear_velocity<DIM> get_linear_velocity() const { 
                    
                    return linear_velocity_; 
                
                }


                constexpr linear_acceleration<DIM> get_linear_acceleration() const { 
                    
                    return linear_acceleration_; 
                
                }


                constexpr matrix<DIM, 2> get_state() const { 

                    return matrix<DIM, 2>({position_.as_vector(), linear_velocity_.as_vector()}); 

                }


                constexpr force<DIM> get_force() const { 
                    
                    return force_; 
                    
                }


                constexpr size_t dim() const {

                    return DIM; 

                }


                constexpr object<DIM> as_object() const { 
                    
                    return *this; 
                    
                }


                constexpr void print() const {

                    std::cout << "\nobject:\n";
                    std::cout << "type = " << type_ << "\n";
                    if (id_ != 0) std::cout << "id = " << id_ << "\n";
                    if (! name_.empty()) std::cout << "name = " << name_ << "\n"; 
                    // shape_->print(); 
                    position_.print(); 
                    linear_velocity_.print(); 
                    linear_acceleration_.print();
                    force_.print(); 

                }               


        }; // class object        


        template <size_t DIM>
        class mass : public object<DIM> {


            protected: 

                // =============================================
                // class member
                // =============================================

                measurement mass_; 

                bool gravitational_field_;


            public: 

                // =============================================
                // constructors and destructor
                // =============================================

                explicit constexpr mass(const measurement& mass,
                                        const position<DIM> pos, 
                                        const uint32_t& id,
                                        const std::string& name,
                                        const bool& gravity = true) noexcept :

                    object<DIM>("mass", pos, linear_velocity<DIM>(), linear_acceleration<DIM>(), force<DIM>(), id, name),
                    mass_{mass},
                    gravitational_field_{gravity} {}


                explicit constexpr mass(const measurement& mass,
                                        // shape shape = geometry::shapes::circle(0.0 * m),
                                        const position<DIM>& pos, 
                                        const linear_velocity<DIM>& lin_vel = linear_velocity<DIM>(), 
                                        const linear_acceleration<DIM>& lin_acc = linear_acceleration<DIM>(),
                                        const force<DIM>& F = force<DIM>(),
                                        const uint32_t& id = 0,
                                        const std::string& name = "",
                                        const std::string& type = "mass",
                                        const bool& gravity = true) noexcept : 
                    
                    object<DIM>(type, pos, lin_vel, lin_acc, F, id, name),
                    mass_(mass),
                    gravitational_field_(gravity) {

                        std::cout << "mass constructor called\n";

                    }


                explicit constexpr mass(const measurement& mass,
                                        const object<DIM>& obj,
                                        const bool& gravity = true) noexcept : 
                                 
                    object<DIM>(obj),
                    mass_(mass),
                    gravitational_field_(gravity) {

                        std::cout << "mass constructor called\n";

                    }


                constexpr mass(const mass& other) noexcept : 

                    object<DIM>(other.as_object()),
                    mass_(other.mass_),
                    gravitational_field_(other.gravitational_field_) {

                        std::cout << "mass copy constructor called\n";

                    }


                ~mass() noexcept {
                        
                    std::cout << "mass destructor called\n";
    
                }


                // =============================================
                // operators
                // =============================================

                constexpr mass operator=(const mass& other) {

                    mass_ = other.mass_; 
                    gravitational_field_ = other.gravitational_field_; 
                    this->id_ = other.id_;
                    this->type_ = other.type_;
                    this->name_ = other.name_; 
                    // this->shape_ = other.shape_; 
                    this->position_ = other.position_; 
                    this->linear_velocity_ = other.linear_velocity_; 
                    this->linear_acceleration_ = other.linear_acceleration_; 
                    return *this; 

                }


                constexpr mass operator=(mass&& other) {

                    mass_ = other.mass_; 
                    gravitational_field_ = other.gravitational_field_; 
                    this->id_ = other.id_;
                    this->type_ = other.type_;
                    this->name_ = other.name_;
                    // this->shape_ = other.shape_; 
                    this->position_ = other.position_; 
                    this->linear_velocity_ = other.linear_velocity_; 
                    this->linear_acceleration_ = other.linear_acceleration_; 
                    return *this; 

                }


                // =============================================
                // set and get methods
                // =============================================

                constexpr void set_mass(const measurement& mass) { 
                    
                    assert(mass.value() >= 0.0);
                    assert(mass.units().base() == base::kilogram);
                    mass_ = mass; 
                    
                }
                

                constexpr measurement get_mass() const { 
                    
                    return mass_; 
                    
                }


                constexpr void gravitational_field(const bool& gravity) { 
                    
                    gravitational_field_ = gravity; 
                    
                }


                constexpr bool gravitational_field() const { 
                    
                    return gravitational_field_; 
                    
                }


                constexpr vector<DIM> momentum() const { 
                    
                    return mass_ * this->linear_velocity_.as_vector();
                    
                }


                constexpr std::vector<measurement> angular_momentum() const { 
                    
                    return mass_ * cross(this->position_.as_vector(), this->linear_velocity_.as_vector()); 
                    
                }


                constexpr measurement kinetic_energy() const { 
                    
                    return 0.5 * mass_ * this->linear_velocity_.norm2(); 
                    
                }


                constexpr mass<DIM> as_mass() const { 
                    
                    return *this; 
                    
                }


                constexpr linear_acceleration<DIM> get_gravitational_pull(const position<DIM>& pos) const { 
                    
                    return linear_acceleration<DIM>(this->position_.direction(pos) * (- constants::G * mass_ / this->position_.distance2(pos))); 
                    
                }


                constexpr void gravitate(const mass& other) {

                    if (gravitational_field_ == true && other.gravitational_field_ == true) {
                        add_force(this->position_.direction(other.get_position()) * (constants::G * mass_ * other.mass_) / this->position_.distance2(other.position_));
                    }

                }


                constexpr void update_linear_acceleration() {

                    this->linear_acceleration_ = this->force_ / mass_; 
                    //  + this->solver_.RK4(equations_, get_linear_acceleration(), t.as_measurement()); 

                }
                

                // =============================================
                // print methods
                // =============================================

                void print_momentum() const {

                    std::cout << "momentum = {\n"; 
                    for (auto p : momentum()) {
                        std::cout << "\t"; 
                        p.print(); 
                    }
                    std::cout << "}\n";

                }


                void print_kinetic_energy() const {

                    std::cout << "kinetic energy = "; 
                    kinetic_energy().print(); 

                }


                constexpr void print() const {

                    std::cout << "\nobject: \n"; 
                    std::cout << "type = " << this->type_ << "\n"; 
                    if (this->id_ != 0) std::cout << "id = " << this->id_ << "\n";
                    if (! this->name_.empty()) std::cout << "name = " << this->name_ << "\n";
                    // shape_->print(); 
                    std::cout << "mass = "; 
                    mass_.print();
                    this->position_.print(); 
                    this->linear_velocity_.print(); 
                    this->linear_acceleration_.print();
                    this->force_.print();

                }     


        }; // class mass


        template <class T>
        class system {


            protected:
            
                // =============================================
                // class members
                // =============================================     

                std::vector<T> bodies_; 


            public:

                // =============================================
                // constructors and destructor
                // =============================================     

                explicit constexpr system() noexcept {}


                constexpr system(const std::vector<T>& objs) noexcept : 
                
                    bodies_{objs} {}


                ~system() = default; 


                // =============================================
                // objects methods
                // =============================================     
                
                constexpr void add(const T& other) noexcept { 
                    
                    bodies_.emplace_back(other); 
                
                } 


                constexpr void reset() noexcept { 
                    
                    bodies_.clear(); 
                
                }


                constexpr size_t count() const { 
                    
                    return bodies_.size(); 
                    
                }


                constexpr std::vector<T> objects() const { 
                    
                    return bodies_; 
                    
                }


                constexpr T element(const uint& pos) const { 
                    
                    return bodies_[pos]; 
                    
                }
                

                constexpr T operator[](const uint& pos) const { 
                    
                    return bodies_[pos]; 
                    
                }


                constexpr T& operator[](const uint& pos) { 
                    
                    return bodies_[pos]; 
                    
                }


                constexpr virtual void update() {} 


                constexpr virtual void evolve(const tools::time& dt) {}


        }; // class system

        
        template <size_t DIM>
        class mass_system : public system<mass<DIM>> {


            private: 

                // =============================================
                // class members
                // =============================================     

                hamiltonian<DIM, measurement, measurement, position<DIM>> hamiltonian_{dynamic_cast<ode<DIM, measurement, measurement, position<DIM>>*>(new gravitational_potential<DIM>())};


                constexpr position<DIM> get_center_of_mass(const position<DIM>& SR_center, const measurement& initial_mass) const {

                    vector<DIM> center_of_mass(base::metre * base::kilogram);
                    for (auto& i : system<mass<DIM>>::bodies_) if (i.get_position() != SR_center) center_of_mass += i.get_position() * i.get_mass();
                    return center_of_mass / (get_total_mass() - initial_mass);

                }


            public:


                constexpr mass_system() noexcept : 
                    
                    system<mass<DIM>>{} {}

                
                constexpr mass_system(const std::vector<mass<DIM>>& objs) noexcept : 
                    
                    system<mass<DIM>>{objs} {}


                constexpr measurement get_total_mass() const {

                    measurement total_mass(0.0 * kg);
                    for (auto& i : system<mass<DIM>>::bodies_) total_mass += i.get_mass();
                    return total_mass;

                }


                constexpr position<DIM> get_center_of_mass() const {

                    vector<DIM> center_of_mass(base::metre * base::kilogram);
                    for (auto& i : system<mass<DIM>>::bodies_) center_of_mass += i.get_position() * i.get_mass();
                    return center_of_mass / get_total_mass();

                }


                constexpr void evolve(const tools::time& dt = (1 / 90.) * s) override {

                    for (auto& obj : system<mass<DIM>>::bodies_) {
                        
                        obj.set_state(hamiltonian_.solve(obj.get_mass(), obj.get_state(), obj.get_mass(), get_total_mass() - obj.get_mass(), get_center_of_mass(obj.get_position(), obj.get_mass()), dt));

                    }

                }


                constexpr void print() const {

                    for (auto& i : system<mass<DIM>>::bodies_) i.print();

                }


        }; // class mass_system


    }


} // namespace physics


namespace math {
    

    // namespace defining some tools for math
    namespace tools {


        using namespace physics::tools;


        // namespace defining some descriptive statistic functions
        namespace descriptive_statistics {
                

            // using measurement = physics::measurements::measurement;


            double mean(const std::vector<double>& v) {

                if (v.size() == 0) throw std::invalid_argument("Can't operate a descriptive statistic funtion on an empty vector"); 
                return std::accumulate(v.begin(), v.end(), 0.) / v.size();

            }


            // double median(const std::vector<double>& v) {

            //     if (v.size() == 0) throw std::invalid_argument("Can't operate a descriptive statistic funtion on an empty vector"); 
            //     std::vector<double> appo = v; 
            //     if (std::is_sorted(appo.begin(), appo.end()) == false) std::sort(appo.begin(), appo.end());
            //     if (v.size() % 2 != 0) return appo[v.size() / 2];
            //     else return (appo[v.size() / 2] + appo[(v.size() / 2) - 1]) / 2; 

            // }


            double variance(const std::vector<double>& v) {

                double average = mean(v);
                double accu{}; 
                for (auto x : v) accu += std::pow(x - average, 2); 
                return accu / v.size();

            }


            inline double sd(const std::vector<double>& v) { 
                
                return std::sqrt(variance(v)); 
                
            }


            // inline double sdom(const std::vector<double>& v) { 
                
            //     return std::sqrt(variance(v) / v.size()); 
                
            // }


            /**
             * @brief Compute the mean value of a vector of measurements, 
             *        the uncertainty is computed as the standard deviation of mean (sdom)
             * 
             * @param vec: vector of measurements
             * 
             * @note 
             * @return uncertain_measurement
             */
            uncertain_measurement mean(const std::vector<measurement>& vec) {

                if (vec.size() == 0) throw std::invalid_argument("Can't operate a descriptive statistic funtion on an empty vector"); 
                
                measurement average = std::accumulate(vec.begin(), vec.end(), measurement(0., vec[0].units())) / vec.size();
                measurement sigma_sq = measurement(0., vec[0].units().square()); 
                for (auto x : vec) sigma_sq += (x - average).square(); 

                return uncertain_measurement(average, sigma_sq.sqrt() / vec.size());                 

            }


            /**
             * @brief Compute the median of a vector of measurements
             * 
             * @param vec: vector of measurements
             * 
             * @return measurement 
             */
            measurement median(const std::vector<measurement>& v) {

                if (v.size() == 0) throw std::invalid_argument("Can't operate a descriptive statistic funtion on an empty vector"); 
                std::vector<measurement> copy{v}; 
                if (std::is_sorted(copy.begin(), copy.end()) == false) std::sort(copy.begin(), copy.end());
                if (v.size() % 2 != 0) return copy[v.size() / 2];
                else return (copy[v.size() / 2] + copy[(v.size() / 2) - 1]) / 2; 

            }


            /**
             * @brief Compute the variance of a vector of measurements
             * 
             * @param vec: vector of measurements
             * 
             * @return measurement 
             */
            measurement variance(const std::vector<measurement>& vec) {

                measurement average = mean(vec).as_measurement();
                measurement sigma_sq = measurement(0., vec[0].units().square()); 
                for (auto x : vec) sigma_sq += (x - average).square(); 
                return sigma_sq / vec.size();

            }


            /**
             * @brief Compute the standard deviation of a vector of measurements
             * 
             * @param vec: vector of measurements
             * 
             * @return measurement 
             */
            inline measurement standard_dev(const std::vector<measurement>& vec) { 
                
                return variance(vec).sqrt(); 
            
            }


            /**
             * @brief Compute the standard error of the mean of a vector of measurements
             * 
             * @param vec: vector of measurements
             * 
             * @return measurement 
             */
            inline measurement sdom(const std::vector<measurement>& vec) { 
                
                return (variance(vec) / vec.size()).sqrt(); 
            
            }

            
            /**
             * @brief Compute the weighted mean of a vector of measurements
             * 
             * @param vec 
             * 
             * @return uncertain_measurement 
             */
            uncertain_measurement wmean(const std::vector<uncertain_measurement>& vec) {

                if (vec.size() == 0) throw std::invalid_argument("Can't operate a descriptive statistic funtion on an empty vector"); 
                
                measurement weighted = measurement(0., vec[0].units().inv());
                measurement weights = measurement(0., vec[0].units().inv().square());
                for (auto& x : vec) {
                    weighted += x.as_measurement() * x.weight(); 
                    weights += x.weight();
                }

                return uncertain_measurement(weighted / weights, weights.inv().sqrt());

            }


            /**
             * @brief Compute the weighted variance of a vector of measurements
             * 
             * @param vec 
             * 
             * @return measurement 
             */
            measurement wvariance(const std::vector<uncertain_measurement>& vec) {

                if (vec.size() == 0) throw std::invalid_argument("Can't operate a descriptive statistic funtion on an empty vector"); 
                measurement weights = measurement(0., vec[0].units().inv().square());
                for (auto& x : vec) weights += x.weight(); 
                return weights.inv(); 

            }


            /**
             * @brief Compute the weighted standard deviation of a vector of measurements
             * 
             * @param vec 
             * 
             * @return measurement 
             */
            inline measurement wsd(const std::vector<uncertain_measurement>& vec) {

                return wvariance(vec).sqrt();
            
            }


            measurement chi_sq(const std::vector<measurement>& v, 
                               const std::vector<measurement>& expected) {

                if (v.size() != expected.size()) throw std::invalid_argument("Can't operate a chi square funtion on vectors of different size"); 
                measurement accu = measurement(0., v[0].units()); 
                for (size_t i{}; i < v.size(); ++i) accu += (v[i] - expected[i]).square() / expected[i]; 
                return accu; 

            }         


            // chi squared reduced of an std::vector<measurement>
            inline measurement chi_sq_r(const std::vector<measurement>& v, 
                                        const std::vector<measurement>& expected, 
                                        const int& gdl) {

                return chi_sq(v, expected) / gdl; 

            }


            class linear_regression {


                public:


                    linear_regression() noexcept : 
                        
                        intercept_{},
                        slope_{} {}


                    void train(const std::vector<measurement>& xData, 
                               const std::vector<measurement>& yData) {

                        if (xData.size() != yData.size() || (xData.size() == 0 && yData.size() == 0)) 
                            throw std::runtime_error("Can't operate a linear regression training session with empty data sets or data sets of different sizes.");

                        measurement sum_X(0.0, xData.front().units()); 
                        measurement sum_XX(0.0, xData.front().units().square()); 
                        measurement sum_Y(0.0, yData.front().units()); 
                        measurement sum_XY(0.0, xData.front().units() * yData.front().units());
                        measurement delta(0.0, xData.front().units().square()); 
                        measurement sigma_y(0.0, yData.front().units()); 

                        size_t N{xData.size()};      
                        for (size_t i{}; i < N; ++i) {

                            sum_X += xData[i];
                            sum_XX += xData[i].square(); 
                            sum_Y += yData[i];
                            sum_XY += xData[i] * yData[i]; 

                        }

                        delta = (N * sum_XX - sum_X.square());
                        slope_ = (N * sum_XY - sum_X * sum_Y) / delta;
                        intercept_ = (sum_XX * sum_Y - sum_X * sum_XY) / delta;

                        for (size_t i{}; i < N; ++i) 
                            sigma_y += (yData[i] - intercept_.as_measurement() - slope_.as_measurement() * xData[i]).square(); 

                        sigma_y /= (N - 2);

                        slope_.uncertainty() = (N * sigma_y / delta).sqrt().value();
                        intercept_.uncertainty() = (sigma_y * sum_XX / delta).sqrt().value();

                    }


                    void train(const std::vector<measurement>& xData, 
                               const std::vector<uncertain_measurement>& yData) {

                        if (xData.size() != yData.size() || (xData.size() == 0 && yData.size() == 0)) 
                            throw std::runtime_error("Can't operate a linear regression training session with empty data sets or data sets of different sizes.");

                        measurement weight(0.0, yData.front().units().inv().square());
                        measurement wsum(0.0, weight.units());
                        measurement wsum_X(0.0, xData.front().units()); 
                        measurement wsum_Y(0.0, yData.front().units()); 
                        measurement wsum_XX(0.0, xData.front().units().square()); 
                        measurement wsum_XY(0.0, xData.front().units() * yData.front().units());
                        measurement delta(0.0, xData.front().units().square()); 

                        size_t N{xData.size()};      
                        for (size_t i{}; i < N; ++i) {

                            weight = yData[i].uncertainty_as_measurement().inv().square();
                            wsum += weight; 
                            wsum_X += xData[i].as_measurement() * weight;
                            wsum_Y += yData[i].as_measurement() * weight;
                            wsum_XX += xData[i].as_measurement().square() * weight; 
                            wsum_XY += xData[i].as_measurement() * yData[i].as_measurement() * weight; 

                        }

                        delta = (wsum * wsum_XX - wsum_X.square());
                        slope_ = uncertain_measurement((wsum * wsum_XY - wsum_X * wsum_Y) / delta, (wsum / delta).sqrt().value());
                        intercept_ = uncertain_measurement((wsum_XX * wsum_Y - wsum_X * wsum_XY) / delta, (wsum_XX / delta).sqrt().value());

                    }


                    void train(const std::vector<uncertain_measurement>& xData, 
                               const std::vector<uncertain_measurement>& yData, 
                               const measurement& sigma_y_from_x = measurement(0.0, unitless)) {

                        if (xData.size() != yData.size() || (xData.size() == 0 && yData.size() == 0)) 
                            throw std::runtime_error("Can't operate a linear regression training session with empty data sets or data sets of different sizes.");
                       
                        measurement weight(0.0, (xData.front().units() * yData.front().units()).inv().square());
                        measurement wsum(0.0, weight.units());
                        measurement wsum_X(0.0, xData.front().units()); 
                        measurement wsum_Y(0.0, yData.front().units()); 
                        measurement wsum_XX(0.0, xData.front().units().square()); 
                        measurement wsum_XY(0.0, xData.front().units() * yData.front().units());
                        measurement delta(0.0, xData.front().units().square()); 

                        size_t N{xData.size()};      
                        for (size_t i{}; i < N; ++i) {

                            weight = ((xData[i].uncertainty_as_measurement() * sigma_y_from_x).square() + yData[i].uncertainty_as_measurement().square()).inv().square();
                            wsum += weight; 
                            wsum_X += xData[i].as_measurement() * weight;
                            wsum_Y += yData[i].as_measurement() * weight;
                            wsum_XX += xData[i].as_measurement().square() * weight; 
                            wsum_XY += xData[i].as_measurement() * yData[i].as_measurement() * weight; 

                        }

                        delta = (wsum * wsum_XX - wsum_X.square());
                        slope_ = uncertain_measurement((wsum * wsum_XY - wsum_X * wsum_Y) / delta, (wsum / delta).sqrt().value());
                        intercept_ = uncertain_measurement((wsum_XX * wsum_Y - wsum_X * wsum_XY) / delta, (wsum_XX / delta).sqrt().value());

                    }


                    constexpr uncertain_measurement predict(const double& x) const noexcept {

                        return intercept_ + slope_ * x;

                    }


                    constexpr uncertain_measurement intercept() const noexcept {

                        return intercept_;
                    
                    }


                    constexpr uncertain_measurement slope() const noexcept {

                        return slope_;
                    
                    }


                private:


                    uncertain_measurement intercept_, slope_;


            };
            

        } // namespace descriptive_statistics 


        // class random_generator for generating pseudo-casual numbers
        class random_generator {


            private: 

                // =============================================
                // class members
                // =============================================        
                
                uint32_t m1{}, m2{}, m3{}, m4{}; 

                uint32_t l1{}, l2{}, l3{}, l4{}; 
                
                uint32_t n1{}, n2{}, n3{}, n4{};


            public: 
        
                // =============================================
                // constructor & destructor
                // =============================================

                constexpr random_generator() { 
                    
                    this->set_up(); 
                    
                }

                ~random_generator() {

                    this->save_seed(); 

                }


                // =============================================
                // set & get methods
                // =============================================

                constexpr void set_seed(const uint32_t * s, const uint32_t& p1, const uint32_t& p2) {
                    m1 = 502;
                    m2 = 1521;
                    m3 = 4071;
                    m4 = 2107;
                    l1 = s[0];
                    l2 = s[1];
                    l3 = s[2];
                    l4 = s[3];
                    l4 = 2 * (l4 / 2) + 1;
                    n1 = 0;
                    n2 = 0;
                    n3 = p1;
                    n4 = p2;
                }


                void set_up() {

                    uint32_t seed[4];
                    uint32_t p1, p2;

                    std::ifstream file_in("../random/primes.in");
                    if (file_in.is_open()) file_in >> p1 >> p2 ;
                    else throw std::runtime_error("Unable to open file '../random/primes.in'");
                    file_in.close();

                    file_in.open("../random/seed.in");
                    std::string property;
                    if (file_in.is_open()) {

                        while (!file_in.eof()) {
                            
                            file_in >> property;
                            if (property == "RANDOMSEED") {

                                file_in >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                                this->set_seed(seed, p1, p2);

                            }

                        }
                        file_in.close();
                    } else std::cerr << "PROBLEM: Unable to open seed.in\n";

                }


                void save_seed() const {  

                    std::ofstream file_out("../random/seed.out");
                    if (file_out.is_open()) file_out << l1 << "\t" << l2 << "\t" << l3 << "\t" << l4 << "\n";
                    else std::cerr << "PROBLEM: Unable to open seed.out\n";
                    file_out.close();

                }


                // =============================================
                // distributions methods
                // =============================================

                constexpr double rannyu() {

                    const double twom12{0.000244140625};
                    int i1{}, i2{}, i3{}, i4{};
                    i1 = l1 * m4 + l2 * m3 + l3 * m2 + l4 * m1 + n1;
                    i2 = l2 * m4 + l3 * m3 + l4 * m2 + n2;
                    i3 = l3 * m4 + l4 * m3 + n3;
                    i4 = l4 * m4 + n4;
                    l4 = i4 % 4096;
                    i3 = i3 + i4 / 4096;
                    l3 = i3 % 4096;
                    i2 = i2 + i3 / 4096;
                    l2 = i2 % 4096;
                    l1 = (i1 + i2 / 4096) % 4096;
                    return twom12 * (l1 + twom12 * (l2 + twom12 * (l3 + twom12 * l4)));

                }


                constexpr double unif(const double& min, 
                                      const double& max) {

                    return min + std::fabs(max - min) * rannyu(); 

                }


                constexpr double exp(const double& lambda) {

                    return - std::log(1 - rannyu()) / lambda; 

                }


                constexpr double lorentzian(const double& mean, 
                                            const double& gamma) {

                    return mean + gamma * std::tan(constants::pi * (rannyu() - 0.5));

                }


                constexpr double gauss_box_muller(const double& mean, 
                                                  const double& sigma) {

                    return mean + sigma * std::sqrt(-2 * std::log(rannyu())) * std::cos(2 * constants::pi * rannyu());

                }


                constexpr double gauss_accept_reject(const double& mean, 
                                                     const double& sigma) {

                    double x{}, y{}, g{}; 
                    while (true) {
                        x = unif(mean - 3. * sigma, mean + 3. * sigma); 
                        y = rannyu(); 
                        g = exp(std::pow(x, 2) / 2); 
                        if (y >= g) break;
                    }
                    return x;

                }


        }; // class tools::random_generator


        // class integral for evaluating the integrals
        class integral {


            private: 

                // =============================================
                // class members
                // =============================================

                double a_{}, b_{}, h_{};
                
                uint32_t steps_{};

                int32_t sign_{}; 

                double sum_{}, integral_{}, old_integral_{}, error_{};  

                random_generator rg_;


                // =============================================
                // set methods
                // =============================================

                constexpr void steps(const uint32_t& n) { 

                    steps_ = n; 
                    h_ = std::fabs(b_ - a_) / steps_;

                }        

                
                constexpr void check_range() { 
                    
                    sign_ = (a_ == b_ ? 0. : (b_ > a_ ? 1. : -1)); 
                    
                }
                
                
                constexpr void sum(const double& sum) { 
                    
                    sum_ = sum; 
                
                }

                
                constexpr void reset_integral() { 
                    
                    integral_ = 0; 
                
                }    

                
                constexpr void begin_integration(const double& a, 
                                                 const double& b, 
                                                 const uint32_t& n = 1000, 
                                                 const double& sum0 = 0) {

                    a_ = a; 
                    b_ = b; 
                    check_range(); 
                    steps(n);
                    reset_integral(); 
                    sum(sum0); 

                }


            public: 

                // =============================================
                // constructors
                // =============================================

                constexpr integral() noexcept {

                    rg_.set_up(); 

                }

                
                ~integral() {} 

            
                // =============================================
                // get methods
                // =============================================

                constexpr double a() const { 
                    
                    return a_; 
                    
                }


                constexpr double b() const { 
                    
                    return b_; 
                    
                }


                constexpr int sign() const { 
                    
                    return sign_; 
                    
                }


                constexpr uint32_t steps() const { 
                    
                    return steps_; 
                    
                }


                constexpr double h() const { 
                    
                    return h_; 
                    
                }


                constexpr double sum() const { 
                    
                    return sum_; 
                    
                }


                constexpr double value() const { 
                    
                    return integral_; 
                    
                }


                constexpr double error() const { 
                    
                    return error_; 
                    
                }



                // =============================================
                // print methods
                // =============================================

                void print_value(const double& precision = 1.e-6) {

                    std::cout << "integral of f(x) in [" << a_ << ", " << b_ << "] = " << std::setprecision(precision) << integral_ << "\n";

                }

                void print_error(const double& precision = 1.e-6) {

                    std::cout << "error = " << std::setprecision(precision) << error_ << "\n";

                }        

                void print_integral(const double& precision = 1.e-6) {

                    print_value(precision); 
                    print_error(precision); 

                }
                
                
                // =============================================
                // integration methods
                // =============================================

                constexpr void midpoint(const double& a, 
                                        const double& b, 
                                        const std::function<double(double)>& f, 
                                        const uint32_t& n = 1000) {

                    begin_integration(a, b, n); 
                    for (uint32_t i{}; i < steps_; i++) { sum_ += f(a_ + (i + 0.5) * h_); }
                    integral_ = sum_ * h_; 

                }


                constexpr void midpoint_fixed(const double& a, 
                                              const double& b, 
                                              const std::function<double(double)>& f, 
                                              const double& prec = 1.e-6) {

                    double old_integral_2{}, old_integral_3{};
                    begin_integration(a, b, 1); 
                    while (true) {
                        old_integral_3 = old_integral_2;
                        old_integral_2 = old_integral_;
                        old_integral_ = integral_; 
                        midpoint(a_, b_, f, steps_ * 2);
                        error_ = 64 * std::fabs(64 * integral_ - 84 * old_integral_ + 21 * old_integral_2 - old_integral_3) / 2835; // errore al sesto ordine
                        if (error_ < prec) break;
                    }  
                    integral_ = (4096 * integral_ - 1344 * old_integral_ + 84 * old_integral_2 - old_integral_3) / 2835; // integrale all'ottavo ordine

                }

                
                void trapexoid(const double& a, 
                                         const double& b, 
                                         const std::function<double(double)>& f, 
                                         const uint32_t& n = 1000) {

                    begin_integration(a, b, n, (f(a) + f(b)) / 2.);
                    for (uint32_t i{1}; i < steps_; i++) sum_ += f(a_ + i * h_); 
                    integral_ = sum_ * h_; 

                } 


                void trapexoid_fixed(const double& a, 
                                               const double& b, 
                                               const std::function<double(double)>& f, 
                                               const double& prec = 1.e-6) {

                    double old_integral_2{}, old_integral_3{};
                    begin_integration(a, b, 2, f(a) + f(b) / 2. + f((a + b) / 2.)); 
                    while (true) {
                        old_integral_3 = old_integral_2;
                        old_integral_2 = old_integral_;
                        old_integral_ = integral_; 
                        trapexoid(a_, b_, f, steps_ * 2);
                        error_ = 64 * std::fabs(64 * integral_ - 84 * old_integral_ + 21 * old_integral_2 - old_integral_3) / 2835; // errore al sesto ordine
                        if (error_ < prec) break;
                    }
                    integral_ = (4096 * integral_ - 1344 * old_integral_ + 84 * old_integral_2 - old_integral_3) / 2835; // integrale all'ottavo ordine

                }


                void simpson(const double& a, 
                                       const double& b, 
                                       const std::function<double(double)>& f, 
                                       const uint32_t& n = 1000) {

                    if (n % 2 == 0) begin_integration(a, b, n, (f(a) + f(b)) / 3.);
                    else begin_integration(a, b, n + 1);  
                    for (uint32_t i{1}; i < steps_; i++) sum_ += 2 * (1 + i % 2) * (f(a_ + i * h_)) / 3.; 
                    integral_ = sum_ * h_; 

                }


                void simpson_fixed(const double& a, 
                                             const double& b, 
                                             const std::function<double(double)>& f, 
                                             const double& prec = 1.e-6) {

                    double old_integral_2{}, old_integral_3{};
                    begin_integration(a, b, 2, (f(a) + f(b)) / 3.); 
                    while (true) {
                        old_integral_3 = old_integral_2;
                        old_integral_2 = old_integral_;
                        old_integral_ = integral_; 
                        simpson(a_, b_, f, steps_ * 2);
                        error_ = 256 * std::fabs(1024 * integral_ - 1104 * old_integral_ + 81 * old_integral_2 - old_integral_3) / 240975;
                        if (error_ < prec) break; 
                    }
                    integral_ = (1024 * integral_ - 80 * old_integral_ + old_integral_2) / 945; 
                    
                }


                void mean(const double& a, 
                                    const double& b, 
                                    const std::function<double(double)>& f, 
                                    const uint32_t& n = 1000) {

                    begin_integration(a, b, n); 
                    for (uint32_t i{}; i < n; i ++) sum_ += f(rg_.unif(a, b)); 
                    integral_ = (b_ - a_) * sum_ / steps_; 

                }


                void mean_fixed(const double& a, 
                                const double& b, 
                                const std::function<double(double)>& f, 
                                const double& prec = 1.e-6) {

                    std::vector<double> k{};
                    for (unsigned i{}; i < 10000; i++) {
                        mean(a, b, f);
                        k.emplace_back(integral_); 
                    }
                    double k_mean = sqrt(100) * descriptive_statistics::sd(k); 
                    mean(a, b, f, static_cast<uint>(std::pow(k_mean / prec, 2))); 

                }

        
                constexpr void hit_or_miss(const double& a, 
                                           const double& b, 
                                           const std::function<double(double)>& f, 
                                           const double& fmax, 
                                           const uint32_t& n = 1000) {

                    begin_integration(a, b, n); 
                    double x{}, y{}; 
                    uint32_t hits{};
                    for (uint32_t i{}; i < n; i ++) {
                        x = rg_.unif(a, b); 
                        y = rg_.unif(0., fmax);  
                        if (y <= f(x)) hits++; 
                    }
                    integral_ = (b_ - a_) * fmax * hits / n; 

                }


                void hit_or_miss_fixed(const double& a, 
                                       const double& b, 
                                       const std::function<double(double)>& f, 
                                       const double& fmax, 
                                       const double& prec = 1.e-6) {

                    std::vector<double> k{};
                    for (unsigned i{}; i < 10000; i++) {
                        hit_or_miss(a, b, f, fmax);
                        k.emplace_back(integral_); 
                    }
                    double k_mean = sqrt(100) * descriptive_statistics::sd(k); 
                    uint32_t N = static_cast<uint32_t>(std::pow(k_mean / prec, 2)); 
                    hit_or_miss(a, b, f, fmax, N); 

                }


        }; // class integral


    } // namespace tools


} // namespace math
