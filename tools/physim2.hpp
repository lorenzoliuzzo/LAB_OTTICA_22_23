
/**
 * @brief   PHYSIM is a WIP c++ header file that constains the basic tools for computational physics
 * 
 * @author: Lorenzo Liuzzo
 * 
 * @email:  lorenzoliuzzo@outlook.com
 * 
 * @github: https://github.com/lorenzoliuzzo/physim 
 * 
 * @updated: 29/11/2022
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


using scalar = double;


/// @brief Namespace contains numerical constants, algorithms, data structures and tools for maths
namespace math {
            

    /// @brief Namespace contains some usefull tools for working with double precision floating point numbers
    namespace tools {


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

            if (v1 == 0.0 || std::fpclassify(v1) == FP_SUBNORMAL) 
                return true; 

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


    } // namespace tools


    /// @brief Namespace contains some math constants
    namespace constants {

        
        constexpr scalar infinity = std::numeric_limits<scalar>::infinity();
        
        constexpr scalar invalid_conversion = std::numeric_limits<scalar>::signaling_NaN();

        constexpr scalar pi = 3.14159265358979323846;
        
        constexpr scalar e = 2.7182818284590452353603;

        constexpr scalar sqrt2 = std::sqrt(2.0); 

        constexpr scalar sqrt3 = std::sqrt(3.0);


    } // namespace constants
    

} // namespace math


/// @brief Namespace contains physical constants, data structures and tools for computational physics
namespace physics {


    /// @brief Namespace contains some usefull tools for working with physical measurements
    namespace measurements {


        /// @brief Namespace contains the necessary tools for working with units of measurement
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


            /// @brief Class represents an unit base using powers of the seven SI unit bases 
            class unit_base {


                public:
                
                // =============================================
                // constructor and destructor
                // =============================================  

                    /// @brief Construct a new default unit_base object
                    explicit constexpr unit_base() noexcept :

                        metre_{0}, second_{0}, kilogram_{0}, 
                        ampere_{0}, kelvin_{0}, mole_{0}, candela_{0} {}


                    /**
                     * @brief Construct a new unit_base object with powers of the seven SI unit bases
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
                     * @brief Construct a new unit_base object from a string
                     * 
                     * @param unit_string: string represents the unit base
                     */
                    constexpr unit_base(const std::string& unit_string) noexcept : 
                        
                        metre_{0}, second_{0}, kilogram_{0},
                        ampere_{0}, kelvin_{0}, mole_{0}, candela_{0} {
                        
                        if (!unit_string.empty()) {

                            size_t finder = unit_string.find("m");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') 
                                    metre_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') 
                                        metre_ = - std::stoi(unit_string.substr(finder + 3));
                                    else 
                                        metre_ = std::stoi(unit_string.substr(finder + 2));
                                }
                                
                            }

                            finder = unit_string.find("s");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') 
                                    second_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') 
                                        second_ = - std::stoi(unit_string.substr(finder + 3));
                                    else 
                                        second_ = std::stoi(unit_string.substr(finder + 2));
                                }
                                
                            }

                            finder = unit_string.find("kg");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') 
                                    kilogram_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') 
                                        kilogram_ = - std::stoi(unit_string.substr(finder + 3));
                                    else 
                                        kilogram_ = std::stoi(unit_string.substr(finder + 2));
                                }

                            }

                            finder = unit_string.find("A");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') 
                                    ampere_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') 
                                        ampere_ = - std::stoi(unit_string.substr(finder + 3));
                                    else 
                                        ampere_ = std::stoi(unit_string.substr(finder + 2));
                                }                            
                            }

                            finder = unit_string.find("K");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') 
                                    kelvin_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') 
                                        kelvin_ = - std::stoi(unit_string.substr(finder + 3));
                                    else 
                                        kelvin_ = std::stoi(unit_string.substr(finder + 2));
                                }

                            }

                            finder = unit_string.find("mol");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') 
                                    mole_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') 
                                        mole_ = - std::stoi(unit_string.substr(finder + 3));
                                    else 
                                        mole_ = std::stoi(unit_string.substr(finder + 2));
                                }

                            }

                            finder = unit_string.find("cd");
                            if (finder != std::string::npos) {

                                if (finder == unit_string.size() - 1 || unit_string.at(finder + 1) != '^') 
                                    candela_ = 1; 
                                else {
                                    if (unit_string.at(finder + 2) == '-') 
                                        candela_ = - std::stoi(unit_string.substr(finder + 3));
                                    else 
                                        candela_ = std::stoi(unit_string.substr(finder + 2));
                                }
                                
                            }

                        } 

                    }


                    /**
                     * @brief Copy construct a new unit_base from an other unit base object
                     * 
                     * @param other: unit_base object to copy as l-value const reference
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
                     * @param other: unit_base to copy as l-value const reference
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
                     * @param other: unit_base object to multiply with as l-value const reference
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
                     * @param other: unit_base object to divide with as l-value const reference
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
                     * @param other: unit_base object to multiply with as l-value const reference
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
                     * @param other: unit_base object to divide with as l-value const reference
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
                     * @param other: unit_base object to compare with as l-value const reference
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
                     * @param other: unit_base object to !compare with as l-value const reference
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

                        if (this->has_valid_root(power)) 
                            return unit_base(metre_ / power,
                                             second_ / power,
                                             kilogram_ / power,
                                             ampere_ / power,
                                             kelvin_ / power,
                                             mole_ / power,
                                             candela_ / power);

                        else throw std::invalid_argument("Invalid root power");

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
                        
                        if (metre_ == 1) 
                            unit_base_string.append("m");
                        else if (metre_ != 0) 
                            unit_base_string.append("m^" + std::to_string(metre_)); 

                        if (second_ == 1) 
                            unit_base_string.append("s"); 
                        else if (second_ != 0) 
                            unit_base_string.append("s^" + std::to_string(second_)); 

                        if (kilogram_ == 1) 
                            unit_base_string.append("kg"); 
                        else if (kilogram_ != 0)    
                            unit_base_string.append("kg^" + std::to_string(kilogram_)); 

                        if (ampere_ == 1) 
                            unit_base_string.append("A"); 
                        else if (ampere_ != 0)  
                            unit_base_string.append("A^" + std::to_string(ampere_)); 

                        if (kelvin_ == 1) 
                            unit_base_string.append("K");
                        else if (kelvin_ != 0) 
                            unit_base_string.append("K^" + std::to_string(kelvin_)); 

                        if (mole_ == 1) 
                            unit_base_string.append("mol"); 
                        else if (mole_ != 0) 
                            unit_base_string.append("mol^" + std::to_string(mole_)); 

                        if (candela_ == 1) 
                            unit_base_string.append("cd"); 
                        else if (candela_ != 0) 
                            unit_base_string.append("cd^" + std::to_string(candela_)); 

                        return unit_base_string; 
                    
                    }


                    /// @brief Print the unit_base to the standard output
                    inline void print() const noexcept {

                        std::cout << this->to_string(); 

                    }


                private: 

                // =============================================
                // class members
                // =============================================


                    signed int metre_ : bitwidth::metre; ///< Metre exponent
                    
                    signed int second_ : bitwidth::second; ///< Second exponent
                    
                    signed int kilogram_ : bitwidth::kilogram; ///< Kilogram exponent
                    
                    signed int ampere_ : bitwidth::ampere; ///< Ampere exponent
                    
                    signed int kelvin_ : bitwidth::kelvin; ///< Kelvin exponent
                    
                    signed int mole_ : bitwidth::mole; ///< Mole exponent
                    
                    signed int candela_ : bitwidth::candela; ///< Candela exponent
        

                    /// @brief static array with the SI exponents 
                    static constexpr uint32_t bits[7] = { bitwidth::metre, bitwidth::second, bitwidth::kilogram, bitwidth::ampere, bitwidth::kelvin, bitwidth::mole, bitwidth::candela };  


            }; // class unit_base


            /// @brief Class represents an unit prefix using a multiplier (double) and a symbol (char)
            class unit_prefix {


                public: 

                // =============================================
                // constructor and destructor
                // ============================================= 

                    /// @brief Construct a new default unit prefix object
                    constexpr unit_prefix() noexcept : 

                        multiplier_{1.}, 
                        symbol_{'\0'} {}


                    /**
                     * @brief Create a new unit_prefix object from a multiplier and a symbol
                     * 
                     * @param mult: scalar multiplier for scaling the measurement 
                     * @param symbol: char symbol for the string reppresentation
                     * 
                     * @note mult must be positive (> 0)
                     */
                    constexpr unit_prefix(const scalar& mult, 
                                          const char& symbol) : 

                        multiplier_{mult}, 
                        symbol_{symbol} {

                            if (mult <= 0) 
                                throw std::invalid_argument("unit_prefix multiplier must be positive");

                        }
                

                    /// @brief Default destructor
                    ~unit_prefix() = default; 


                // =============================================
                // operators
                // ============================================= 

                    /**
                     * @brief Copy assignment operator
                     * 
                     * @param other: unit_prefix to copy as l-value const reference
                     * 
                     * @return constexpr unit_prefix&
                     */
                    constexpr unit_prefix& operator=(const unit_prefix& other) noexcept {

                        multiplier_ = other.multiplier_; 
                        symbol_ = other.symbol_; 

                        return *this; 

                    }


                    /**
                     * @brief Multiply this unit_prefix to an unit_prefix by multiplying the multipliers together
                     * 
                     * @param other: unit_prefix object to multiply with as l-value const reference
                     * 
                     * @return constexpr unit_prefix&
                     */
                    constexpr unit_prefix& operator*=(const unit_prefix& other) noexcept {

                        multiplier_ *= other.multiplier_; 

                        return *this; 

                    }


                    /**
                     * @brief Divide this unit_prefix to an unit_prefix by dividing the multipliers together
                     * 
                     * @param other: unit_prefix object to divide with as l-value const reference
                     * 
                     * @return constexpr unit_prefix&
                     */
                    constexpr unit_prefix& operator/=(const unit_prefix& other) noexcept {

                        multiplier_ /= other.multiplier_; 

                        return *this; 

                    }


                    /**
                     * @brief Perform a multiplication between unit_prefixes 
                     * 
                     * @param other: unit_prefix object to multiply with as l-value const reference
                     * 
                     * @return constexpr unit_prefix 
                     */
                    constexpr unit_prefix operator*(const unit_prefix& other) const noexcept {

                        return unit_prefix(multiplier_ * other.multiplier_, symbol_); 

                    }


                    /**
                     * @brief Perform a division between unit_prefixes 
                     * 
                     * @param other: unit_prefix object to muldividetiply with as l-value const reference
                     * 
                     * @return constexpr unit_prefix 
                     */
                    constexpr unit_prefix operator/(const unit_prefix& other) const noexcept {

                        return unit_prefix(multiplier_ / other.multiplier_, symbol_); 

                    }


                    /**
                     * @brief Equality operator
                     * 
                     * @param other: unit_prefix object to compare with as l-value const reference
                     * 
                     * @return bool
                     */
                    constexpr bool operator==(const unit_prefix& other) const noexcept {
                    
                        return (multiplier_ == other.multiplier_ && symbol_ == other.symbol_);   

                    }


                    /**
                     * @brief Inequality operator
                     * 
                     * @param other: unit_prefix object to !compare with as l-value const reference
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
                     */
                    constexpr unit_prefix pow(const int& power) const noexcept { 
                        
                        return unit_prefix(std::pow(multiplier_, power), symbol_); 
                        
                    }


                    /**
                     * @brief Take the square of the unit_prefix
                     * 
                     * @return constexpr unit_prefix 
                     */
                    constexpr unit_prefix square() const noexcept { 
                        
                        return unit_prefix(std::pow(multiplier_, 2), symbol_); 
                        
                    }


                    /**
                     * @brief Take the cube of the unit_prefix
                     * 
                     * @return constexpr unit_prefix 
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
                     */
                    constexpr unit_prefix root(const int& power) const noexcept { 
                        
                        return unit_prefix(std::pow(multiplier_, power), symbol_); 
                        
                    }


                    /**
                     * @brief Take the square root of the unit_prefix
                     * 
                     * @return constexpr unit_prefix 
                     */
                    constexpr unit_prefix sqrt() const noexcept { 
                        
                        return unit_prefix(std::pow(multiplier_, 1. / 2.), symbol_); 
                        
                    }


                    /**
                     * @brief Take the cube root of the unit_prefix
                     * 
                     * @return constexpr unit_prefix 
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
                     * @return constexpr scalar 
                     */
                    constexpr scalar multiplier() const noexcept {

                        return multiplier_; 

                    }


                    /**
                     * @brief Get the multiplier of the unit_prefix
                     * 
                     * @return constexpr scalar& 
                     */
                    constexpr scalar& multiplier() noexcept {

                        return multiplier_; 

                    }


                    /**
                     * @brief Get the symbol of the unit_prefix
                     * 
                     * @return constexpr char
                     */
                    constexpr char symbol() const noexcept {

                        return symbol_; 

                    }


                    /**
                     * @brief Get the symbol of the unit_prefix
                     * 
                     * @return constexpr char&
                     */
                    constexpr char& symbol() noexcept {

                        return symbol_; 

                    }


                private:    

                    scalar multiplier_; ///< multiplier of the unit_prefix

                    char symbol_; ///< symbol of the unit_prefix


            }; // class unit_prefix


            /// @brief Class represents an unit of measurement as an unit_base and an unit_prefix
            class unit {


                public:

                // =============================================
                // constructors & destructor
                // ============================================= 

                    /// @brief Construct a new default unit object
                    explicit constexpr unit() noexcept : 
                        
                        base_(), 
                        prefix_() {}


                    /**
                     * @brief Construct a new unit object from a prefix and a base 
                     * 
                     * @param prefix: unit_prefix as l-value const reference
                     * @param base: unit_base as l-value const reference
                     */
                    explicit constexpr unit(const unit_prefix& prefix, 
                                            const unit_base& base) noexcept : 
                        
                        base_(base), 
                        prefix_(prefix) {}


                    /**
                     * @brief Construct a new unit object from a unit_prefix and an unit_base
                     * 
                     * @param prefix: unit_prefix as r-value reference
                     * @param base: unit_base as r-value reference
                     */
                    explicit constexpr unit(unit_prefix&& prefix, 
                                            unit_base&& base) noexcept : 
                        
                        base_(std::move(base)), 
                        prefix_(std::move(prefix)) {}


                    /**
                     * @brief Construct a new unit object from an unit_base and an unit_prefix
                     * 
                     * @param base: unit_base object as l-value const reference
                     * @param prefix: unit_prefix as l-value const reference
                     */
                    explicit constexpr unit(const unit_base& base, 
                                            const unit_prefix& prefix = unit_prefix()) noexcept : 
                        
                        base_(base), 
                        prefix_(prefix) {}


                    /**
                     * @brief Construct a new unit object from a unit_base and unit_prefix
                     * 
                     * @param base: unit_base as r-value reference
                     * @param prefix: unit_prefix as r-value reference
                     */
                    explicit constexpr unit(unit_base&& base, 
                                            unit_prefix&& prefix) noexcept : 
                        
                        base_(std::move(base)), 
                        prefix_(std::move(prefix)) {}


                    /**
                     * @brief Construct a new unit from an unit_prefix and an unit object
                     * 
                     * @param prefix: unit_prefix as l-value const reference
                     * @param unit: unit as l-value const reference
                     */
                    explicit constexpr unit(const unit_prefix& prefix, 
                                            const unit& unit) noexcept : 
                        
                        base_(unit.base_), 
                        prefix_(prefix * unit.prefix_) {}


                    /**
                     * @brief Copy construct a new unit from another unit object
                     * 
                     * @param other: unit object to copy from as l-value const reference
                     */
                    constexpr unit(const unit& other) noexcept : 

                        base_(other.base_), 
                        prefix_(other.prefix_) {}


                    /**
                     * @brief Move construct a new unit from another unit object
                     * 
                     * @param other: unit object to move from as r-value reference
                     */
                    constexpr unit(unit&& other) noexcept : 

                        base_(std::move(other.base_)), 
                        prefix_(std::move(other.prefix_)) {}


                    /// @brief Default destructor
                    ~unit() = default; 


                // =============================================
                // operators
                // ============================================= 

                    /**
                     * @brief Copy assignment operator
                     * 
                     * @param other: unit object to copy from as l-value const reference
                     * 
                     * @return constexpr unit&
                     */
                    constexpr unit& operator=(const unit& other) noexcept {

                        base_ = other.base_;
                        prefix_ = other.prefix_;

                        return *this;

                    }


                    /**
                     * @brief Move assignment operator
                     * 
                     * @param other: unit object to move from as r-value reference
                     * 
                     * @return constexpr unit& 
                     */
                    constexpr unit& operator=(unit&& other) noexcept {

                        base_ = std::move(other.base_);
                        prefix_ = std::move(other.prefix_);

                        return *this;

                    }


                    /**
                     * @brief Multiply this unit with an unit
                     * 
                     * @param other: unit object to multiply with as l-value const reference
                     *                      
                     * @return constexpr unit&
                     */
                    constexpr unit& operator*=(const unit& other) noexcept { 
                        
                        base_ *= other.base_; 
                        prefix_ *= other.prefix_; 

                        return *this; 
                        
                    }


                    /**
                     * @brief Multiply this unit with an unit
                     * 
                     * @param other: unit object to multiply with as r-value reference
                     *                      
                     * @return constexpr unit&
                     */
                    constexpr unit& operator*=(unit&& other) noexcept { 
                        
                        base_ *= std::move(other.base_); 
                        prefix_ *= std::move(other.prefix_); 

                        return *this; 
                        
                    }


                    /**
                     * @brief Divide this unit with an unit
                     * 
                     * @param other: unit object to divide with as l-value const reference
                     *                      
                     * @return constexpr unit&
                     */
                    constexpr unit& operator/=(const unit& other) noexcept { 
                        
                        base_ /= other.base_; 
                        prefix_ /= other.prefix_; 

                        return *this; 
                        
                    }                 


                    /**
                     * @brief Divide this unit with an unit
                     * 
                     * @param other: unit object to divide with as r-value reference
                     *                      
                     * @return constexpr unit&
                     */
                    constexpr unit& operator/=(unit&& other) noexcept { 
                        
                        base_ /= std::move(other.base_); 
                        prefix_ /= std::move(other.prefix_); 

                        return *this; 
                        
                    }  


                    /**
                     * @brief Perform a multiplication by multiply the bases and the prefixes
                     * 
                     * @param other: unit object to multiply with as l-value const reference
                     *                      
                     * @return constexpr unit
                     */
                    constexpr unit operator*(const unit& other) const noexcept { 
                        
                        return unit(base_ * other.base_, prefix_ * other.prefix_); 
                        
                    }


                    /**
                     * @brief Perform a division by divide the bases and the prefixes
                     * 
                     * @param other: unit object to divide with as l-value const reference
                     *                      
                     * @return constexpr unit
                     */
                    constexpr unit operator/(const unit& other) const noexcept { 
                        
                        return unit(base_ / other.base_, prefix_ / other.prefix_);
                        
                    }                    


                    /**
                     * @brief Equality operator
                     * 
                     * @param other: unit object to compare with as l-value const reference
                     * 
                     * @return constexpr bool 
                     */
                    constexpr bool operator==(const unit& other) const noexcept {

                        if (base_ != other.base_) 
                            return false; 
                        else 
                            return prefix_ == other.prefix_; 

                    }


                    /**
                     * @brief Inequality operator
                     * 
                     * @param other: unit object to !compare with as l-value const reference
                     * 
                     * @return constexpr bool 
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
                     */
                    constexpr unit pow(const int& power) const noexcept { 
                        
                        return unit(base_.pow(power), prefix_.pow(power)); 
                        
                    }


                    /**
                     * @brief Take the square of the unit
                     * 
                     * @return constexpr unit 
                     */
                    constexpr unit square() const noexcept { 
                        
                        return unit(base_.square(), prefix_.square()); 
                        
                    }


                    /**
                     * @brief Take the cube of the unit
                     * 
                     * @return constexpr unit 
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
                     */
                    constexpr unit root(const int& power) const { 
                        
                        return unit(base_.root(power), prefix_.root(power)); 
                        
                    }


                    /**
                     * @brief Take the square root of the unit
                     * 
                     * @return constexpr unit 
                     */
                    constexpr unit sqrt() const { 
                        
                        return unit(base_.sqrt(), prefix_.sqrt()); 
                        
                    }


                    /**
                     * @brief Take the cube root of the unit
                     * 
                     * @return constexpr unit 
                     */
                    constexpr unit cbrt() const { 
                        
                        return unit(base_.cbrt(), prefix_.cbrt()); 
                        
                    }


                // =============================================
                // setters & getters
                // ============================================= 

                    /**
                     * @brief Get the unit base object
                     * 
                     * @return constexpr unit_base 
                     */
                    constexpr unit_base base() const noexcept { 
                        
                        return base_; 
                        
                    }


                    /**
                     * @brief Get the unit base object
                     * 
                     * @return constexpr unit_base& 
                     */
                    constexpr unit_base& base() noexcept { 
                        
                        return base_; 
                        
                    }


                    /**
                     * @brief Get the unit prefix object
                     * 
                     * @return constexpr unit_prefix 
                     */
                    constexpr unit_prefix prefix() const noexcept { 
                        
                        return prefix_; 
                        
                    }


                    /**
                     * @brief Get the unit prefix object
                     * 
                     * @return constexpr unit_prefix&
                     */
                    constexpr unit_prefix& prefix() noexcept { 
                        
                        return prefix_; 
                        
                    }


                    /**
                     * @brief Get the unit object
                     * 
                     * @return constexpr unit
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
                     * @return constexpr scalar 
                     * @note the units will only convert if they have the same base unit
                     */
                    constexpr scalar convertion_factor(const unit& other) const noexcept { 
                        
                        return (base_ == other.base_) ? prefix_.multiplier() / other.prefix_.multiplier() : math::constants::invalid_conversion; 
                        
                    }


                    /**
                     * @brief Convert a value in the current unit to the target unit 
                     * 
                     * @param value: value to convert 
                     * @param other: desired unit unit
                     * 
                     * @return constexpr scalar 
                     * @note the units will only convert if they have the same base unit
                     */
                    constexpr scalar convert(const scalar& value, const unit& other) const noexcept { 
                        
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


                    /// @brief Print the unit to the standard output
                    inline void print() const noexcept {

                        std::cout << this->to_string();

                    }


                private:
            
                // =============================================
                // class members
                // ============================================= 

                    unit_base base_; ///< unit base

                    unit_prefix prefix_; ///< unit prefix


            }; // class unit     


            #define CREATE_UNIT_BASE(NAME, METRE, SECOND, KILOGRAM, AMPERE, KELVIN, MOLE, CANDELA) \
                constexpr unit_base NAME(METRE, SECOND, KILOGRAM, AMPERE, KELVIN, MOLE, CANDELA); 
            

            #define CREATE_UNIT_PREFIX(NAME, MULTIPLIER, SYMBOL) \
                constexpr unit_prefix NAME(MULTIPLIER, SYMBOL); 


            #define CREATE_UNIT(NAME, PREFIX, BASE) \
                constexpr unit NAME(PREFIX, BASE); 


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
                constexpr unit K(prefix::default_type, base::Kelvin);
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

                constexpr unit newton(kg * m / s.square());
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

                /// @brief Construct a new default measurement object (value = 0, unit = unitless)
                explicit constexpr measurement() noexcept : 

                    value_{0.0}, 
                    units_() {}


                /**
                 * @brief Construct a new measurement object from a value and an unit
                 * 
                 * @param value: scalar value of the measurement
                 * @param unit: unit of measurement as l-value const reference
                 */
                explicit constexpr measurement(const scalar& value, 
                                               const unit& units = unit()) noexcept : 
                                               
                    value_{value}, 
                    units_(units) {}


                /**
                 * @brief Copy constuct a new measurement object from another measurement
                 * 
                 * @param other: the measurement to copy
                 */
                constexpr measurement(const measurement& other) noexcept : 
                
                    value_{other.value_}, 
                    units_(other.units_) {}
                    
                
                /**
                 * @brief Move construct a new measurement object from another measurement
                 * 
                 * @param other: the measurement to move from
                 */
                constexpr measurement(measurement&& other) noexcept :
                    
                    value_{std::move(other.value_)}, 
                    units_(std::move(other.units_)) {}


                /// @brief Default destructor of the measurement object
                ~measurement() = default;


            // =============================================                                                                                         
            // operators
            // =============================================  
                
                /**
                 * @brief Copy assign a measurement from another measurement
                 * 
                 * @param other: the measurement to copy as l-value const reference
                 * 
                 * @return measurement& 
                 */
                constexpr measurement& operator=(const measurement& other) noexcept {

                    value_ = other.value_;
                    units_ = other.units_;

                    return *this;

                }     


                /**
                 * @brief Add a measurement to this measurement
                 * 
                 * @param other: the measurement to add as l-value const reference
                 * 
                 * @return measurement& 
                 */
                constexpr measurement& operator+=(const measurement& other) { 
                    
                    if (units_.base() != other.units_.base()) 
                        throw std::invalid_argument("Cannot add measurements with different unit bases");
                    
                    if (units_ != unitless) 
                        value_ += other.value_as(units_); 
                    else {
                        value_ += other.value_; 
                        units_ = other.units_; 
                    }  

                    return *this; 
                
                }


                /**
                 * @brief Subtract a measurement to this measurement
                 * 
                 * @param other: the measurement to subtract as l-value const reference
                 * 
                 * @return measurement& 
                 */
                constexpr measurement& operator-=(const measurement& other) { 

                    if (units_.base() != other.units_.base()) 
                        throw std::invalid_argument("Cannot subtract measurements with different unit bases");
                    
                    if (units_ != unitless) 
                        value_ -= other.value_as(units_);   
                    else {
                        value_ -= other.value_;    
                        units_ = other.units_; 
                    }   

                    return *this; 
                
                }


                /**
                 * @brief Multiply this measurement and a measurement
                 * 
                 * @param other: the measurement to multiply as l-value const reference
                 * 
                 * @return measurement& 
                 */
                constexpr measurement& operator*=(const measurement& meas) noexcept { 

                    value_ *= meas.value_;
                    units_ *= meas.units_;    

                    return *this; 
                
                }


                /**
                 * @brief Divide this measurement and a measurement
                 * 
                 * @param other: the measurement to divide as l-value const reference
                 * 
                 * @return measurement& 
                 */
                constexpr measurement& operator/=(const measurement& meas) { 

                    if (meas.value_ == 0.0)
                        throw std::runtime_error("Cannot divide measurement by a zero measurement");
                   
                    value_ /= meas.value_;
                    units_ /= meas.units_;                                    
                    
                    return *this; 
                
                }


                /**
                 * @brief Move assign a measurement from another measurement
                 * 
                 * @param other: the measurement to move from as r-value reference
                 * 
                 * @return measurement& 
                 */
                constexpr measurement& operator=(measurement&& other) noexcept {

                    value_ = std::move(other.value_);
                    units_ = std::move(other.units_);

                    return *this;

                }


                /**
                 * @brief Add a measurement to this measurement
                 * 
                 * @param other: the measurement to add as r-value reference
                 * 
                 * @return measurement& 
                 */
                constexpr measurement& operator+=(measurement&& other) { 

                    if (units_.base() != other.units_.base()) 
                        throw std::invalid_argument("Cannot add measurements with different unit bases");
                    
                    if (units_ != unitless) 
                        value_ += std::move(other.value_as(units_));   
                    else {
                        value_ += std::move(other.value_);  
                        units_ = std::move(other.units_); 
                    }   

                    return *this; 
                
                }


                /**
                 * @brief Subtract a measurement to this measurement
                 * 
                 * @param other: the measurement to subtract as r-value reference
                 * 
                 * @return measurement& 
                 */
                constexpr measurement& operator-=(measurement&& other) { 

                    if (units_.base() != other.units_.base()) 
                        throw std::invalid_argument("Cannot subtract measurements with different unit bases");
                   
                    if (units_ != unitless)
                        value_ -= std::move(other.value_as(units_));   
                    else {
                        value_ -= std::move(other.value_);    
                        units_ = std::move(other.units_); 
                    }     

                    return *this;   

                }


                /**
                 * @brief Multiply this measurement and a measurement
                 * 
                 * @param other: the measurement to multiply as r-value reference
                 * 
                 * @return measurement& 
                 */
                constexpr measurement& operator*=(measurement&& meas) noexcept { 

                    value_ *= std::move(meas.value_);
                    units_ *= std::move(meas.units_);   

                    return *this; 
                
                }


                /**
                 * @brief Divide this measurement and a measurement
                 * 
                 * @param other: the measurement to divide as r-value reference
                 * 
                 * @return measurement& 
                 */
                constexpr measurement& operator/=(measurement&& meas) { 
                    
                    if (meas.value_ == 0.0) 
                        throw std::runtime_error("Cannot divide measurement by a zero measurement");
                   
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
                 */
                constexpr measurement& operator*=(const scalar& value) noexcept { 

                    value_ *= value;    

                    return *this; 
                
                }


                /**
                 * @brief Divide this measurement and a scalar
                 * 
                 * @param value 
                 * 
                 * @return constexpr measurement& 
                 */
                constexpr measurement& operator/=(const scalar& value) { 
                    
                    if (value == 0.0) 
                        throw std::runtime_error("Cannot divide measurement by 0");
                    
                    value_ /= value;                    
                    
                    return *this; 
                
                }

                
                /**
                 * @brief Return the opposite of this measurement
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement operator-() const noexcept { 
                    
                    return measurement(-value_, units_); 
                
                }


                /**
                 * @brief Sum a measurement to this measurement
                 * 
                 * @param other: measurement to add as l-value const reference
                 * 
                 * @return measurement 
                 */
                constexpr measurement operator+(const measurement& other) const { 
                    
                    if (units_.base() != other.units_.base()) 
                        throw std::invalid_argument("Cannot sum measurements with different unit bases");
                   
                    return measurement(value_ + other.value_as(units_), units_); 
                
                }


                /**
                 * @brief Subtract a measurement to this measurement
                 * 
                 * @param other: measurement to subtract as l-value const reference
                 *  
                 * @return constexpr measurement 
                 */
                constexpr measurement operator-(const measurement& other) const { 

                    if (units_.base() != other.units_.base()) 
                        throw std::invalid_argument("Cannot subtract measurements with different unit bases");
                    
                    return measurement(value_ - other.value_as(units_), units_); 
                
                }

                
                /**
                 * @brief Multiply this measurement and a measurement
                 * 
                 * @param other: measurement to multiply as l-value const reference
                 * 
                 * @return constexpr measurement 
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
                 */
                constexpr measurement operator/(const measurement& other) const { 
                    
                    if (other.value_ == 0.0) 
                        throw std::runtime_error("Cannot divide measurement by a zero measurement");
                   
                    return measurement(value_ / other.value_, units_ / other.units_); 
                
                }

                
                /**
                 * @brief Multiply this measurement and a scalar
                 * 
                 * @param value 
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement operator*(const scalar& value) const noexcept { 
                    
                    return measurement(value_ * value, units_); 
                    
                }

                
                /**
                 * @brief Divide this measurement and a scalar
                 * 
                 * @param value
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement operator/(const scalar& value) const { 
                    
                    if (value_ == 0.0) 
                        throw std::runtime_error("Cannot divide measurement by 0");
                  
                    return measurement(value_ / value, units_); 
                    
                } 


                /**
                 * @brief Perform a product between a scalar and a measurement
                 * 
                 * @param val: scalar
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement 
                 */
                friend constexpr measurement operator*(const scalar& val, 
                                                       const measurement& meas) noexcept { 
                                                    
                    return meas * val; 
                    
                }
                

                /**
                 * @brief Perform a division between a scalar and a measurement
                 * 
                 * @param val: scalar
                 * @param meas: measurement
                 * 
                 * @return constexpr measurement 
                 */
                friend constexpr measurement operator/(const scalar& val, 
                                                       const measurement& meas) { 

                    if (meas.value_ == 0.0) 
                        throw std::runtime_error("Cannot divide a scalar by a zero measurement");          
                    
                    return measurement(val / meas.value_, meas.units_.inv()); 
                    
                }


                /**
                 * @brief Equality operator
                 * 
                 * @param other: measurement to compare as l-value const reference
                 * 
                 * @return bool
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
                constexpr bool operator==(const scalar& val) const { 
                    
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
                constexpr bool operator!=(const scalar& val) const { 
                    
                    return !operator==(val); 
                    
                }
                
               
                /**
                 * @brief More than with value
                 * 
                 * @param val 
                 * 
                 * @return bool
                 */
                constexpr bool operator>(const scalar& val) const { 
                    
                    return value_ > val; 
                    
                }
               
               
                /**
                 * @brief Less than with value
                 * 
                 * @param val 
                 * 
                 * @return bool
                 */
                constexpr bool operator<(const scalar& val) const { 
                    
                    return value_ < val; 
                    
                }
                
                
                /**
                 * @brief More than or equal with value
                 * 
                 * @param val 
                 * 
                 * @return bool
                 */
                constexpr bool operator>=(const scalar& val) const { 
                    
                    return (value_ >= val) ? true : operator==(val); 
                    
                }
                
                
                /**
                 * @brief Less than or equal with value
                 * 
                 * @param val 
                 * 
                 * @return bool
                 */
                constexpr bool operator<=(const scalar& val) const { 
                    
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
                 * @brief Invert the measurement
                 * 
                 * @return constexpr measurement 
                 * 
                 * @note Cannot invert a measurement with a zero value
                 */
                constexpr measurement inv() const { 
                    
                    if (value_ == 0) 
                        throw std::runtime_error("Cannot invert a measurement with a zero value");
                        
                    return measurement(1 / value_, units_.inv()); 
                
                }
                
                
                /**
                 * @brief Get the absolute measurement object
                 * 
                 * @param meas: measurement as l-value const reference
                 * @return constexpr measurement
                 */
                friend constexpr measurement abs(const measurement& meas) noexcept { 
                    
                    return (meas.value_ < 0.0) ? -meas : meas; 
                
                }

                
                /**
                 * @brief Take the power of the measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * @param power 
                 * 
                 * @return constexpr measurement 
                 */
                friend constexpr measurement pow(const measurement& meas, const int& power) noexcept { 
                    
                    return measurement(std::pow(meas.value_, power), meas.units_.pow(power)); 
                
                }


                /**
                 * @brief Take the root power of the measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * @param power 
                 * 
                 * @return constexpr measurement 
                 */
                friend constexpr measurement root(const measurement& meas, const int& power) { 
                    
                    return measurement(std::pow(meas.value_, 1.0 / power), meas.units_.root(power)); 
                
                }


                /**
                 * @brief Take the square of the measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement 
                 */
                friend constexpr measurement square(const measurement& meas) noexcept { 
                    
                    return measurement(std::pow(meas.value_, 2), meas.units_.square()); 
                
                }

                
                /**
                 * @brief Take the cube of the measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement 
                 */
                friend constexpr measurement cube(const measurement& meas) noexcept { 
                    
                    return measurement(std::pow(meas.value_, 3), meas.units_.cube()); 
                
                }

                
                /**
                 * @brief Take the square root of the measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement sqrt(const measurement& meas) { 
                    
                    if (meas.value_ < 0.0) 
                        throw std::runtime_error("Cannot take the square root of a negative measurement");
                    
                    return measurement(std::sqrt(meas.value_), meas.units_.sqrt()); 
                
                }

                
                /**
                 * @brief Take the cubic root of the measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement cbrt(const measurement& meas) { 
                    
                    return measurement(std::cbrt(meas.value_), meas.units_.cbrt()); 
                
                }


                /**
                 * @brief Take the exponential of a measurement
                 * @note The base of the exponential is e
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement exp(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the exponential of a measurement that is not unitless"); 
                    
                    return measurement(std::exp(meas.value_), unitless); 
                
                }


                /**
                 * @brief Take the logarithm of a measurement
                 * @note The base of the logarithm is e
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement log(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the logarithm of a measurement that is not unitless"); 
                    
                    return measurement(std::log(meas.value_), unitless); 
                
                }
                

                /**
                 * @brief Take the exponential base 10 of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement exp10(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the exponential of a measurement that is not unitless"); 
                    
                    return measurement(std::pow(10, meas.value_), unitless); 
                
                }


                /**
                 * @brief Take the logarithm base 10 of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement log10(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the logarithm of a measurement that is not unitless"); 
                    
                    return measurement(std::log10(meas.value_), unitless); 
                
                }


                /**
                 * @brief Take the sine of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement sin(const measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the sine of a measurement that is not in radians"); 
                    
                    return measurement(std::sin(meas.value_), unitless); 
                
                }


                /**
                 * @brief Take the cosine of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement cos(const measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the cosine of a measurement that is not in radians"); 
                    
                    return measurement(std::cos(meas.value_), unitless); 
                
                }


                /**
                 * @brief Take the tangent of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement tan(const measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the tangent of a measurement that is not in radians"); 
                    
                    return measurement(std::tan(meas.value_), unitless);

                }


                /**
                 * @brief Take the arcsine of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement asin(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the arcsine of a measurement that is not unitless"); 
                    
                    return measurement(std::asin(meas.value_), rad);

                }


                /**
                 * @brief Take the arccosine of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement acos(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the arccosine of a measurement that is not unitless"); 
                    
                    return measurement(std::acos(meas.value_), rad);

                }


                /**
                 * @brief Take the arctangent of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement atan(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the arctangent of a measurement that is not unitless"); 
                    
                    return measurement(std::atan(meas.value_), rad);

                }


                /**
                 * @brief Take the hyperbolic sine of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement sinh(const measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the hyperbolic sine of a measurement that is not in radians"); 
                    
                    return measurement(std::sinh(meas.value_), unitless);

                }


                /**
                 * @brief Take the hyperbolic cosine of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement cosh(const measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the hyperbolic cosine of a measurement that is not in radians"); 
                    
                    return measurement(std::cosh(meas.value_), unitless);
                
                }


                /**
                 * @brief Take the hyperbolic tangent of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement tanh(const measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the hyperbolic tangent of a measurement that is not in radians"); 
                    
                    return measurement(std::tanh(meas.value_), unitless);

                }


                /**
                 * @brief Take the hyperbolic arcsine of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement asinh(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the hyperbolic arcsine of a measurement that is not unitless"); 
                    
                    return measurement(std::asinh(meas.value_), rad);

                }


                /**
                 * @brief Take the hyperbolic arccosine of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement acosh(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the hyperbolic arccosine of a measurement that is not unitless"); 
                    
                    return measurement(std::acosh(meas.value_), rad);
                
                }


                /**
                 * @brief Take the hyperbolic arctangent of a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr measurement atanh(const measurement& meas) { 
                    
                    if (meas.units() != unitless) 
                        throw std::runtime_error("Cannot take the hyperbolic arctangent of a measurement that is not unitless"); 
                    
                    return measurement(std::atanh(meas.value_), rad);

                }


            // =============================================                                                                                         
            // get methods
            // =============================================  

                /**
                 * @brief Get the sign of the measurement 
                 * 
                 * @return constexpr int32_t
                 */
                constexpr int32_t sign() const noexcept { 
                    
                    if (value_ == 0.) 
                        return 0; 
                    else 
                        return (value_ > 0) ? 1 : -1; 
                    
                }


                /**
                 * @brief Get the value of the measurement
                 * 
                 * @return constexpr const scalar
                 */
                constexpr scalar value() const noexcept { 
                    
                    return value_; 
                
                }                        
                
                
                /**
                 * @brief Get the value of the measurement
                 * 
                 * @return constexpr scalar& 
                 */
                constexpr scalar& value() noexcept { 
                    
                    return value_; 
                
                }  


                /**
                 * @brief Get the value of the measurement expressed in another units of the measurement
                 * 
                 * @param desired_units 
                 * 
                 * @return constexpr scalar 
                 */
                constexpr scalar value_as(const unit& desired_units) const { 
                    
                    return (units_ == desired_units) ? value_ : units_.convert(value_, desired_units); 
                
                }


                /**
                 * @brief Get the units of the measurement
                 * 
                 * @return constexpr unit 
                 */
                constexpr unit units() const noexcept { 
                    
                    return units_; 
                
                }


                /**
                 * @brief Get the units of the measurement
                 * 
                 * @return constexpr unit&
                 */
                constexpr unit& units() noexcept { 
                    
                    return units_; 
                
                }


                /**
                 * @brief Get the measurement object
                 * 
                 * @return constexpr measurement 
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
                 */
                constexpr measurement convert_to(const unit& desired_units) const { 
                    
                    return measurement(units_.convert(value_, desired_units), desired_units); 
                
                }


                /**
                 * @brief Print the measurement
                 * 
                 * @param newline: if set to true it prints a newline character at the end of the measurement
                 */
                inline void print(const bool& newline = true) const noexcept { 

                    std::cout << value_ << " " << units_;
                    if (newline) std::cout << "\n"; 

                }

            
            private:

            // =============================================                                                                                         
            // class members
            // =============================================  

                scalar value_; ///< The numerical value of the measurement

                unit units_; ///< The units of the measurement


            // =============================================                                                                                         
            // friend class
            // =============================================  

                friend class uncertain_measurement;
                

        }; // class measurement


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
                    units_() {}


                /**
                 * @brief Construct a new uncertain measurement object
                 * 
                 * @param value: the numerical value of the measurement
                 * @param uncertainty: the uncertainty of the measurement
                 * @param units: the units of the measurement
                 */
                explicit constexpr uncertain_measurement(const scalar& val, 
                                                         const scalar& uncertainty_val, 
                                                         const unit& unit) :

                    value_{val}, 
                    uncertainty_{uncertainty_val}, 
                    units_(unit) {
                        
                        if (uncertainty_ < 0.0) 
                            throw std::invalid_argument("Uncertainty cannot be negative");

                    }


                /**
                 * @brief Construct a new uncertain measurement object with 0.0 uncertainty
                 * 
                 * @param value: the numerical value of the measurement
                 * @param units: the units of the measurement
                 */
                explicit constexpr uncertain_measurement(const scalar& val, 
                                                         const unit& unit) noexcept :

                    value_{val}, 
                    uncertainty_{0.0},
                    units_(unit) {}


                /**
                 * @brief Construct a new uncertain measurement object from a measurement and an uncertainty
                 * 
                 * @param measurement: measurement (value and unit)
                 * @param uncertainty: uncertainty
                 */
                explicit constexpr uncertain_measurement(const measurement& other, 
                                                         const scalar& uncertainty_val) : 

                    value_{other.value_}, 
                    uncertainty_{uncertainty_val}, 
                    units_(other.units_) {

                        if (uncertainty_ < 0.0) 
                            throw std::invalid_argument("Uncertainty cannot be negative");
                            
                        if (uncertainty_ > std::fabs(value_)) 
                            throw std::invalid_argument("Uncertainty cannot be greater than the value");

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
                    units_(value.units_) {

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
                 */
                constexpr uncertain_measurement(const uncertain_measurement& other) noexcept :

                    value_{other.value_},
                    uncertainty_{other.uncertainty_},
                    units_(other.units_) {}


                /**
                 * @brief Move construct a new uncertain measurement object
                 * 
                 * @param other: uncertain measurement to move from as a r-value reference
                 */
                constexpr uncertain_measurement(uncertain_measurement&& other) noexcept :

                    value_{std::move(other.value_)},
                    uncertainty_{std::move(other.uncertainty_)},
                    units_(std::move(other.units_)) {}


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
                 */
                constexpr uncertain_measurement operator*(const uncertain_measurement& other) const noexcept {

                    scalar tval1 = uncertainty_ / value_;
                    scalar tval2 = other.uncertainty_ / other.value_;
                    scalar ntol = std::sqrt(tval1 * tval1 + tval2 * tval2);
                    scalar nval = value_ * other.value_;

                    return uncertain_measurement(nval, nval * ntol, units_ * other.units_);

                }
        

                /**
                 * @brief Compute a product and calculate the new uncertainties using the simple method for uncertainty propagation
                 * 
                 * @param other: uncertain_measurement to multiply with
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement simple_product(const uncertain_measurement& other) const noexcept {
                    
                    scalar ntol = uncertainty_ / value_ + other.uncertainty_ / other.value_;
                    scalar nval = value_ * other.value_;

                    return uncertain_measurement(nval, nval * ntol, units_ * other.units_);
                
                }
                

                /**
                 * @brief Compute a product with a measurement, equivalent to uncertain_measurement multiplication with 0 uncertainty
                 * 
                 * @param other: measurement to multiply with
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement operator*(const measurement& other) const noexcept {
                    
                    return uncertain_measurement(value_ * other.value_, other.value_ * uncertainty_, units_ * other.units_);
                
                }


                /**
                 * @brief Compute a product with a scalar and calculate the new uncertainties using the simple uncertainty multiplication method
                 * 
                 * @param val: the scalar to multiply with
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement operator*(const scalar& val) const noexcept { 
                    
                    return uncertain_measurement(value_ * val, uncertainty_ * val, units_); 
                
                }


                /**
                 * @brief Compute a division and calculate the new uncertainties using the root sum of squares(rss) method
                 * 
                 * @param other: uncertain_measurement to divide by
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement operator/(const uncertain_measurement& other) const {
                    
                    if (other.value_ == 0.0) 
                        throw std::invalid_argument("Cannot divide uncertain_measurement by 0");

                    scalar tval1 = uncertainty_ / value_;
                    scalar tval2 = other.uncertainty_ / other.value_;
                    scalar ntol = std::sqrt(tval1 * tval1 + tval2 * tval2);
                    scalar nval = value_ / other.value_;

                    return uncertain_measurement(nval, nval * ntol, units_ / other.units_);
                
                }


                /**
                 * @brief Compute a division and calculate the new uncertainties using simple method for uncertainty propagation
                 * 
                 * @param other: uncertain_measurement to divide by
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement simple_divide(const uncertain_measurement& other) const {
                    
                    if (other.value_ == 0.0) 
                        throw std::invalid_argument("Cannot divide uncertain_measurement by 0");

                    scalar ntol = uncertainty_ / value_ + other.uncertainty_ / other.value_;
                    scalar nval = value_ / other.value_;

                    return uncertain_measurement(nval, nval * ntol, units_ / other.units_);
                
                }


                /**
                 * @brief Compute a division with a measurement, equivalent to uncertain_measurement division with 0 uncertainty
                 * 
                 * @param other: measurement to divide by
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement operator/(const measurement& other) const {
                    
                    if (other.value_ == 0.0) 
                        throw std::invalid_argument("Cannot divide uncertain_measurement by 0");

                    return uncertain_measurement(value_ / other.value_, uncertainty_ / other.value_, units_ / other.units_);
                
                }


                /**
                 * @brief Divide this uncertain_measurement with a scalar
                 * 
                 * @param val: the scalar to divide with
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement operator/(const scalar& val) const {

                    if (val == 0.0) 
                        throw std::invalid_argument("Cannot divide uncertain_measurement by 0");

                    return uncertain_measurement(value_ / val, uncertainty_ / val, units_);
                
                }


                /**
                 * @brief Compute an addition and calculate the new uncertainties using the root sum of squares(rss) method
                 * 
                 * @param other: uncertain_measurement to sum with
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement operator+(const uncertain_measurement& other) const {
                    
                    if (units_.base() != other.units_.base()) 
                        throw std::invalid_argument("Cannot add uncertain_measurements with different unit bases");

                    scalar cval = other.units_.convertion_factor(units_);
                    scalar ntol = std::sqrt(uncertainty_ * uncertainty_ + cval * cval * other.uncertainty_ * other.uncertainty_);

                    return uncertain_measurement(value_ + cval * other.value_, ntol, units_);
               
                }


                /**
                 * @brief Compute an addition and calculate the new uncertainties using the simple uncertainty summation method
                 * 
                 * @param other: uncertain_measurement to sum with
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement simple_add(const uncertain_measurement& other) const {
                    
                    if (units_.base() != other.units_.base()) 
                        throw std::invalid_argument("Cannot add uncertain_measurements with different unit bases");

                    scalar cval = other.units_.convertion_factor(units_);
                    scalar ntol = uncertainty_ + other.uncertainty_ * cval;

                    return uncertain_measurement(value_ + cval * other.value_, ntol, units_);
                
                }


                /**
                 * @brief Compute an addition with a measurement and calculate the new uncertainties using the simple uncertainty summation method
                 * 
                 * @param other: measurement to sum with
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement operator+(const measurement& other) const {
                    
                    if (units_.base() != other.units_.base()) 
                        throw std::invalid_argument("Cannot add uncertain_measurement and measurement with different unit bases");

                    return uncertain_measurement(value_ + other.value_as(units_), uncertainty_, units_);
                
                }


                /**
                 * @brief Return the opposite of this uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement 
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
                 */
                constexpr uncertain_measurement operator-(const uncertain_measurement& other) const {
                    
                    if (units_.base() != other.units_.base()) 
                        throw std::invalid_argument("Cannot subtract uncertain_measurements with different unit bases");

                    scalar cval = other.units_.convertion_factor(units_);
                    scalar ntol = std::sqrt(uncertainty_ * uncertainty_ + cval * cval * other.uncertainty_ * other.uncertainty_);

                    return uncertain_measurement(value_ - cval * other.value_, ntol, units_);
               
                }


                /**
                 * @brief Compute a subtraction and calculate the new uncertainties using the simple uncertainty summation method
                 * 
                 * @param other: uncertain_measurement to subtract with
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement simple_subtract(const uncertain_measurement& other) const {
                    
                    if (units_.base() != other.units_.base()) 
                        throw std::invalid_argument("Cannot subtract uncertain_measurements with different unit bases");

                    auto cval = other.units_.convertion_factor(units_);
                    scalar ntol = uncertainty_ + other.uncertainty_ * cval;

                    return uncertain_measurement(value_ - cval * other.value_, ntol, units_);
                
                }

                
                /**
                 * @brief Compute a subtraction with a measurement and calculate the new uncertainties using the simple uncertainty summation method
                 * 
                 * @param other: measurement to subtract with
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement operator-(const measurement& other) const {

                    if (units_.base() != other.units_.base()) 
                        throw std::invalid_argument("Cannot subtract uncertain_measurement and measurement with different unit bases");

                    return uncertain_measurement(value_ - other.value_as(units_), uncertainty_, units_);
                
                }
  

                /**
                 * @brief Equality operator between uncertain_measurement and measurement
                 * 
                 * @param other: measurement to compare as l-value const reference
                 * 
                 * @return bool
                 */
                constexpr bool operator==(const measurement& other) const noexcept {
                    
                    if (uncertainty_ == 0.0) 
                        return (value_ == other.value_as(units_)) ? true : math::tools::compare_round_equals(value_, other.value_as(units_)); 
                    else 
                        return (other.value_as(units_) >= (value_ - uncertainty_) && other.value_as(units_) <= (value_ + uncertainty_));
               
                }


                /**
                 * @brief Equality operator
                 * 
                 * @param other: uncertain_measurement to compare as l-value const reference
                 * 
                 * @return bool
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
                 */
                constexpr bool operator>(const measurement& other) const noexcept { 
                    
                    return value_ > other.value_as(units_); 
                    
                }


                /** 
                 * @brief More than operator
                 * 
                 * @param val: scalar to compare as l-value const reference
                 * 
                 * @return bool
                 */
                constexpr bool operator>(const scalar& val) const noexcept { 
                    
                    return value_ > val; 
                    
                }

                
                /**
                 * @brief Less than operator
                 * 
                 * @param other: uncertain_measurement to compare as l-value const reference
                 * 
                 * @return bool
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
                 */
                constexpr bool operator<(const measurement& other) const noexcept { 
                    
                    return value_ < other.value_as(units_); 
                    
                }

                
                /**
                 * @brief Less than operator
                 * 
                 * @param val: scalar to compare as l-value const reference
                 * 
                 * @return bool
                 */
                constexpr bool operator<(const scalar& val) const noexcept { 
                    
                    return value_ < val; 
                    
                }

                
                /**
                 * @brief More than or equal operator
                 * 
                 * @param other: uncertain_measurement to compare as l-value const reference
                 * 
                 * @return bool
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
                 */
                constexpr bool operator>=(const measurement& other) const noexcept {
                    
                    return (value_ >= other.value_as(units_)) ? true : operator==(measurement(other.value_as(units_), units_));
                
                }


                /**
                 * @brief More than or equal operator
                 * 
                 * @param val: scalar to compare as l-value const reference
                 * 
                 * @return bool
                 */
                constexpr bool operator>=(const scalar& val) const noexcept { 

                    return value_ >= val; 
                    
                }


                /**
                 * @brief Less than or equal operator
                 * 
                 * @param other: uncertain_measurement to compare as l-value const reference
                 * 
                 * @return bool
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
                 */
                constexpr bool operator<=(const measurement& other) const noexcept {
                    
                    return (value_ <= other.value_as(units_)) ? true : operator==(measurement(other.value_as(units_), units_));
               
                }


                /**
                 * @brief Less than or equal operator
                 * 
                 * @param val: scalar to compare as l-value const reference
                 * 
                 * @return bool
                 */
                constexpr bool operator<=(const scalar& val) const noexcept { 

                    return value_ <= val; 
                    
                }

                
                /**
                 * @brief Perform a product between a measurement and an uncertain_measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * @param umeas: uncertain_measurement
                 *  
                 * @return constexpr uncertain_measurement 
                 */
                friend constexpr uncertain_measurement operator*(const measurement& meas, 
                                                                 const uncertain_measurement& umeas) noexcept { 
                                    
                    return umeas.operator*(meas); 
                                
                }


                /**
                 * @brief Perform a product between a scalar and an uncertain_measurement
                 * 
                 * @param value: scalar
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                friend constexpr uncertain_measurement operator*(const scalar& value, 
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
                 */
                friend constexpr uncertain_measurement operator/(const measurement& meas, 
                                                                 const uncertain_measurement& umeas) {
                                                            
                    if (umeas.value() == 0.0) 
                        throw std::runtime_error("Cannot divide a measurement by a zero uncertain_measurement");

                    scalar ntol = umeas.uncertainty() / umeas.value();
                    scalar nval = meas.value() / umeas.value();

                    return uncertain_measurement(nval, nval * ntol, meas.units() / umeas.units());
                
                }


                /**
                 * @brief Perform a division between a scalar and an uncertain_measurement
                 * 
                 * @param value: scalar
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                friend constexpr uncertain_measurement operator/(const scalar& v1, 
                                                                 const uncertain_measurement& umeas) {

                    if (umeas.value() == 0.0) 
                        throw std::runtime_error("Cannot divide a scalar by a zero uncertain_measurement");

                    scalar ntol = umeas.uncertainty() / umeas.value();
                    scalar nval = v1 / umeas.value();
                    return uncertain_measurement(nval, nval * ntol, umeas.units().inv());
                
                }

                
                /**
                 * @brief Perform a sum between a measurement and an uncertain_measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * @param umeas: uncertain_measurement
                 *  
                 * @return constexpr uncertain_measurement 
                 */
                friend constexpr uncertain_measurement operator+(const measurement& meas, 
                                                                 const uncertain_measurement& umeas) {
                    
                    if (meas.units().base() != umeas.units().base()) 
                        throw std::invalid_argument("Cannot sum measurement and uncertain_measurement with different unit bases");

                    scalar cval = umeas.units().convertion_factor(meas.units());
                    scalar ntol = umeas.uncertainty() * cval;

                    return uncertain_measurement(meas.value() + cval * umeas.value(), ntol, meas.units());
                
                }


                /**
                 * @brief Perform a subtraction between a measurement and an uncertain_measurement
                 * 
                 * @param meas: measurement 
                 * @param umeas: uncertain_measurement
                 *  
                 * @return constexpr uncertain_measurement 
                 */
                friend constexpr uncertain_measurement operator-(const measurement& meas, 
                                                                 const uncertain_measurement& umeas) {

                    if (meas.units().base() != umeas.units().base()) 
                        throw std::invalid_argument("Cannot subtract measurement and uncertain_measurement with different unit bases");

                    scalar cval = umeas.units().convertion_factor(meas.units());
                    scalar ntol = umeas.uncertainty() * cval;

                    return uncertain_measurement(meas.value() - cval * umeas.value(), ntol, meas.units());
                
                }


                /**
                 * @brief Equality operator between a measurement and an uncertain_measurement
                 * 
                 * @param meas: measurement
                 * @param umeas: uncertain_measurement
                 * 
                 * @return bool
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

                    scalar abs_value = std::fabs(umeas.value_);
                    
                    // first significative digit positions
                    int32_t n_val = ((umeas.uncertainty_ >= 1) ? 
                                        std::ceil(std::log10(abs_value)) : 
                                        ((abs_value >= 1) ? 
                                            std::ceil(std::log10(abs_value)) : 
                                            std::floor(std::log10(abs_value)))); 

                    int32_t n_unc = ((umeas.uncertainty_ >= 1) ? 
                                        std::ceil(std::log10(umeas.uncertainty_)) : 
                                        std::floor(std::log10(umeas.uncertainty_))); 

                    int32_t prec = (n_unc > n_val) ? 0 : n_val - n_unc; 

                    bool scientific_notation_needed = (abs_value >= 1e4) || 
                                                      (abs_value <= 1e-4) || 
                                                      (umeas.uncertainty_ >= 1e4) || 
                                                      (umeas.uncertainty_ <= 1e-4);

                    // check if the uncertainty needs to be printed
                    if (umeas.uncertainty_ == 0.0) 
                        os << umeas.as_measurement(); 
                    
                    // check if scientific notation is needed
                    if (scientific_notation_needed) {

                        os << std::scientific; 
                        os << std::setprecision(prec) << umeas.value_ << " ± "; 
                        os << std::setprecision(0) << umeas.uncertainty_ << " " << umeas.units_;

                    } else {

                        os << std::fixed; 

                        if (umeas.uncertainty_ >= 1.) 
                            os << std::setprecision(0); 
                        else 
                            os << std::setprecision(std::fabs(n_unc)); 
                            
                        os << umeas.value_ << " ± " << umeas.uncertainty_ << " " << umeas.units_;

                    }

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

                    scalar abs_value = std::fabs(umeas.value_);
                    
                    // first significative digit positions
                    int32_t n_val = ((umeas.uncertainty_ >= 1) ? 
                                        std::ceil(std::log10(abs_value)) : 
                                        ((abs_value >= 1) ? 
                                            std::ceil(std::log10(abs_value)) : 
                                            std::floor(std::log10(abs_value)))); 

                    int32_t n_unc = ((umeas.uncertainty_ >= 1) ? 
                                        std::ceil(std::log10(umeas.uncertainty_)) : 
                                        std::floor(std::log10(umeas.uncertainty_))); 

                    bool scientific_notation_needed = (abs_value >= 1e4) || 
                                                      (abs_value <= 1e-4) || 
                                                      (umeas.uncertainty_ >= 1e4) || 
                                                      (umeas.uncertainty_ <= 1e-4);

                    // check if the uncertainty needs to be printed
                    if (umeas.uncertainty_ == 0.0) 
                        file << umeas.as_measurement(); 
                    
                    // check if scientific notation is needed
                    if (scientific_notation_needed) {

                        file << std::scientific; 
                        file << std::setprecision((n_unc > n_val) ? 0 : n_val - n_unc) << umeas.value_ << " ± "; 
                        file << std::setprecision(0) << umeas.uncertainty_ << " " << umeas.units_;

                    } else {

                        file << std::fixed; 
                            file << std::fixed;
                        file << std::fixed; 

                        if (umeas.uncertainty_ >= 1.) 
                            file << std::setprecision(0); 
                        else 
                            file << std::setprecision(std::fabs(n_unc)); 
                            
                        file << umeas.value_ << " ± " << umeas.uncertainty_ << " " << umeas.units_;

                    }

                    file << std::defaultfloat << std::setprecision(6);
                    
                    return file; 
                    
                }


            // =============================================
            // operations
            // ============================================= 

                /**
                 * @brief Get the absolute uncertain measurement object
                 * 
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement
                 */
                friend constexpr uncertain_measurement abs(const uncertain_measurement& umeas) noexcept { 
                    
                    return (umeas.value_ < 0.0) ? -umeas : umeas; 
                
                }


                /**                  
                 * @brief Invert the uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 * @note Cannot invert an uncertain_measurement with a zero value
                 * @note The uncertainty is not inverted
                 */
                constexpr uncertain_measurement inv() const { 
                    
                    if (value_ == 0) 
                        throw std::runtime_error("Cannot invert an uncertain_measurement with a zero value");

                    return uncertain_measurement(1 / value_, uncertainty_ / std::pow(value_, 2), units_.inv());
                    
                } 


                /**
                 * @brief Take the power of the uncertain_measurement
                 * 
                 * @param power 
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement pow(const int& power) const noexcept { 
                    
                    return uncertain_measurement(std::pow(value_, power), std::fabs(power * std::pow(value_, power - 1)) * uncertainty_, units_.pow(power)); 
                    
                }

                
                /**
                 * @brief Take the square of the uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                friend constexpr uncertain_measurement square(const uncertain_measurement& umeas) noexcept { 
                    
                    return uncertain_measurement(std::pow(umeas.value_, 2), 2. * std::fabs(umeas.value_) * umeas.uncertainty_, umeas.units_.square()); 
                    
                }

                
                /**
                 * @brief Take the cube of the uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                friend constexpr uncertain_measurement cube(const uncertain_measurement& umeas) noexcept { 
                    
                    return uncertain_measurement(std::pow(umeas.value_, 3), 3. * std::pow(umeas.value_, 2) * umeas.uncertainty_, umeas.units_.cube()); 
                    
                }

                
                /**
                 * @brief Take the root power of the uncertain_measurement
                 * 
                 * @param power 
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement root(const int& power) const { 
                    
                    return uncertain_measurement(std::pow(value_, 1.0 / power), std::fabs(std::pow(value_, 1.0 / power - 1)) * uncertainty_ / power, units_.root(power)); 
                    
                }

                
                /**
                 * @brief Take the square root of the uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement
                 *  
                 */
                friend constexpr uncertain_measurement sqrt(const uncertain_measurement& umeas) { 
                    
                    return uncertain_measurement(std::sqrt(umeas.value_), umeas.uncertainty_ / (2. * std::sqrt(umeas.value_)), umeas.units_.sqrt()); 
                    
                }


                /**
                 * @brief Take the cubic root of the uncertain_measurement
                 * 
                 * @return constexpr uncertain_measurement
                 *  
                 */                
                friend constexpr uncertain_measurement cbrt(const uncertain_measurement& umeas) { 
                    
                    return uncertain_measurement(std::cbrt(umeas.value_), std::pow(umeas.value_, - 2. / 3.) * umeas.uncertainty_ / 3., umeas.units_.cbrt());
                    
                }


                /**
                 * @brief Take the sine of an uncertain_measurement
                 * 
                 * @param umeas: uncertain_measurement 
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement sin(const uncertain_measurement& umeas) { 
                    
                    if (umeas.units() != rad) 
                        throw std::runtime_error("Cannot take the sine of an uncertain_measurement that is not in radians"); 

                    return uncertain_measurement(std::sin(umeas.value_), std::fabs(std::cos(umeas.value_) * umeas.uncertainty_), unitless); 
                
                }


                /**
                 * @brief Take the cosine of an uncertain_measurement
                 * 
                 * @param umeas: uncertain_measurement 
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement cos(const uncertain_measurement& umeas) { 
                    
                    if (umeas.units() != rad) 
                        throw std::runtime_error("Cannot take the cosine of an uncertain_measurement that is not in radians"); 

                    return uncertain_measurement(std::cos(umeas.value_), std::fabs(-std::sin(umeas.value_) * umeas.uncertainty_), unitless); 
                
                }


                /**
                 * @brief Take the tangent of an uncertain_measurement
                 * 
                 * @param umeas: uncertain_measurement 
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement tan(const uncertain_measurement& meas) { 
                    
                    if (meas.units() != rad) 
                        throw std::runtime_error("Cannot take the tangent of an uncertain_measurement that is not in radians");

                    
                    return uncertain_measurement(std::tan(meas.value_), (1 + std::pow(meas.value_, 2)) * meas.uncertainty_, unitless);

                }


                /**
                 * @brief Take the arcsine of an uncertain_measurement
                 * 
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement asin(const uncertain_measurement& umeas) { 
                    
                    if (umeas.units() != unitless) 
                        throw std::runtime_error("Cannot take the arcsine of an uncertain_measurement that is not unitless"); 
                    
                    return uncertain_measurement(std::asin(umeas.value_), umeas.uncertainty_ / std::sqrt(1 - std::pow(umeas.value_, 2)), rad);

                }


                /**
                 * @brief Take the arccosine of an uncertain_measurement
                 * 
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement acos(const uncertain_measurement& umeas) { 
                    
                    if (umeas.units() != unitless) 
                        throw std::runtime_error("Cannot take the arccosine of an uncertain_measurement that is not unitless"); 
                    
                    return uncertain_measurement(std::acos(umeas.value_), umeas.uncertainty_ / std::sqrt(1 - std::pow(umeas.value_, 2)), rad);

                }


                /**
                 * @brief Take the arctangent of an uncertain_measurement
                 * 
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement atan(const uncertain_measurement& umeas) { 
                    
                    if (umeas.units() != unitless) 
                        throw std::runtime_error("Cannot take the arctangent of an uncertain_measurement that is not unitless"); 
                    
                    return uncertain_measurement(std::atan(umeas.value_), umeas.uncertainty_ / (1 + std::pow(umeas.value_, 2)), rad);

                }


                /**
                 * @brief Take the hyperbolic sine of an uncertain_measurement
                 * 
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement sinh(const uncertain_measurement& umeas) { 
                    
                    if (umeas.units() != rad) 
                        throw std::runtime_error("Cannot take the hyperbolic sine of an uncertain_measurement that is not in radians"); 
                    
                    return uncertain_measurement(std::sinh(umeas.value_), std::cosh(umeas.value_) * umeas.uncertainty_, unitless);

                }


                /**
                 * @brief Take the hyperbolic cosine of an uncertain_measurement
                 * 
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement cosh(const uncertain_measurement& umeas) { 
                    
                    if (umeas.units() != rad) 
                        throw std::runtime_error("Cannot take the hyperbolic cosine of an uncertain_measurement that is not in radians"); 
                    
                    return uncertain_measurement(std::cosh(umeas.value_), std::fabs(std::sinh(umeas.value_)) * umeas.uncertainty_, unitless);
                
                }


                /**
                 * @brief Take the hyperbolic tangent of an uncertain_measurement
                 * 
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement tanh(const uncertain_measurement& umeas) { 
                    
                    if (umeas.units() != rad) 
                        throw std::runtime_error("Cannot take the hyperbolic tangent of an uncertain_measurement that is not in radians"); 
                    
                    return uncertain_measurement(std::tanh(umeas.value_), std::fabs((1 - std::pow(umeas.value_, 2))) * umeas.uncertainty_, unitless);

                }


                /**
                 * @brief Take the hyperbolic arcsine of an uncertain_measurement
                 * 
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement asinh(const uncertain_measurement& umeas) { 
                    
                    if (umeas.units() != unitless) 
                        throw std::runtime_error("Cannot take the hyperbolic arcsine of an uncertain_measurement that is not unitless"); 
                    
                    return uncertain_measurement(std::asinh(umeas.value_), umeas.uncertainty_ / std::sqrt(std::pow(umeas.value_, 2) + 1), rad);

                }


                /**
                 * @brief Take the hyperbolic arccosine of an uncertain_measurement
                 * 
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement acosh(const uncertain_measurement& umeas) { 
                    
                    if (umeas.units() != unitless) 
                        throw std::runtime_error("Cannot take the hyperbolic arccosine of an uncertain_measurement that is not unitless"); 
                    
                    return uncertain_measurement(std::acosh(umeas.value_), umeas.uncertainty_ / std::fabs(std::sqrt(std::pow(umeas.value_, 2) - 1)), rad);
                
                }


                /**
                 * @brief Take the hyperbolic arctangent of an uncertain_measurement
                 * 
                 * @param umeas: uncertain_measurement
                 * 
                 * @return constexpr measurement
                 */
                friend constexpr uncertain_measurement atanh(const uncertain_measurement& umeas) { 
                    
                    if (umeas.units() != unitless) 
                        throw std::runtime_error("Cannot take the hyperbolic arctangent of an uncertain_measurement that is not unitless"); 
                    
                    return uncertain_measurement(std::atanh(umeas.value_), umeas.uncertainty_ / std::fabs(std::sqrt(1 - std::pow(umeas.value_, 2))), rad);

                }


            // =============================================                                                                                         
            // set & get methods
            // =============================================  
                
                /**
                 * @brief Get the value of the measurement
                 * 
                 * @return constexpr const scalar
                 */
                constexpr scalar value() const noexcept { 
                    
                    return value_; 
                
                }                        
                
                
                /**
                 * @brief Get the value of the measurement
                 * 
                 * @return constexpr scalar& 
                 */
                constexpr scalar& value() noexcept { 
                    
                    return value_; 
                
                }  


                /**
                 * @brief Get the value of the measurement expressed in another units of the measurement
                 * 
                 * @param desired_units 
                 * @return constexpr scalar 
                 */
                constexpr scalar value_as(const unit& desired_units) const { 
                    
                    return (units_ == desired_units) ? value_ : units_.convert(value_, desired_units); 
                
                }


                /**
                 * @brief Cast to a measurement
                 * 
                 * @return measurement 
                 */
                constexpr measurement as_measurement() const noexcept { 
                    
                    return measurement(value_, units_); 
                    
                }


                /**
                 * @brief Get the uncertainty of the measurement
                 * 
                 * @return constexpr const scalar
                 */
                constexpr scalar uncertainty() const noexcept { 
                    
                    return uncertainty_; 
                
                }                        
                
                
                /**
                 * @brief Get the uncertainty of the measurement
                 * 
                 * @return constexpr scalar& 
                 */
                constexpr scalar& uncertainty() noexcept { 
                    
                    return uncertainty_; 
                
                }  


                /**
                 * @brief Get the uncertainty of the measurement expressed in another units of the measurement
                 * 
                 * @param desired_units 
                 * @return constexpr scalar 
                 */
                constexpr scalar uncertainty_as(const unit& desired_units) const { 
                    
                    return (units_ == desired_units) ? uncertainty_ : units_.convert(uncertainty_, desired_units); 
                
                }


                /**
                 * @brief Get the relative uncertainty of the measurement
                 * 
                 * @return constexpr scalar
                 */
                constexpr scalar relative_uncertainty() const noexcept { 
                    
                    return uncertainty_ / value_;

                }


                /**
                 * @brief Get the weight of the measurement
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement weight() const {
                        
                    return square(this->uncertainty_as_measurement().inv());
    
                }


                /**
                 * @brief Get the uncertainty as a separate measurement
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement uncertainty_as_measurement() const noexcept { 
                    
                    return measurement(uncertainty_, units_); 
                
                }


                /**
                 * @brief Get the units of the measurement
                 * 
                 * @return constexpr unit 
                 */
                constexpr unit units() const noexcept { 
                    
                    return units_; 
                
                }


                /**
                 * @brief Get the units of the measurement
                 * 
                 * @return constexpr unit&
                 */
                constexpr unit& units() noexcept { 
                    
                    return units_; 
                
                }


                /**
                 * @brief Get the uncertain_measurement object
                 * 
                 * @return constexpr uncertain_measurement 
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
                 * @brief Add uncertainty to the uncertain_measurement
                 * 
                 * @param uncertainty: scalar
                 * 
                 * @return void
                 */
                constexpr void add_uncertainty(const scalar new_uncertainty) noexcept { 
                    
                    uncertainty_ = std::sqrt(std::pow(uncertainty_, 2) + std::pow(new_uncertainty, 2));

                }


                /**
                 * @brief Convert the uncertain_measurement to another units
                 * 
                 * @param desired_units: desired unit of measurement 
                 * 
                 * @return constexpr uncertain_measurement 
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
                 */
                void print(const bool& newline = true) const noexcept { 

                    std::cout << *this; 
                    if (newline) 
                        std::cout << "\n";

                }


            private:

            // =============================================                                                                                         
            // class members
            // =============================================  

                scalar value_; ///< the numerical value of the measurement
                
                scalar uncertainty_; ///< the uncertainty of the measurement

                unit units_; ///< the units of the measurement


        }; // class uncertain_measurement


        // =============================================
        // some more operators
        // =============================================

        /**
         * @brief Create a measurement by multiplying a scalar with an unit
         * 
         * @param val: scalar
         * @param units: unit
         * 
         * @return constexpr measurement 
         * 
         */
        constexpr measurement operator*(const scalar& val, 
                                        const unit& units) noexcept { 
                                            
            return measurement(val, units); 
            
        }


        /**
         * @brief Create a measurement by dividing a scalar by an unit
         * 
         * @param val: scalar
         * @param units: unit
         * 
         * @return constexpr measurement 
         * 
         */
        constexpr measurement operator/(const scalar& val, 
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

            scalar value{}; 
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

            scalar value{}, uncertainty{};
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


        using namespace measurements;


        /**
         * @brief Class expressing a generic vector of measurements in a n-dimentional system
         * 
         * @tparam DIM: the number of dimensions
         */
        template <size_t DIM> 
        class vector {


            public: 

            // =============================================
            // constructors and destructor
            // =============================================


                /**
                 * @brief Construct a new vector specifying the unit_base of the vector
                 * 
                 * @param unit: unit_base (default = unit_base())
                 */
                explicit constexpr vector(const unit_base& base = base::default_type) noexcept {

                    data_.fill(measurement(0, unit(base)));

                }


                /**
                 * @brief Construct a new vector from a std::array of measurements
                 * 
                 * @tparam DIM: the number of dimensions
                 * @param data: std::array<measurement, DIM> as l-value const reference
                 */
                constexpr vector(const std::array<measurement, DIM>& data) noexcept : 
                    
                    data_(data) {}


                /**
                 * @brief Construct a new vector from an std::array of measurements
                 * 
                 * @tparam DIM: the number of dimensions
                 * @param data: std::array<measurement, DIM> as r-value reference
                 */
                constexpr vector(std::array<measurement, DIM>&& data) noexcept : 
                    
                    data_(std::move(data)) {}


                /**
                 * @brief Copy construct a new vector from a vector object 
                 * 
                 * @param other: vector as l-value const reference
                 */
                constexpr vector(const vector& other) noexcept : 
                    
                    data_(other.data_) {}


                /**
                 * @brief Move construct a new vector from a vector object
                 * 
                 * @param other: vector as r-value reference
                 */
                constexpr vector(vector&& other) noexcept : 
                    
                    data_(std::move(other.data_)) {}


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
                 * @return constexpr vector&
                 */
                constexpr vector& operator=(const vector& other) noexcept {

                    data_ = other.data_; 

                    return *this; 

                }


                /**
                 * @brief Add a vector to the current vector
                 * 
                 * @note: the two vectors must have the same unit of measurement
                 * 
                 * @param other: vector to add as l-value const reference
                 * 
                 * @return constexpr vector& 
                 */
                constexpr vector& operator+=(const vector& other) {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] += other.data_[i];

                    return *this; 

                }

                
                /**
                 * @brief Subtract a vector to the current vector
                 * 
                 * @param other: vector to subtract as l-value const reference
                 * 
                 * @return constexpr vector& 
                 * 
                 * @note: the two vectors must have the same unit of measurement
                 */
                constexpr vector& operator-=(const vector& other) {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] -= other.data_[i];

                    return *this; 

                }


                /**
                 * @brief Move assignment operator
                 * 
                 * @param other: vector as r-value reference
                 * 
                 * @return constexpr vector& 
                 */
                constexpr vector& operator=(vector&& other) noexcept {
                    
                    data_ = std::move(other.data_); 

                    return *this; 

                }


                /**
                 * @brief Add a vector to the current vector
                 * @note: the two vectors must have the same unit of measurement
                 * 
                 * @param other: vector to add as r-value reference
                 * 
                 * @return constexpr vector& 
                 */
                constexpr vector& operator+=(vector&& other) {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] += std::move(other.data_[i]);

                    return *this; 

                }


                /**
                 * @brief Subtract a vector to the current vector
                 * 
                 * @param other: vector to subtract as r-value reference
                 * 
                 * @return constexpr vector& 
                 * 
                 * @note: the two vectors must have the same unit of measurement
                 */
                constexpr vector& operator-=(vector&& other) {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] -= std::move(other.data_[i]);

                    return *this; 

                }


                /**
                 * @brief Multiply the current vector by a measurement
                 * 
                 * @param meas: measurement to multiply with as l-value const reference
                 * 
                 * @return constexpr vector& 
                 */
                constexpr vector& operator*=(const measurement& meas) noexcept {

                    for (size_t i{}; i < DIM; ++i)  
                        data_[i] *= meas;

                    return *this; 

                }


                /**
                 * @brief Divide the current vector by a measurement
                 * 
                 * @param meas: measurement to divide by as l-value const reference
                 *  
                 * @return constexpr vector& 
                 */
                constexpr vector& operator/=(const measurement& meas) {
                    
                    if (meas.value() == 0.0) 
                        throw std::runtime_error("Cannot divide a vector by a zero measurement");

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] /= meas;
                    
                    return *this; 

                }


                /**
                 * @brief Multiply the current vector by a measurement
                 * 
                 * @param meas: measurement to multiply with as r-value reference
                 * 
                 * @return constexpr vector&
                 */
                constexpr vector& operator*=(measurement&& meas) noexcept {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] *= std::move(meas);

                    return *this; 

                }


                /**
                 * @brief Divide the current vector by a measurement
                 * 
                 * @param meas: measurement to divide by as r-value reference
                 * 
                 * @return constexpr vector&
                 */
                constexpr vector& operator/=(measurement&& meas) {

                    if (meas.value() == 0.0) 
                        throw std::runtime_error("Cannot divide a vector by a zero measurement");
                    
                    for (size_t i{}; i < DIM; ++i)
                        data_[i] /= meas;

                    return *this; 

                }


                /**
                 * @brief Multiply the current vector by a scalar
                 * 
                 * @param scalar: scalar as l-value const reference
                 * 
                 * @return constexpr vector& 
                 */
                constexpr vector& operator*=(const scalar& scalar) noexcept {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] *= scalar;

                    return *this; 

                }


                /**
                 * @brief Divide the current vector by a scalar
                 * 
                 * @param scalar: scalar as l-value const reference
                 *  
                 * @return constexpr vector& 
                 */
                constexpr vector& operator/=(const scalar& scalar) {

                    if (scalar == 0.0) 
                        throw std::runtime_error("Cannot divide a vector by zero");

                    for (size_t i{}; i < DIM; ++i)
                        data_[i] /= scalar;
                    
                    return *this; 
                    
                }


                /**
                 * @brief Return the opposite of the current vector
                 * 
                 * @return constexpr vector 
                 */
                constexpr vector operator-() const noexcept {

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = -data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Sum the current vector and another vector
                 * 
                 * @param other: vector to add as l-value const reference
                 * 
                 * @return constexpr vector 
                 * 
                 * @note: the two vectors must have the same unit of measurement and the same size
                 */
                constexpr vector operator+(const vector& other) const {

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] + other.data_[i]; 
                    
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
                 */
                constexpr vector operator-(const vector& other) const {

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] - other.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Multiply the current vector by a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr vector 
                 */
                constexpr vector operator*(const measurement& meas) const noexcept {

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] * meas; 
                    
                    return result;

                }


                /**
                 * @brief Divide the current vector by a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr vector 
                 */
                constexpr vector operator/(const measurement& meas) const {

                    if (meas.value() == 0.0) 
                        throw std::runtime_error("Cannot divide a vector by a zero measurement");

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] / meas; 
                    
                    return result;

                }


                /**
                 * @brief Multiply the a measurement and a vector
                 * 
                 * @param meas: measurement as l-value const reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 */
                friend constexpr vector operator*(const measurement& meas, const vector& vec) noexcept {

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = meas * vec.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Divide the a measurement and a vector
                 * 
                 * @param meas: measurement as l-value const reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 */
                friend constexpr vector operator/(const measurement& meas, const vector& other) {

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = meas / other.data_[i];
                    
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

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] * scalar; 
                    
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
                constexpr vector operator/(const scalar& scalar) const {

                    if (scalar == 0.0) 
                        throw std::runtime_error("Cannot divide a vector by zero");

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] / scalar; 
                    
                    return result;

                }


                /**
                 * @brief Multiply the a scalar and a vector
                 * 
                 * @param scalar: scalar as l-value const reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 */
                friend constexpr vector operator*(const scalar& scalar, const vector& vec) noexcept {

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = scalar * vec.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Divide the a scalar and a vector
                 * 
                 * @param scalar: scalar as l-value const reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 */                
                friend constexpr vector operator/(const scalar& scalar, const vector& vec) {

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = scalar / vec.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Multiply the current vector with a std::array<scalar>
                 * 
                 * @param scalar_arr: std::array<scalar> as l-value const reference
                 * 
                 * @return constexpr vector 
                 */
                constexpr vector operator*(const std::array<scalar, DIM>& scalar_arr) const noexcept {

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] * scalar_arr[i]; 
                    
                    return result;

                }


                /**
                 * @brief Divide the current vector by a std::array<scalar>
                 * 
                 * @param scalar_arr: std::array<scalar> as l-value const reference
                 * 
                 * @return constexpr vector 
                 */
                constexpr vector operator/(const std::array<scalar, DIM>& scalar_arr) const {

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        if (scalar_arr[i] != 0)
                            result[i] = data_[i] / scalar_arr[i];
                        else throw std::runtime_error("Cannot divide a vector by a vector with a zero component");
                    
                    return result;

                }


                /**
                 * @brief Multiply the a std::array<scalar> and a vector
                 * 
                 * @param scalar_arr: std::array<scalar> as l-value const reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return constexpr vector
                 */
                friend constexpr vector operator*(const std::array<scalar, DIM>& scalar_vec, const vector& other) noexcept {

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = scalar_vec[i] * other.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Divide the a std::array<scalar> and a vector
                 * 
                 * @param scalar_arr: std::array<scalar> as l-value const reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 */
                friend constexpr vector operator/(const std::array<scalar, DIM>& scalar_vec, const vector& other) noexcept {

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = scalar_vec[i] / other.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Equality operator
                 * 
                 * @param other: vector as l-value const reference
                 * 
                 * @return bool
                 */
                constexpr bool operator==(const vector& other) const noexcept {

                    for (size_t i{}; i < DIM; ++i)
                        if (data_[i] != other.data_[i]) 
                            return false;
                    
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

                    for (size_t i{}; i < DIM; ++i) 
                        if (data_[i] != other.data_[i]) 
                            return true;

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
                    
                    if (index >= DIM) 
                        throw std::out_of_range("Cannot access a vector with an index out of range");

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
                 */
                constexpr measurement& operator[](const size_t& index) { 
                    
                    if (index >= DIM) 
                        throw std::out_of_range("Cannot access a vector with an index out of range");

                    return data_[index]; 
                    
                }


                /**
                 * @brief Print the vector to the standard output
                 * 
                 * @param os: ostream as l-value reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return std::ostream&
                 */
                friend std::ostream& operator<<(std::ostream& os, const vector& vec) noexcept {

                    os << "{ ";
                    for (size_t i{}; i < DIM; ++i) 
                        os << vec.data_[i] << (i != DIM - 1 ? ", " : " }");

                    return os;

                }


                /**
                 * @brief Print the vector to file
                 * 
                 * @param file: ofstream as l-value reference
                 * @param vec: vector as l-value const reference
                 * 
                 * @return std::ofstream&
                 */
                friend std::ofstream& operator<<(std::ofstream& file, const vector& vec) noexcept {

                    file << "{ ";
                    for (size_t i{}; i < DIM; ++i) 
                        file << vec.data_[i] << (i != DIM - 1 ? ", " : " }");

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
                 */
                constexpr friend vector cross(const vector& v1, const vector& v2) {

                    std::array<measurement, DIM> cross_vec;
                    for (size_t i{}; i < DIM; ++i) 
                        cross_vec[i] = v1[(i + 1) % v1.size()] * v2[(i + 2) % v1.size()] - v1[(i + 2) % v1.size()] * v2[(i + 1) % v1.size()]; 

                    return cross_vec;
                
                }

                
                /**
                 * @brief Compute the dot product between two vectors
                 * 
                 * @param v1: vector as l-value const reference
                 * @param v2: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 */
                constexpr friend measurement dot(const vector& v1, const vector& v2) noexcept {

                    measurement result(v1[0].units() * v2[0].units());
                    for (size_t i{}; i < v1.size(); ++i)
                        result += v1[i] * v2[i]; 
                    
                    return result;
                
                }


                /**
                 * @brief Invert the vector
                 * 
                 * @return constexpr vector 
                 */
                constexpr vector inv() const {

                    std::array<measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i].inv(); 
                    
                    return result;

                }


                /**
                 * @brief Take the power of the vector
                 * 
                 * @param vec: vector as l-value const reference
                 * @param power: int
                 *  
                 * @return constexpr vector 
                 */
                constexpr vector pow(const vector& vec, const int& power) noexcept {

                    std::array<measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = pow(vec.data_[i], power);
                    
                    return result;

                }


                /**
                 * @brief Take the square of the vector
                 * 
                 * @param vec: vector as l-value const reference

                 * @return constexpr vector 
                 */
                friend constexpr vector square(const vector& vec) noexcept {

                    std::array<measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = square(vec.data_[i]);
                    
                    return result;

                }


                /**
                 * @brief Take the cube of the vector
                 * 
                 * @param vec: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 */
                friend constexpr vector cube(const vector& vec) noexcept {

                    std::array<measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = cube(vec.data_[i]);
                    
                    return result;

                }


                /**
                 * @brief Take the root of the vector
                 * 
                 * @param vec: vector as l-value const reference
                 * @param power: int
                 * 
                 * @return constexpr vector 
                 */
                friend constexpr vector root(const vector& vec, const int& power) {

                    std::array<measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = root(vec.data_[i], power);
                    
                    return result;

                }


                /**
                 * @brief Take the square root of the vector
                 * 
                 * @param vec: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 */
                friend constexpr vector sqrt(const vector& vec) {

                    std::array<measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = sqrt(vec.data_[i]);
                    
                    return result;

                }


                /**
                 * @brief Take the cube root of the vector
                 * 
                 * @param vec: vector as l-value const reference
                 * 
                 * @return constexpr vector 
                 */
                friend constexpr vector cbrt(const vector& vec) {

                    std::array<measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = cbrt(vec.data_[i]);
                    
                    return result;

                }


            // =============================================
            // set & get methods
            // =============================================

                /**
                 * @brief Get the size of the vector
                 * 
                 * @return constexpr size_t 
                 */
                constexpr size_t size() const noexcept { 
                    
                    return DIM; 
                    
                }


                /**
                 * @brief Get the first element of the vector
                 * 
                 * @note the vector must have at least one element
                 * 
                 * @return constexpr measurement
                 */
                constexpr measurement x() const noexcept requires (DIM >= 1) { 
                    
                    return data_[0]; 
                
                }
                                

                /**
                 * @brief Get the second element of the vector
                 * 
                 * @note the vector must have at least two elements
                 * 
                 * @return constexpr measurement
                 */
                constexpr measurement y() const noexcept requires (DIM >= 2) { 
                    
                    return data_[1]; 
                
                }


                /**
                 * @brief Get the third element of the vector
                 * 
                 * @note the vector must have at least three elements
                 *
                 * @return constexpr measurement
                 */
                constexpr measurement z() const noexcept requires (DIM >= 3) { 
                    
                    return data_[2]; 
                
                }


                /**
                 * @brief Get the first element of the vector
                 * 
                 * @note the vector must have at least one element
                 * 
                 * @return constexpr measurement&
                 */
                constexpr measurement& x() noexcept requires (DIM >= 1) { 
                    
                    return data_[0]; 
                
                }
                                

                /**
                 * @brief Get the second element of the vector
                 * 
                 * @note the vector must have at least two elements
                 * 
                 * @return constexpr measurement&
                 */
                constexpr measurement& y() noexcept requires (DIM >= 2) { 
                    
                    return data_[1]; 
                
                }


                /**
                 * @brief Get the third element of the vector
                 * 
                 * @note the vector must have at least three elements
                 * 
                 * @return constexpr measurement&
                 */
                constexpr measurement& z() noexcept requires (DIM >= 3) { 
                    
                    return data_[2]; 
                
                }


                /**
                 * @brief Get the unit of the vector
                 * 
                 * @return constexpr unit 
                 */
                constexpr unit units() const noexcept { 

                    return data_.front().units(); 
                
                }


                /**
                 * @brief Get the data of the vector
                 * 
                 * @return constexpr std::array<measurement> 
                 */
                constexpr std::array<measurement, DIM> data() const noexcept { 
                    
                    return data_; 
                
                }


                /**
                 * @brief Get the data of the vector
                 * 
                 * @return constexpr std::array<measurement, DIM>& 
                 */
                constexpr std::array<measurement, DIM>& data() noexcept { 
                    
                    return data_;
                
                }


                /**
                 * @brief Get the vector object
                 * 
                 * @return constexpr vector 
                 */
                constexpr vector as_vector() const noexcept { 
                    
                    return *this; 
                
                }


                /**
                 * @brief Get the norm of the vector
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement norm() const noexcept { 

                    if constexpr (DIM == 1) 
                        return data_[0];

                    vector squared = square(*this);

                    return sqrt(std::accumulate(squared.data().begin(), squared.data().end(), measurement{0.0, squared.units()}));

                }


                /**
                 * @brief Get the squared norm of the vector
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement norm2() const noexcept { 

                    if constexpr (DIM == 1) 
                        return data_[0]; 

                    vector squared(square(*this));

                    return std::accumulate(squared.data().begin(), squared.data().end(), measurement{0.0, squared.units()});
                    
                }


                /**
                 * @brief Get the normalization of the vector
                 * 
                 * @return constexpr vector 
                 */
                constexpr vector versor() const {

                    return *this / this->norm(); 

                } 


                /**
                 * @brief Get the polar angle
                 * 
                 * @note the vector must have at least two elements
                 * 
                 * @return constexpr measurement
                 */
                constexpr measurement phi() const noexcept requires (DIM >= 2) { 
                    
                    return atan(data_[1] / data_[0]); 
                
                }     
                

                /**
                 * @brief Get the azimuthal angle
                 * 
                 * @note the vector must have at least three elements
                 * 
                 * @return constexpr measurement
                 */
                constexpr measurement theta() const requires (DIM >= 3) { 

                    if (data_[2] == 0.0 * m) 
                        return 0.0 * rad;
                    else    
                        return acos((data_[2] / norm()));

                }
                

                /**
                 * @brief Print the vector to the standard output
                 * 
                 * @param newline: if true, print a newline character at the end
                 * 
                 * @return void
                 */
                constexpr void print(const bool& newline = true) const noexcept {

                    std::cout << "{\n";
                    for (size_t i{}; i < DIM; ++i) 
                        std::cout << "\t" << data_[i] << "\n";
                    std::cout << "}";  
                    if (newline) 
                        std::cout << "\n";

                }   


                /**
                 * @brief Save the vector to a file
                 * 
                 * @param file_name: the name of the file
                 * @param units: desired units for the output
                 * 
                 * @return void
                 */
                void save(const std::string& file_name, const unit& units) const {

                    std::ofstream file_out(file_name, std::ios::app);
                    if (file_out.is_open()) 

                        for (size_t i{}; i < DIM; ++i) {

                            file_out << data_[i].value_as(units); 
                            if (i < DIM - 1) 
                                file_out << ' '; 

                        }

                    else 
                        throw std::invalid_argument("Unable to open '" + file_name + "'");

                    file_out << '\n';
                    file_out.close();

                }


            protected:

            // =============================================
            // class members
            // =============================================

                std::array<measurement, DIM> data_; ///< array of measurement


        }; // class vector


        /**
         * @brief Class expressing a generic vector of uncertain_measurements in a n-dimentional system
         * 
         * @tparam DIM: the number of dimensions
         */
        template <size_t DIM> 
        class uvector {


            public: 

            // =============================================
            // constructors and destructor
            // =============================================


                /// @brief Construct a new uvector object
                explicit constexpr uvector() noexcept : 

                    data_() {}


                /**
                 * @brief Construct a new uvector specifying the unit_base of the uvector
                 * 
                 * @param unit: unit_base (default = unit_base())
                 */
                explicit constexpr uvector(const unit_base& base = base::default_type) noexcept {

                    data_.fill(uncertain_measurement(0., 0., unit(base)));

                }


                /**
                 * @brief Construct a new uvector from a std::array of uncertain_measurements
                 * 
                 * @tparam DIM: the number of dimensions
                 * @param data: std::array<uncertain_measurement, DIM> as l-value const reference
                 */
                explicit constexpr uvector(const std::array<uncertain_measurement, DIM>& data) noexcept : 
                    
                    data_(data) {}


                /**
                 * @brief Construct a new uvector from an std::array of uncertain_measurements
                 * 
                 * @tparam DIM: the number of dimensions
                 * @param data: std::array<uncertain_measurement, DIM> as r-value reference
                 */
                constexpr uvector(std::array<uncertain_measurement, DIM>&& data) noexcept : 
                    
                    data_(std::move(data)) {}


                /**
                 * @brief Copy construct a new uvector from a uvector object 
                 * 
                 * @param other: uvector as l-value const reference
                 */
                constexpr uvector(const uvector& other) noexcept : 
                    
                    data_(other.data_) {}


                /**
                 * @brief Move construct a new uvector from a uvector object
                 * 
                 * @param other: uvector as r-value reference
                 */
                constexpr uvector(uvector&& other) noexcept : 
                    
                    data_(std::move(other.data_)) {}


                /// @brief Default destructor
                ~uvector() = default;


            // =============================================
            // operators
            // =============================================

                /**
                 * @brief Copy assignment operator
                 * 
                 * @param other: uvector as l-value const reference
                 * 
                 * @return constexpr uvector&
                 */
                constexpr uvector& operator=(const uvector& other) noexcept {

                    data_ = other.data_; 

                    return *this; 

                }


                /**
                 * @brief Add a uvector to the current uvector
                 * 
                 * @note: The uvectors must have the same unit of measurement
                 * 
                 * @param other: uvector to add as l-value const reference
                 * 
                 * @return constexpr uvector& 
                 */
                constexpr uvector& operator+=(const uvector& other) {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] += other.data_[i];

                    return *this; 

                }

                
                /**
                 * @brief Subtract a uvector to the current uvector
                 * 
                 * @note: The uvectors must have the same unit of measurement
                 * 
                 * @param other: uvector to subtract as l-value const reference
                 * 
                 * @return constexpr uvector& 
                 */
                constexpr uvector& operator-=(const uvector& other) {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] -= other.data_[i];

                    return *this; 

                }


                /**
                 * @brief Copy assignment operator from a vector object
                 * 
                 * @param other: vector as l-value const reference
                 * 
                 * @return constexpr uvector&
                 */
                constexpr uvector& operator=(const vector<DIM>& other) noexcept {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] = uncertain_measurement(other.data_[i].value(), 0.0, other.data_[i].unit());

                    return *this; 

                }


                /**
                 * @brief Add a vector to the current uvector
                 * 
                 * @note: The uvector and the vector must have the same unit of measurement
                 * 
                 * @param other: vector to add as l-value const reference
                 * 
                 * @return constexpr uvector& 
                 */
                constexpr uvector& operator+=(const vector<DIM>& other) {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] += other.data_[i];

                    return *this; 

                }

                
                /**
                 * @brief Subtract a vector to the current uvector
                 * 
                 * @note: The uvector and the vector must have the same unit of measurement
                 * 
                 * @param other: vector to subtract as l-value const reference
                 * 
                 * @return constexpr uvector& 
                 */
                constexpr uvector& operator-=(const vector<DIM>& other) {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] -= other.data_[i];

                    return *this; 

                }


                /**
                 * @brief Multiply the current uvector by a measurement
                 * 
                 * @param meas: measurement to multiply with as l-value const reference
                 * 
                 * @return constexpr uvector& 
                 */
                constexpr uvector& operator*=(const measurement& meas) noexcept {

                    for (size_t i{}; i < DIM; ++i)  
                        data_[i] *= meas;

                    return *this; 

                }


                /**
                 * @brief Divide the current uvector by a measurement
                 * 
                 * @param meas: measurement to divide by as l-value const reference
                 *  
                 * @return constexpr uvector& 
                 */
                constexpr uvector& operator/=(const measurement& meas) {
                    
                    if (meas.value() == 0.0) 
                        throw std::runtime_error("Cannot divide a uvector by a zero measurement");

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] /= meas;
                    
                    return *this; 

                }


                /**
                 * @brief Multiply the current uvector by a uncertain_measurement
                 * 
                 * @param meas: uncertain_measurement to multiply with as l-value const reference
                 * 
                 * @return constexpr uvector& 
                 */
                constexpr uvector& operator*=(const uncertain_measurement& meas) noexcept {

                    for (size_t i{}; i < DIM; ++i)  
                        data_[i] *= meas;

                    return *this; 

                }


                /**
                 * @brief Divide the current uvector by a uncertain_measurement
                 * 
                 * @param meas: uncertain_measurement to divide by as l-value const reference
                 *  
                 * @return constexpr uvector& 
                 */
                constexpr uvector& operator/=(const uncertain_measurement& meas) {
                    
                    if (meas.value() == 0.0) 
                        throw std::runtime_error("Cannot divide a uvector by a zero uncertain_measurement");

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] /= meas;
                    
                    return *this; 

                }


                /**
                 * @brief Multiply the current uvector by a scalar
                 * 
                 * @param scalar: scalar as l-value const reference
                 * 
                 * @return constexpr uvector& 
                 */
                constexpr uvector& operator*=(const scalar& scalar) noexcept {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] *= scalar;

                    return *this; 

                }


                /**
                 * @brief Divide the current uvector by a scalar
                 * 
                 * @param scalar: scalar as l-value const reference
                 *  
                 * @return constexpr uvector& 
                 */
                constexpr uvector& operator/=(const scalar& scalar) {

                    if (scalar == 0.0) 
                        throw std::runtime_error("Cannot divide a uvector by zero");

                    for (size_t i{}; i < DIM; ++i)
                        data_[i] /= scalar;
                    
                    return *this; 
                    
                }


                /**
                 * @brief Move assignment operator
                 * 
                 * @param other: uvector as r-value reference
                 * 
                 * @return constexpr uvector& 
                 */
                constexpr uvector& operator=(uvector&& other) noexcept {
                    
                    data_ = std::move(other.data_); 

                    return *this; 

                }


                /**
                 * @brief Add a uvector to the current uvector
                 * @note: the two uvectors must have the same unit of measurement
                 * 
                 * @param other: uvector to add as r-value reference
                 * 
                 * @return constexpr uvector& 
                 */
                constexpr uvector& operator+=(uvector&& other) {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] += std::move(other.data_[i]);

                    return *this; 

                }


                /**
                 * @brief Subtract a uvector to the current uvector
                 * 
                 * @param other: uvector to subtract as r-value reference
                 * 
                 * @return constexpr uvector& 
                 * 
                 * @note: the two uvectors must have the same unit of measurement
                 */
                constexpr uvector& operator-=(uvector&& other) {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] -= std::move(other.data_[i]);

                    return *this; 

                }


                /**
                 * @brief Multiply the current uvector by a measurement
                 * 
                 * @param meas: measurement to multiply with as r-value reference
                 * 
                 * @return constexpr uvector&
                 */
                constexpr uvector& operator*=(measurement&& meas) noexcept {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] *= std::move(meas);

                    return *this; 

                }


                /**
                 * @brief Divide the current uvector by a measurement
                 * 
                 * @param meas: measurement to divide by as r-value reference
                 * 
                 * @return constexpr uvector&
                 */
                constexpr uvector& operator/=(measurement&& meas) {

                    if (meas.value() == 0.0) 
                        throw std::runtime_error("Cannot divide a uvector by a zero measurement");
                    
                    for (size_t i{}; i < DIM; ++i)
                        data_[i] /= meas;

                    return *this; 

                }


                /**
                 * @brief Multiply the current uvector by a uncertain_measurement
                 * 
                 * @param meas: uncertain_measurement to multiply with as r-value reference
                 * 
                 * @return constexpr uvector&
                 */
                constexpr uvector& operator*=(uncertain_measurement&& meas) noexcept {

                    for (size_t i{}; i < DIM; ++i) 
                        data_[i] *= std::move(meas);

                    return *this; 

                }


                /**
                 * @brief Divide the current uvector by a uncertain_measurement
                 * 
                 * @param meas: uncertain_measurement to divide by as r-value reference
                 * 
                 * @return constexpr uvector&
                 */
                constexpr uvector& operator/=(uncertain_measurement&& meas) {

                    if (meas.value() == 0.0) 
                        throw std::runtime_error("Cannot divide a uvector by a zero uncertain_measurement");
                    
                    for (size_t i{}; i < DIM; ++i)
                        data_[i] /= meas;

                    return *this; 

                }


                /**
                 * @brief Return the opposite of the current uvector
                 * 
                 * @return constexpr uvector 
                 */
                constexpr uvector operator-() const noexcept {

                    std::array<measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = -data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Sum the current uvector and another uvector
                 * 
                 * @param other: uvector to add as l-value const reference
                 * 
                 * @return constexpr uvector 
                 * 
                 * @note: the two uvectors must have the same unit of measurement and the same size
                 */
                constexpr uvector operator+(const uvector& other) const {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] + other.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Subtract the current uvector and another uvector
                 * 
                 * @param other: uvector to subtract as l-value const reference
                 * 
                 * @return constexpr uvector 
                 * 
                 * @note: the two uvectors must have the same unit of measurement and the same size
                 */
                constexpr uvector operator-(const uvector& other) const {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] - other.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Multiply the current uvector by a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                constexpr uvector operator*(const measurement& meas) const noexcept {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] * meas; 
                    
                    return result;

                }


                /**
                 * @brief Divide the current uvector by a measurement
                 * 
                 * @param meas: measurement as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                constexpr uvector operator/(const measurement& meas) const {

                    if (meas.value() == 0.0) 
                        throw std::runtime_error("Cannot divide a uvector by a zero measurement");

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] / meas; 
                    
                    return result;

                }


                /**
                 * @brief Multiply the current uvector by a scalar
                 * 
                 * @param scalar: scalar as l-value const reference
                 * 
                 * @return constexpr uvector
                 *  
                 */
                constexpr uvector operator*(const scalar& scalar) const noexcept {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] * scalar; 
                    
                    return result;

                }


                /**
                 * @brief Divide the current uvector by a scalar
                 * 
                 * @param scalar: scalar as l-value const reference
                 * 
                 * @return constexpr uvector
                 *  
                 */
                constexpr uvector operator/(const scalar& scalar) const {

                    if (scalar == 0.0) 
                        throw std::runtime_error("Cannot divide a uvector by zero");

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] / scalar; 
                    
                    return result;

                }


                /**
                 * @brief Multiply the current uvector with a std::array<scalar>
                 * 
                 * @param scalar_arr: std::array<scalar> as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                constexpr uvector operator*(const std::array<scalar, DIM>& scalar_arr) const noexcept {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i] * scalar_arr[i]; 
                    
                    return result;

                }


                /**
                 * @brief Divide the current uvector by a std::array<scalar>
                 * 
                 * @param scalar_arr: std::array<scalar> as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                constexpr uvector operator/(const std::array<scalar, DIM>& scalar_arr) const {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        if (scalar_arr[i] != 0)
                            result[i] = data_[i] / scalar_arr[i];
                        else throw std::runtime_error("Cannot divide a uvector by a uvector with a zero component");
                    
                    return result;

                }


                /**
                 * @brief Multiply the a measurement and a uvector
                 * 
                 * @param meas: measurement as l-value const reference
                 * @param vec: uvector as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                friend constexpr uvector operator*(const measurement& meas, const uvector& vec) noexcept {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = meas * vec.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Divide the a measurement and a uvector
                 * 
                 * @param meas: measurement as l-value const reference
                 * @param vec: uvector as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                friend constexpr uvector operator/(const measurement& meas, const uvector& other) {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = meas / other.data_[i];
                    
                    return result;

                }


                /**
                 * @brief Multiply the a uncertain_measurement and a uvector
                 * 
                 * @param meas: uncertain_measurement as l-value const reference
                 * @param vec: uvector as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                friend constexpr uvector operator*(const uncertain_measurement& meas, const uvector& vec) noexcept {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = meas * vec.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Divide the a uncertain_measurement and a uvector
                 * 
                 * @param meas: uncertain_measurement as l-value const reference
                 * @param vec: uvector as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                friend constexpr uvector operator/(const uncertain_measurement& meas, const uvector& other) {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = meas / other.data_[i];
                    
                    return result;

                }


                /**
                 * @brief Multiply the a scalar and a uvector
                 * 
                 * @param scalar: scalar as l-value const reference
                 * @param vec: uvector as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                friend constexpr uvector operator*(const scalar& scalar, const uvector& vec) noexcept {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = scalar * vec.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Divide the a scalar and a uvector
                 * 
                 * @param scalar: scalar as l-value const reference
                 * @param vec: uvector as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */                
                friend constexpr uvector operator/(const scalar& scalar, const uvector& vec) {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = scalar / vec.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Multiply the a std::array<scalar> and a uvector
                 * 
                 * @param scalar_arr: std::array<scalar> as l-value const reference
                 * @param vec: uvector as l-value const reference
                 * 
                 * @return constexpr uvector
                 */
                friend constexpr uvector operator*(const std::array<scalar, DIM>& scalar_vec, const uvector& other) noexcept {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = scalar_vec[i] * other.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Divide the a std::array<scalar> and a uvector
                 * 
                 * @param scalar_arr: std::array<scalar> as l-value const reference
                 * @param vec: uvector as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                friend constexpr uvector operator/(const std::array<scalar, DIM>& scalar_vec, const uvector& other) noexcept {

                    std::array<uncertain_measurement, DIM> result; 
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = scalar_vec[i] / other.data_[i]; 
                    
                    return result;

                }


                /**
                 * @brief Equality operator
                 * 
                 * @param other: uvector as l-value const reference
                 * 
                 * @return bool
                 */
                constexpr bool operator==(const uvector& other) const noexcept {

                    for (size_t i{}; i < DIM; ++i)
                        if (data_[i] != other.data_[i]) 
                            return false;
                    
                    return true;

                }


                /**
                 * @brief Inequality operator
                 * 
                 * @param other: uvector as l-value const reference
                 * 
                 * @return bools
                 *  
                 */
                constexpr bool operator!=(const uvector& other) const noexcept {

                    for (size_t i{}; i < DIM; ++i) 
                        if (data_[i] != other.data_[i]) 
                            return true;

                    return false;

                }
            

                /**
                 * @brief Access the i-th element of the uvector
                 * 
                 * @param index: size_t
                 * 
                 * @return constexpr uncertain_measurement 
                 * 
                 * @note: index must be in the range [0, DIM)
                 */
                constexpr const uncertain_measurement& operator[](const size_t& index) const { 
                    
                    if (index >= DIM) 
                        throw std::out_of_range("Cannot access a uvector with an index out of range");

                    return data_[index]; 
                    
                }


                /**
                 * @brief Access the i-th element of the uvector
                 * 
                 * @param index: size_t 
                 * 
                 * @return constexpr uncertain_measurement& 
                 * 
                 * @note: index must be in the range [0, DIM)
                 */
                constexpr uncertain_measurement& operator[](const size_t& index) { 
                    
                    if (index >= DIM) 
                        throw std::out_of_range("Cannot access a uvector with an index out of range");

                    return data_[index]; 
                    
                }


                /**
                 * @brief Print the uvector to the standard output
                 * 
                 * @param os: ostream as l-value reference
                 * @param vec: uvector as l-value const reference
                 * 
                 * @return std::ostream&
                 */
                friend std::ostream& operator<<(std::ostream& os, const uvector& vec) noexcept {

                    os << "{ ";
                    for (size_t i{}; i < DIM; ++i) 
                        os << vec.data_[i] << (i != DIM - 1 ? ", " : " }");

                    return os;

                }


                /**
                 * @brief Print the uvector to file
                 * 
                 * @param file: ofstream as l-value reference
                 * @param vec: uvector as l-value const reference
                 * 
                 * @return std::ofstream&
                 */
                friend std::ofstream& operator<<(std::ofstream& file, const uvector& vec) noexcept {

                    file << "{ ";
                    for (size_t i{}; i < DIM; ++i) 
                        file << vec.data_[i] << (i != DIM - 1 ? ", " : " }");

                    return file;

                }


            // =============================================
            // operations
            // =============================================

                /**
                 * @brief Compute the cross product between two uvectors
                 * 
                 * @param v1: uvector as l-value const reference
                 * @param v2: uvector as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                constexpr friend uvector cross(const uvector& v1, const uvector& v2) {

                    std::array<uncertain_measurement, DIM> cross_vec;
                    for (size_t i{}; i < DIM; ++i) 
                        cross_vec[i] = v1[(i + 1) % v1.size()] * v2[(i + 2) % v1.size()] - v1[(i + 2) % v1.size()] * v2[(i + 1) % v1.size()]; 

                    return cross_vec;
                
                }

                
                /**
                 * @brief Compute the dot product between two uvectors
                 * 
                 * @param v1: uvector as l-value const reference
                 * @param v2: uvector as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                constexpr friend uncertain_measurement dot(const uvector& v1, const uvector& v2) noexcept {

                    uncertain_measurement result(v1[0].units() * v2[0].units());
                    for (size_t i{}; i < v1.size(); ++i)
                        result += v1[i] * v2[i]; 
                    
                    return result;
                
                }


                /**
                 * @brief Invert the uvector
                 * 
                 * @return constexpr uvector 
                 */
                constexpr uvector inv() const {

                    std::array<uncertain_measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = data_[i].inv(); 
                    
                    return result;

                }


                /**
                 * @brief Take the power of the uvector
                 * 
                 * @param vec: uvector as l-value const reference
                 * @param power: int
                 *  
                 * @return constexpr uvector 
                 */
                constexpr uvector pow(const uvector& vec, const int& power) noexcept {

                    std::array<uncertain_measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = pow(vec.data_[i], power);
                    
                    return result;

                }


                /**
                 * @brief Take the square of the uvector
                 * 
                 * @param vec: uvector as l-value const reference

                 * @return constexpr uvector 
                 */
                friend constexpr uvector square(const uvector& vec) noexcept {

                    std::array<uncertain_measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = square(vec.data_[i]);
                    
                    return result;

                }


                /**
                 * @brief Take the cube of the uvector
                 * 
                 * @param vec: uvector as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                friend constexpr uvector cube(const uvector& vec) noexcept {

                    std::array<uncertain_measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = cube(vec.data_[i]);
                    
                    return result;

                }


                /**
                 * @brief Take the root of the uvector
                 * 
                 * @param vec: uvector as l-value const reference
                 * @param power: int
                 * 
                 * @return constexpr uvector 
                 */

                friend constexpr uvector root(const uvector& vec, const int& power) {

                    std::array<uncertain_measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = root(vec.data_[i], power);
                    
                    return result;

                }


                /**
                 * @brief Take the square root of the uvector
                 * 
                 * @param vec: uvector as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                friend constexpr uvector sqrt(const uvector& vec) {

                    std::array<uncertain_measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = sqrt(vec.data_[i]);
                    
                    return result;

                }


                /**
                 * @brief Take the cube root of the uvector
                 * 
                 * @param vec: uvector as l-value const reference
                 * 
                 * @return constexpr uvector 
                 */
                friend constexpr uvector cbrt(const uvector& vec) {

                    std::array<uncertain_measurement, DIM> result;
                    for (size_t i{}; i < DIM; ++i)
                        result[i] = cbrt(vec.data_[i]);
                    
                    return result;

                }


            // =============================================
            // set & get methods
            // =============================================

                /**
                 * @brief Get the size of the uvector
                 * 
                 * @return constexpr size_t 
                 */
                constexpr size_t size() const noexcept { 
                    
                    return DIM; 
                    
                }


                /**
                 * @brief Get the first element of the uvector
                 * 
                 * @return constexpr uncertain_measurement
                 * 
                 * @note the uvector must have at least one element
                 */
                constexpr uncertain_measurement x() const noexcept requires (DIM >= 1) { 
                    
                    return data_[0]; 
                
                }
                                

                /**
                 * @brief Get the second element of the uvector
                 * 
                 * @return constexpr uncertain_measurement
                 * 
                 * @note the uvector must have at least two elements
                 */
                constexpr uncertain_measurement y() const noexcept requires (DIM >= 2) { 
                    
                    return data_[1]; 
                
                }


                /**
                 * @brief Get the third element of the uvector
                 * 
                 * @return constexpr uncertain_measurement
                 * 
                 * @note the uvector must have at least three elements
                 */
                constexpr uncertain_measurement z() const noexcept requires (DIM >= 3) { 
                    
                    return data_[2]; 
                
                }


                /**
                 * @brief Get the first element of the uvector
                 * 
                 * @return constexpr uncertain_measurement&
                 * 
                 * @note the uvector must have at least one element
                 */
                constexpr uncertain_measurement& x() noexcept requires (DIM >= 1) { 
                    
                    return data_[0]; 
                
                }
                                

                /**
                 * @brief Get the second element of the uvector
                 * 
                 * @return constexpr uncertain_measurement&
                 * 
                 * @note the uvector must have at least two elements
                 */
                constexpr uncertain_measurement& y() noexcept requires (DIM >= 2) { 
                    
                    return data_[1]; 
                
                }


                /**
                 * @brief Get the third element of the uvector
                 * 
                 * @return constexpr uncertain_measurement&
                 * 
                 * @note the uvector must have at least three elements
                 */
                constexpr uncertain_measurement& z() noexcept requires (DIM >= 3) { 
                    
                    return data_[2]; 
                
                }


                constexpr unit units() const noexcept { 

                    return data_.front().units(); 
                
                }


                /**
                 * @brief Get the data of the uvector
                 * 
                 * @return constexpr std::array<uncertain_measurement> 
                 */
                constexpr std::array<uncertain_measurement, DIM> data() const noexcept { 
                    
                    return data_; 
                
                }


                /**
                 * @brief Get the data of the uvector
                 * 
                 * @return constexpr std::array<uncertain_measurement>& 
                 */
                constexpr std::array<uncertain_measurement, DIM>& data() noexcept { 
                    
                    return data_;
                
                }


                /**
                 * @brief Get the uvector object
                 * 
                 * @return constexpr uvector 
                 */
                constexpr uvector as_uvector() const noexcept { 
                    
                    return *this; 
                
                }


                /**
                 * @brief Get the norm of the uvector
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement norm() const noexcept { 

                    if constexpr (DIM == 1) 
                        return data_[0];

                    uvector squared = square(*this);

                    return sqrt(std::accumulate(squared.data().begin(), squared.data().end(), uncertain_measurement{0.0, 0.0, squared.units()}));

                }


                /**
                 * @brief Get the squared norm of the uvector
                 * 
                 * @return constexpr uncertain_measurement 
                 */
                constexpr uncertain_measurement norm2() const noexcept { 

                    if constexpr (DIM == 1) 
                        return data_[0]; 

                    uvector squared(square(*this));

                    return std::accumulate(squared.data().begin(), squared.data().end(), uncertain_measurement{0.0, 0.0, squared.units()});
                    
                }


                /**
                 * @brief Get the normalization of the uvector
                 * 
                 * @return constexpr uvector 
                 */
                constexpr uvector versor() const {

                    uncertain_measurement norm = this->norm();                    
                    uvector result;
                    for (size_t i{}; i < DIM; ++i) 
                        result[i] = data_[i] / norm;

                    return result; 

                } 


                /**
                 * @brief Get the polar angle
                 * 
                 * @note the uvector must have at least two elements
                 * 
                 * @return constexpr uncertain_measurement
                 */
                constexpr uncertain_measurement phi() const noexcept requires (DIM >= 2) { 
                    
                    return atan(data_[1] / data_[0]); 
                
                }     
                

                /**
                 * @brief Get the azimuthal angle
                 * 
                 * @note the uvector must have at least three elements
                 * 
                 * @return constexpr uncertain_measurement
                 */
                constexpr uncertain_measurement theta() const requires (DIM >= 3) { 

                    if (data_[2] == 0.0 * m) 
                        return uncertain_measurement(0.0, data_[2].uncertainty_, rad);
                    else    
                        return acos((data_[2] / norm()));

                }
                

                /// @brief Print the uvector to the standard output
                constexpr void print(const bool& newline) const noexcept {

                    std::cout << "{\n";
                    for (size_t i{}; i < DIM; ++i) 
                        std::cout << "\t" << data_[i] << "\n";
                    std::cout << "}";
                    if (newline) std::cout << "\n";

                }   


                /**
                 * @brief Save the uvector to a file
                 * 
                 * @param file_name: the name of the file
                 * @param units: desired units for the output
                 */
                void save(const std::string& file_name, const unit& units) const {

                    std::ofstream file_out(file_name, std::ios::app);
                    if (file_out.is_open()) 

                        for (size_t i{}; i < DIM; ++i) {

                            file_out << data_[i].value_as(units); 
                            if (i < DIM - 1) 
                                file_out << ' '; 

                        }

                    else throw std::invalid_argument("Unable to open '" + file_name + "'");

                    file_out << '\n';
                    file_out.close();

                }


            protected:

            // =============================================
            // class members
            // =============================================

                std::array<uncertain_measurement, DIM> data_; ///< array of uncertain_measurements


            // =============================================
            // friends
            // =============================================

                friend class vector<DIM>; 
                 

        }; // class uvector


        /**
         * @brief Class expressing a generic matrix of vectors of measurements in a n-dimentional system
         * @see vector
         * 
         * @tparam rows: size_t number of rows
         * @tparam cols: size_t number of columns
         * 
         * @note the matrix is stored in row-major order
         * @note cols is set equal to rows by default 
         */
        template <size_t rows, size_t cols = rows> 
        class matrix {


            public:

                // =============================================
                // class members
                // =============================================

                size_t rows_; 

                size_t cols_;

                using mat = std::array<vector<rows>, cols>;

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


                constexpr matrix operator*(const scalar& scalar) const noexcept {
                    
                    assert(scalar != 0.0); 
                    matrix result;
                    for (size_t i{}; i < cols_; i++) result.data_[i] = data_[i] * scalar;
                    return result;

                }


                constexpr matrix operator/(const scalar& scalar) const noexcept {
                    
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

                /// @brief Construct a new time object as default = measurement(0.0, s)
                explicit constexpr time() noexcept : 

                    time_(0.0, s) {}


                /**
                 * @brief Construct a new time object from value and units of time
                 * 
                 * @param time: scalar value
                 * @param unit: unit as l-value const reference
                 * 
                 * @note the unit.base() must be second
                 */
                explicit constexpr time(const scalar& time, const unit& unit = s) {

                    if (time < 0.0) 
                        throw std::invalid_argument("Time cannot be negative");
                    if (unit.base() != base::second) 
                        throw std::invalid_argument("Wrong time unit, the unit_base must be seconds");
                    time_ = measurement(time, unit);
                    
                }


                /**
                 * @brief Construct a new time object from a time measurement
                 * 
                 * @param time: measurement as l-value const reference
                 * 
                 * @note time.units().base() must be second
                 */
                constexpr time(const measurement& time) {

                    if (time.value() < 0.0) 
                        throw std::invalid_argument("Time cannot be negative");
                    if (time.units().base() != base::second) 
                        throw std::invalid_argument("Wrong time unit, the unit_base must be seconds");
                    time_ = time;

                }


                /**
                 * @brief Construct a new time object from a time measurement
                 * 
                 * @param time: measurement as r-value reference
                 * 
                 * @note time.units().base() must be second
                 */
                constexpr time(measurement&& time) {

                    if (time.value() < 0.0) 
                        throw std::invalid_argument("Time cannot be negative");
                    if (time.units().base() != base::second) 
                        throw std::invalid_argument("Wrong time unit, the unit_base must be seconds");
                    time_ = std::move(time);

                }


                /**
                 * @brief Copy construct a new time object from another time
                 * 
                 * @param other: time object to copy from as l-value const reference
                 */
                constexpr time(const time& other) noexcept :
                    
                    time_{other.time_} {}   

                
                /**
                 * @brief Move construct a new time object from another time
                 * 
                 * @param other: time object to move from as r-value reference
                 */
                constexpr time(time&& other) noexcept :
                    
                    time_{std::move(other.time_)} {}
                

                /// @brief Default destructor
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
         * @see template <size_t> vector 
         * 
         * @note The position unit_base must be metre
         * @note The position is expressed in a cartesian coordinate system
         * @todo Add more coordinate system (e.g. polar, spherical, cylindrical)
         * 
         * @tparam DIM: size_t dimentions of the space
         */
        template <size_t DIM> 
        class position : public vector<DIM> {


            public: 

            // =============================================
            // constructors & destructor
            // =============================================

                /**
                 * @brief Construct a new default position object
                 * 
                 * @note The position is set to zero metres
                 */
                explicit constexpr position() noexcept : 
                
                    vector<DIM>(base::metre) {}


                /**
                 * @brief Construct a new position object from an std::array of measurements
                 * 
                 * @param pos: std::array of measurements as l-value const reference
                 */
                constexpr position(const std::array<measurement, DIM>& pos) noexcept : 
                    
                    vector<DIM>(pos) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (pos[i].units().base() != base::metre) 
                                throw std::invalid_argument("Wrong position unit, the unit_base must be metre");

                        }
                    
                    }


                /**
                 * @brief Construct a new position object from an std::array of measurements
                 * 
                 * @param pos: std::array of measurements as r-value reference
                 */
                constexpr position(std::array<measurement, DIM>&& pos) noexcept : 
                    
                    vector<DIM>(std::move(pos)) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (pos[i].units().base() != base::metre) 
                                throw std::invalid_argument("Wrong position unit, the unit_base must be metres");

                        }
                    
                    }


                /**
                 * @brief Construct a new position object from a vector of measurements
                 * 
                 * @param pos: vector of measurements as l-value const reference
                 */
                constexpr position(const vector<DIM>& pos) noexcept : 
                    
                    vector<DIM>(pos) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (pos[i].units().base() != base::metre) 
                                throw std::invalid_argument("Wrong position unit, the unit_base must be metres");
                        
                        }
                    
                    }


                /**
                 * @brief Construct a new position object from a vector of measurements
                 * 
                 * @param pos: vector of measurements as r-value reference
                 */
                constexpr position(vector<DIM>&& pos) noexcept : 
                    
                    vector<DIM>(std::move(pos)) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (pos[i].units().base() != base::metre) 
                                throw std::invalid_argument("Wrong position unit, the unit_base must be metres");
                        
                        }
                    
                    }

                
                /// @brief Default destructor
                ~position() = default; 


            // =============================================
            // operators
            // =============================================

                /**
                 * @brief Output stream operator
                 * 
                 * @param os: std::ostream&
                 * @param pos: position to be printed as l-value const reference
                 * 
                 * @return constexpr std::ostream& 
                 */
                friend constexpr std::ostream& operator<<(std::ostream& os, const position& pos) noexcept { 

                    os << "position = " << pos.as_vector(); 

                    return os;

                }


                /**
                 * @brief Output stream operator
                 * 
                 * @param file: std::ofstream&
                 * @param pos: position to be printed as l-value const reference
                 * 
                 * @return constexpr std::ofstream& 
                 */
                friend constexpr std::ofstream& operator<<(std::ofstream& file, const position& pos) noexcept { 

                    file << "position = " << pos.as_vector(); 

                    return file;

                }


            // =============================================
            // get methods
            // =============================================

                /**
                 * @brief Get the distance between this position and another position
                 * 
                 * @param other: position to be compared as l-value const reference
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement distance(const position& other) const noexcept {    
                    
                    return (other - *this).norm();
                    
                }       


                /**
                 * @brief Get the square distance between this position and another position
                 * 
                 * @param other: position to be compared as l-value const reference
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement distance2(const position& other) const noexcept {    
                    
                    return (other - *this).norm2();
                    
                }       


                /**
                 * @brief Get the polar angle of this position
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement phi() const noexcept {

                    return vector<DIM>::phi();

                }

                    
                /**
                 * @brief Get the polar angle between this position and another position
                 * 
                 * @param other: position to be compared as l-value const reference
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement phi(const position& other) const noexcept requires (DIM >= 2) { 

                    if (other != *this) 
                        return (other - *this).phi();
                    else 
                        return 0.0 * rad;
                                        
                }


                /**
                 * @brief Get the azimuthal angle of this position
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement theta() const noexcept {

                    return vector<DIM>::theta();

                }


                /**
                 * @brief Get the azimuthal angle between this position and another position
                 * 
                 * @param other: position to be compared as l-value const reference
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement theta(const position& other) const noexcept requires (DIM >= 3) {

                    if (other.pos_[2] != this->operator[](2)) 
                        return (other - *this).theta(); 
                    else 
                        return 0.0 * rad;
                
                }


                /**
                 * @brief Get the versor between this position and another position
                 * 
                 * @param other: position to be compared as l-value const reference
                 * 
                 * @return constexpr verctor
                 */
                constexpr vector<DIM> direction(const position& other) const noexcept requires (DIM >= 2) {

                    if (other != *this)
                        return vector<DIM>((other - *this) / this->distance(other));
                    else 
                        return vector<DIM>(0 * unitless);

                } 


                /**
                 * @brief Print the position to the standard output
                 * 
                 * @param newline: if true, print a newline character at the end of the output
                 * 
                 * @return void
                 */
                constexpr void print(const bool& newline = true) const noexcept {

                    std::cout << "position = "; 
                    this->vector<DIM>::print(newline); 

                }   


        }; // class position


        /**
         * @brief Class expressing the linear velocity of a generic object as a vector of measurements
         * @see template <size_t> vector 
         * 
         * @note The linear_velocity unit_base must be metre / second
         * @note The linear_velocity is expressed in a cartesian coordinate system
         * @todo Add more coordinate system (e.g. polar, spherical, cylindrical)
         * 
         * @tparam DIM: size_t dimentions of the space
         */
        template <size_t DIM> 
        class linear_velocity : public vector<DIM> {

            
            public: 

            // =============================================
            // constructors & destructor
            // =============================================

                /**
                 * @brief Construct a new default linear_velocity object
                 * 
                 * @note The linear_velocity is set to zero metre / second
                 */
                explicit constexpr linear_velocity() noexcept : 
                
                    vector<DIM>(base::metre / base::second) {}


                /**
                 * @brief Construct a new linear_velocity object from an std::array of measurements
                 * 
                 * @param vel: std::array of measurements as l-value const reference
                 */
                constexpr linear_velocity(const std::array<measurement, DIM>& vel) noexcept : 
                    
                    vector<DIM>(vel) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (vel[i].units().base() != base::metre / base::second) 
                                throw std::invalid_argument("Wrong linear_velocity unit, the unit_base must be metre / second");

                        }
                    
                    }


                /**
                 * @brief Construct a new linear_velocity object from an std::array of measurements
                 * 
                 * @param vel: std::array of measurements as r-value reference
                 */
                constexpr linear_velocity(std::array<measurement, DIM>&& vel) noexcept : 
                    
                    vector<DIM>(std::move(vel)) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (vel[i].units().base() != base::metre / base::second) 
                                throw std::invalid_argument("Wrong linear_velocity unit, the unit_base must be metre / second");

                        }
                    
                    }


                /**
                 * @brief Construct a new linear_velocity object from a vector of measurements
                 * 
                 * @param vel: vector of measurements as l-value const reference
                 */
                constexpr linear_velocity(const vector<DIM>& vel) noexcept : 
                    
                    vector<DIM>(vel) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (vel[i].units().base() != base::metre / base::second) 
                                throw std::invalid_argument("Wrong linear_velocity unit, the unit_base must be metre / second");
                        
                        }
                    
                    }


                /**
                 * @brief Construct a new linear_velocity object from a vector of measurements
                 * 
                 * @param vel: vector of measurements as r-value reference
                 */
                constexpr linear_velocity(vector<DIM>&& vel) noexcept : 
                    
                    vector<DIM>(std::move(vel)) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (vel[i].units().base() != base::metre / base::second) 
                                throw std::invalid_argument("Wrong linear_velocity unit, the unit_base must be metre / second");
                        
                        }
                    
                    }

                
                /// @brief Default destructor
                ~linear_velocity() = default; 


            // =============================================
            // operators
            // =============================================

                /**
                 * @brief Output stream operator
                 * 
                 * @param os: std::ostream&
                 * @param pos: linear_velocity to be printed as l-value const reference
                 * 
                 * @return constexpr std::ostream& 
                 */
                friend constexpr std::ostream& operator<<(std::ostream& os, const linear_velocity& vel) noexcept { 

                    os << "linear velocity = " << vel.as_vector(); 

                    return os;

                }


                /**
                 * @brief Output stream operator
                 * 
                 * @param file: std::ofstream&
                 * @param vel: linear_velocity to be printed as l-value const reference
                 * 
                 * @return constexpr std::ofstream& 
                 */
                friend constexpr std::ofstream& operator<<(std::ofstream& file, const linear_velocity& vel) noexcept { 

                    file << "linear velocity = " << vel.as_vector(); 

                    return file;

                }


            // =============================================
            // print
            // =============================================

                /**
                 * @brief Print the linear_velocity to the standard output
                 * 
                 * @param newline: if true, print a newline character at the end of the output
                 * 
                 * @return void
                 */
                constexpr void print(const bool& newline = true) const {

                    std::cout << "linear velocity = ";
                    vector<DIM>::print(newline);

                }   


        }; // class linear_velocity
    

        /**
         * @brief Class expressing the linear acceleration of a generic object as a vector of measurements
         * @see template <size_t> vector 
         * 
         * @note The linear_acceleration unit_base must be metre / second.square()
         * @note The linear_acceleration is expressed in a cartesian coordinate system
         * @todo Add more coordinate system (e.g. polar, spherical, cylindrical)
         * 
         * @tparam DIM: size_t dimentions of the space
         */
        template <size_t DIM> 
        class linear_acceleration : public vector<DIM> {

            
            public: 

            // =============================================
            // constructors & destructor
            // =============================================

                /**
                 * @brief Construct a new default linear_acceleration object
                 * 
                 * @note The linear_acceleration is set to zero metre / second.square()
                 */
                explicit constexpr linear_acceleration() noexcept : 
                
                    vector<DIM>(base::metre / base::second.square()) {}


                /**
                 * @brief Construct a new linear_acceleration object from an std::array of measurements
                 * 
                 * @param acc: std::array of measurements as l-value const reference
                 */
                constexpr linear_acceleration(const std::array<measurement, DIM>& acc) noexcept : 
                    
                    vector<DIM>(acc) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (acc[i].units().base() != base::metre / base::second.square()) 
                                throw std::invalid_argument("Wrong linear_acceleration unit, the unit_base must be metre / second.square()");

                        }
                    
                    }


                /**
                 * @brief Construct a new linear_acceleration object from an std::array of measurements
                 * 
                 * @param acc: std::array of measurements as r-value reference
                 */
                constexpr linear_acceleration(std::array<measurement, DIM>&& acc) noexcept : 
                    
                    vector<DIM>(std::move(acc)) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (acc[i].units().base() != base::metre / base::second.square()) 
                                throw std::invalid_argument("Wrong linear_acceleration unit, the unit_base must be metre / second.square()");

                        }
                    
                    }


                /**
                 * @brief Construct a new linear_acceleration object from a vector of measurements
                 * 
                 * @param acc: vector of measurements as l-value const reference
                 */
                constexpr linear_acceleration(const vector<DIM>& acc) noexcept : 
                    
                    vector<DIM>(acc) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (acc[i].units().base() != base::metre / base::second.square()) 
                                throw std::invalid_argument("Wrong linear_acceleration unit, the unit_base must be metre / second.square()");
                        
                        }
                    
                    }


                /**
                 * @brief Construct a new linear_acceleration object from a vector of measurements
                 * 
                 * @param acc: vector of measurements as r-value reference
                 */
                constexpr linear_acceleration(vector<DIM>&& acc) noexcept : 
                    
                    vector<DIM>(std::move(acc)) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (acc[i].units().base() != base::metre / base::second.square()) 
                                throw std::invalid_argument("Wrong linear_acceleration unit, the unit_base must be metre / second.square()");
                        
                        }
                    
                    }

                
                /// @brief Default destructor
                ~linear_acceleration() = default; 


            // =============================================
            // operators
            // =============================================

                /**
                 * @brief Output stream operator
                 * 
                 * @param os: std::ostream&
                 * @param pos: linear_acceleration to be printed as l-value const reference
                 * 
                 * @return constexpr std::ostream& 
                 */
                friend constexpr std::ostream& operator<<(std::ostream& os, const linear_acceleration& acc) noexcept { 

                    os << "linear acceleration = " << acc.as_vector(); 

                    return os;

                }


                /**
                 * @brief Output stream operator
                 * 
                 * @param file: std::ofstream&
                 * @param acc: linear_acceleration to be printed as l-value const reference
                 * 
                 * @return constexpr std::ofstream& 
                 */
                friend constexpr std::ofstream& operator<<(std::ofstream& file, const linear_acceleration& acc) noexcept { 

                    file << "linear acceleration = " << acc.as_vector(); 

                    return file;

                }


            // =============================================
            // print
            // =============================================

                /**
                 * @brief Print the linear_acceleration to the standard output
                 * 
                 * @param newline: if true, print a newline character at the end of the output
                 * 
                 * @return void
                 */
                constexpr void print(const bool& newline = true) const {

                    std::cout << "linear acceleration = ";
                    vector<DIM>::print(newline);

                }   


        }; // class linear_acceleration


            // /**
            //  * @brief Create a new position vector multiplying a linear velocity with a time measurement
            //  * 
            //  * @tparam DIM: dimentions of the space
            //  * @param vel: linear_velocity as l-value const reference
            //  * @param t: time as l-value const reference 
            //  * 
            //  * @return constexpr position<DIM> 
            //  */
            // template <size_t DIM> 
            // constexpr position<DIM> operator*(const linear_velocity<DIM>& vel, const time& t) noexcept {

            //     return vel * t.as_measurement();

            // }


            // /**
            //  * @brief Create a new position vector multiplying a linear velocity with a time measurement
            //  * 
            //  * @tparam DIM: dimentions of the space
            //  * @param vel: linear_velocity as r-value reference
            //  * @param t: time as r-value reference 
            //  * 
            //  * @return constexpr position<DIM> 
            //  */
            // template <size_t DIM> 
            // constexpr position<DIM> operator*(linear_velocity<DIM>&& vel, time&& t) noexcept {

            //     return vel * t.as_measurement();

            // }


            // /**
            //  * @brief Create a new position vector multiplying a time measurement with a linear velocity 
            //  * 
            //  * @tparam DIM: dimentions of the space
            //  * @param t: time as l-value const reference 
            //  * @param vel: linear_velocity as l-value const reference
            //  * 
            //  * @return constexpr position<DIM> 
            //  */
            // template <size_t DIM> 
            // constexpr position<DIM> operator*(const time& t, const linear_velocity<DIM>& vel) noexcept {

            //     return vel * t.as_measurement();

            // }


            // /**
            // * @brief Create a new position vector multiplying a linear velocity with a time measurement
            // * 
            // * @tparam DIM: dimentions of the space
            // * @param t: time as r-value reference 
            // * @param vel: linear_velocity as r-value reference
            // * 
            // * @return constexpr position<DIM> 
            // */
            // template <size_t DIM> 
            // constexpr position<DIM> operator*(time&& t, linear_velocity<DIM>&& vel) noexcept {

            //     return vel * t.as_measurement();

            // }


            // /**
            //  * @brief Create a new linear_velocity object dividing a position by time 
            //  * 
            //  * @tparam DIM: dimentions of the space
            //  * @param pos: position<DIM> as l-value const reference
            //  * @param t: time as l-value const reference
            //  * 
            //  * @return linear_velocity<DIM>
            //  */
            // template <size_t DIM> 
            // constexpr linear_velocity<DIM> operator/(const position<DIM>& pos, const time& t) {

            //     return pos / t.as_measurement();

            // }


            // /**
            //  * @brief Create a new linear_velocity object dividing a position by time 
            //  * 
            //  * @tparam DIM: dimentions of the space
            //  * @param pos: position<DIM> as r-value reference
            //  * @param t: time as r-value reference
            //  * 
            //  * @return linear_velocity<DIM>
            //  */
            // template <size_t DIM> 
            // constexpr linear_velocity<DIM> operator/(position<DIM>&& pos, time&& t) {

            //     return pos / t.as_measurement();

            // }


            // /**
            //  * @brief Create a new linear_velocity object multiplying a linear_acceleration with a time measurement
            //  * 
            //  * @tparam DIM: dimentions of the space
            //  * @param acc: linear_acceleration as l-value const reference
            //  * @param t: time as l-value const reference 
            //  * 
            //  * @return constexpr linear_velocity<DIM> 
            //  */
            // template <size_t DIM> 
            // constexpr linear_velocity<DIM> operator*(const linear_acceleration<DIM>& acc, const time& t) noexcept {

            //     return acc * t.as_measurement();

            // }


            // /**
            //  * @brief Create a new linear_velocity object multiplying a linear_acceleration with a time measurement
            //  * 
            //  * @tparam DIM: dimentions of the space
            //  * @param acc: linear_acceleration<DIM> as r-value reference
            //  * @param t: time as r-value reference
            //  * 
            //  * @return linear_velocity<DIM>
            //  */
            // template <size_t DIM> 
            // constexpr linear_velocity<DIM> operator*(linear_acceleration<DIM>&& acc, time&& t) noexcept {

            //     return acc * t.as_measurement();

            // }


            // /**
            //  * @brief Create a new linear_velocity object multiplying a time measurement with a linear_acceleration
            //  * 
            //  * @tparam DIM: dimentions of the space
            //  * @param t: time as l-value const reference
            //  * @param acc: linear_acceleration as l-value const reference
            //  * 
            //  * @return linear_velocity<DIM>
            //  */
            // template <size_t DIM> 
            // constexpr linear_velocity<DIM> operator*(const time& t, const linear_acceleration<DIM>& acc) noexcept {

            //     return acc * t.as_measurement();

            // }


            // /**
            //  * @brief Create a new linear_velocity object multiplying a linear_acceleration with a time measurement
            //  * 
            //  * @tparam DIM: dimentions of the space
            //  * @param t: time as r-value reference
            //  * @param acc: linear_acceleration<DIM> as r-value reference
            //  * 
            //  * @return linear_velocity<DIM>
            //  */
            // template <size_t DIM> 
            // constexpr linear_velocity<DIM> operator*(time&& t, linear_acceleration<DIM>&& acc) noexcept {

            //     return acc * t.as_measurement();

            // }


            // /**
            //  * @brief Create a new linear_acceleration object dividing a linear_velocity by time 
            //  * 
            //  * @tparam DIM: dimentions of the space
            //  * @param pos: linear_velocity<DIM> as l-value const reference
            //  * @param t: time as l-value const reference
            //  * 
            //  * @return linear_acceleration<DIM>
            //  */        
            // template <size_t DIM> 
            // constexpr linear_acceleration<DIM> operator/(const linear_velocity<DIM>& vel, const time& t) {

            //     return vel / t.as_measurement();

            // }   


            // /**
            //  * @brief Create a new linear_acceleration object dividing a linear_velocity by time 
            //  * 
            //  * @tparam DIM: dimentions of the space
            //  * @param pos: linear_velocity<DIM> as r-value reference
            //  * @param t: time as r-value reference
            //  * 
            //  * @return linear_acceleration<DIM>
            //  */       
            // template <size_t DIM> 
            // constexpr linear_acceleration<DIM> operator/(linear_velocity<DIM>&& vel, time&& t) {

            //     return vel / t.as_measurement();

            // }   


        /**
         * @brief Class expressing the force as a vector of measurements
         * @see template <size_t> vector 
         * 
         * @note The force unit_base must be kilogram * kilogram * metre / second.square()
         * @note The force is expressed in a cartesian coordinate system
         * @todo Add more coordinate system (e.g. polar, spherical, cylindrical)
         * 
         * @tparam DIM: size_t dimentions of the space
         */
        template <size_t DIM> 
        class force : public vector<DIM> {

            
            public: 

            // =============================================
            // constructors & destructor
            // =============================================

                /**
                 * @brief Construct a new default force object
                 * 
                 * @note The force is set to zero kilogram * metre / second.square()
                 */
                explicit constexpr force() noexcept : 
                
                    vector<DIM>(base::kilogram * base::metre / base::second.square()) {}


                /**
                 * @brief Construct a new force object from an std::array of measurements
                 * 
                 * @param F: std::array of measurements as l-value const reference
                 */
                constexpr force(const std::array<measurement, DIM>& F) noexcept : 
                    
                    vector<DIM>(F) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (F[i].units().base() != base::kilogram * base::metre / base::second.square()) 
                                throw std::invalid_argument("Wrong force unit, the unit_base must be kilogram * metre / second.square()");

                        }
                    
                    }


                /**
                 * @brief Construct a new force object from an std::array of measurements
                 * 
                 * @param F: std::array of measurements as r-value reference
                 */
                constexpr force(std::array<measurement, DIM>&& F) noexcept : 
                    
                    vector<DIM>(std::move(F)) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (F[i].units().base() != base::kilogram * base::metre / base::second.square()) 
                                throw std::invalid_argument("Wrong force unit, the unit_base must be kilogram * metre / second.square()");

                        }
                    
                    }


                /**
                 * @brief Construct a new force object from a vector of measurements
                 * 
                 * @param F: vector of measurements as l-value const reference
                 */
                constexpr force(const vector<DIM>& F) noexcept : 
                    
                    vector<DIM>(F) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (F[i].units().base() != base::kilogram * base::metre / base::second.square()) 
                                throw std::invalid_argument("Wrong force unit, the unit_base must be kilogram * metre / second.square()");
                        
                        }
                    
                    }


                /**
                 * @brief Construct a new force object from a vector of measurements
                 * 
                 * @param F: vector of measurements as r-value reference
                 */
                constexpr force(vector<DIM>&& F) noexcept : 
                    
                    vector<DIM>(std::move(F)) {

                        for (size_t i{}; i < DIM; ++i) {

                            if (F[i].units().base() != base::kilogram * base::metre / base::second.square()) 
                                throw std::invalid_argument("Wrong force unit, the unit_base must be kilogram * metre / second.square()");
                        
                        }
                    
                    }

                
                /// @brief Default destructor
                ~force() = default; 


            // =============================================
            // print methods
            // =============================================

                /**
                 * @brief Output stream operator
                 * 
                 * @param os: std::ostream&
                 * @param pos: force to be printed as l-value const reference
                 * 
                 * @return constexpr std::ostream& 
                 */
                friend constexpr std::ostream& operator<<(std::ostream& os, const force& F) noexcept { 

                    os << "force = " << F.as_vector(); 

                    return os;

                }


                /**
                 * @brief Output stream operator
                 * 
                 * @param file: std::ofstream&
                 * @param F: force to be printed as l-value const reference
                 * 
                 * @return constexpr std::ofstream& 
                 */
                friend constexpr std::ofstream& operator<<(std::ofstream& file, const force& F) noexcept { 

                    file << "force = " << F.as_vector(); 

                    return file;

                }


                /**
                 * @brief Print the force to the standard output
                 * 
                 * @param newline: if true, print a newline character at the end of the output
                 * 
                 * @return void
                 */
                constexpr void print(const bool& newline = true) const {

                    std::cout << "force = ";
                    vector<DIM>::print(newline);

                }   


        }; // class force


        /// @brief Class for timing the execution of a generic function/process
        class timer {

            
            public:

                // =============================================
                // constructor and destructor
                // =============================================   

                /**
                 * @brief Construct a new timer object
                 */
                constexpr timer() = default;


                /**
                 * @brief Destroy the timer object
                 */
                ~timer() = default;

                
                // =============================================
                // timer methods
                // =============================================   

                /**
                 * @brief start the timer
                 */
                inline void start() { 
                    
                    start_ = std::chrono::high_resolution_clock::now(); 
                    
                }


                /**
                 * @brief stop the timer
                 */
                inline void pause() { 
                    
                    pause_ = std::chrono::high_resolution_clock::now(); 
                    
                }
                

                constexpr measurement elapsed(const unit& units = s) const {

                    if (units != s) return measurement(s.convert(static_cast<std::chrono::duration<scalar>>(pause_ - start_).count(), units), units);
                    else return measurement(static_cast<std::chrono::duration<scalar>>(pause_ - start_).count(), units);

                }


                /**
                 * @brief print the elapsed time
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


    /// @brief Namespace constains some physical constants
    namespace constants {


        using namespace physics::measurements; 


        constexpr measurement G(6.67430e-17, N * km.square() / kg.square()); ///< gravitational constant
        
        constexpr measurement c(299792458, m / s); ///< speed of light in vacuum

        constexpr measurement h(6.62607015e-34, J * s); ///< Planck constant

        constexpr measurement hbar(1.054571817e-34, J * s); ///< reduced Planck constant

        constexpr measurement mu0(1.25663706212e-6, N / A.square()); ///< magnetic permeability of vacuum

        constexpr measurement eps0 = 1 / (mu0 * square(c)); ///< electric permittivity of vacuum

        constexpr measurement e(1.602176634e-19, C); ///< elementary charge

        constexpr measurement m_e(9.1093837015e-31, kg); ///< electron mass

        constexpr measurement m_p(1.67262192369e-27, kg); ///< proton mass

        constexpr measurement m_n(1.67492749804e-27, kg); ///< neutron mass

        constexpr measurement k_B(1.380649e-23, J / K); ///< Boltzmann constant

        constexpr measurement N_A(6.02214076e23, mol.inv()); ///< Avogadro constant

        constexpr measurement R = N_A * k_B; ///< ideal gas constant


    } // namespace constants


    /// @brief Namespace contains some objects 
    namespace objects {


        using namespace physics::tools;

        
        /**
         * @brief Class for a generic mass_point
         * 
         * @tparam DIM: dimension of the space
         */
        template <size_t DIM>
        class mass {


            public: 

                // =============================================
                // constructors and destructor 
                // =============================================

                /**
                 * @brief Construct a new mass object
                 * 
                 * @param m: measurement of mass
                 * @param pos: position<DIM> vector
                 * @param vel: linear_velocity vector
                 * @param acc: linear_acceleration vector
                 * @param gravity: bool (if true, the object is affected by the gravitational field)
                 * @param id: uint32_t id
                 * @param type: std::string type
                 * @param name: std::string name
                 */
                explicit constexpr mass(const measurement& m = measurement(0.0, kg), 
                                        const position<DIM>& pos = position<DIM>(), 
                                        const linear_velocity<DIM>& vel = linear_velocity<DIM>(), 
                                        const linear_acceleration<DIM>& acc = linear_acceleration<DIM>(), 
                                        const bool& gravity = true,
                                        const uint32_t& id = 0, 
                                        const std::string& type = nullptr,
                                        const std::string& name = nullptr
                                        // shape shape = geometry::shapes::circle(0.0 * m)
                                        ) : 

                    // shape_(std::move(&shape)), 
                    mass_(m),
                    position_(pos),
                    linear_velocity_(vel), 
                    linear_acceleration_(acc), 
                    gravitational_field_{gravity},
                    id_{id}, 
                    type_{type},
                    name_{name} {

                        assert(mass_ >= 0.0 * kg);

                    }


                /**
                 * @brief Copy construct a new mass object
                 * 
                 * @param other: mass object to copy from as l-value const reference
                 */
                constexpr mass(const mass& other) : 
                

                    // shape_{other.shape_},
                    mass_(other.mass_),
                    position_(other.position_), 
                    linear_velocity_(other.linear_velocity_), 
                    linear_acceleration_(other.linear_acceleration_),
                    gravitational_field_{other.gravitational_field_},
                    id_{other.id_}, 
                    type_{other.type_},
                    name_{other.name_} {

                        std::cout << "object copy constructor called\n";

                    } 


                /**
                 * @brief Move construct a new mass object
                 * 
                 * @param other: mass object to move from as r-value reference
                 */
                constexpr mass(mass&& other) : 
                
                    // shape_{std::move(other.shape_)},
                    mass_(std::move(other.mass_)),
                    position_{std::move(other.position_)}, 
                    linear_velocity_{std::move(other.linear_velocity_)}, 
                    linear_acceleration_{std::move(other.linear_acceleration_)},
                    gravitational_field_{other.gravitational_field_},
                    id_{std::move(other.id_)}, 
                    type_{std::move(other.type_)},
                    name_{std::move(other.name_)} {

                        std::cout << "object move constructor called\n";

                    } 


                /// @brief Default destructor
                ~mass() noexcept = default;


            // =============================================
            // operators
            // =============================================

                /**
                 * @brief Copy assignment operator
                 * 
                 * @param other: mass object to copy from as l-value const reference
                 * 
                 * @return constexpr mass& 
                 */
                constexpr mass& operator=(const mass& other) noexcept {

                    std::cout << "object copy assignment operator called\n";

                    if (*this != other) {

                        id_ = other.id_;
                        type_ = other.type_;
                        name_ = other.name_;
                        // shape_ = other.shape_;
                        position_ = other.position_;
                        linear_velocity_ = other.linear_velocity_;
                        linear_acceleration_ = other.linear_acceleration_;

                    }

                    return *this;

                }


                /**
                 * @brief Move assignment operator
                 * 
                 * @param other: mass object to move from as r-value reference
                 * 
                 * @return constexpr mass& 
                 */
                constexpr mass& operator=(mass&& other) noexcept {

                    std::cout << "object move assignment operator called\n";

                    if (*this != other) {

                        // shape_ = std::move(other.shape_);
                        mass_ = std::move(other.mass_); 
                        position_ = std::move(other.position_);
                        linear_velocity_ = std::move(other.linear_velocity_);
                        linear_acceleration_ = std::move(other.linear_acceleration_);
                        gravitational_field_ = std::move(other.gravitational_field_);
                        id_ = std::move(other.id_);
                        type_ = std::move(other.type_);
                        name_ = std::move(other.name_);

                    }

                    return *this;

                }


                /**
                 * @brief Equality operator
                 * 
                 * @param other: mass object to compare with as l-value const reference
                 * 
                 * @return bool
                 */
                constexpr bool operator==(const mass& other) const noexcept {

                    return (id_ == other.id_ && 
                            type_ == other.type_ && 
                            name_ == other.name_ && 
                            // shape_ == other.shape_ && 
                            position_ == other.position_ && 
                            linear_velocity_ == other.linear_velocity_ && 
                            linear_acceleration_ == other.linear_acceleration_);

                }


                /**
                 * @brief Inequality operator
                 * 
                 * @param other: mass object to !compare with as l-value const reference
                 * 
                 * @return bool
                 */
                constexpr bool operator!=(const mass& other) const noexcept {

                    return !(*this == other);

                }


            // =============================================
            // object methods
            // =============================================


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
            // get methods
            // =============================================

                /**
                 * @brief Get the mass measurement 
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement mass_measurement() const noexcept { 
                    
                    return mass_; 
                    
                }


                /**
                 * @brief Get the mass measurement
                 * 
                 * @return constexpr measurement& 
                 */
                constexpr measurement& mass_measurement() noexcept {

                    return mass_; 

                }


                /**
                 * @brief Get the position vector 
                 * 
                 * @return constexpr position<DIM> 
                 */
                constexpr position<DIM> as_position() const noexcept { 
                    
                    return position_; 
                
                }


                /**
                 * @brief Get the position vector
                 * 
                 * @return constexpr position<DIM>& 
                 */
                constexpr position<DIM>& as_position() noexcept { 
                    
                    return position_; 
                
                }


                /**
                 * @brief Get the linear velocity vector 
                 * 
                 * @return constexpr linear_velocity<DIM> 
                 */
                constexpr linear_velocity<DIM> as_linear_velocity() const noexcept { 
                    
                    return linear_velocity_; 
                
                }


                /**
                 * @brief Get the linear velocity vector
                 * 
                 * @return constexpr linear_velocity<DIM>& 
                 */
                constexpr linear_velocity<DIM>& as_linear_velocity() noexcept { 
                    
                    return linear_velocity_; 
                
                }


                /**
                 * @brief Get the linear acceleration vector 
                 * 
                 * @return constexpr linear_acceleration<DIM> 
                 */
                constexpr linear_acceleration<DIM> as_linear_acceleration() const noexcept { 
                    
                    return linear_acceleration_; 
                
                }


                /**
                 * @brief Get the linear acceleration vector
                 * 
                 * @return constexpr linear_acceleration<DIM>& 
                 */
                constexpr linear_acceleration<DIM>& as_linear_acceleration() noexcept { 
                    
                    return linear_acceleration_; 
                
                }


                /**
                 * @brief Get the matrix of the state (position, linear_velocity, linear_acceleration)
                 * 
                 * @return constexpr matrix<DIM, 3> 
                 */
                constexpr matrix<DIM, 3> state() const noexcept {

                    return matrix<DIM, 2>({position_.as_vector(), linear_velocity_.as_vector(), linear_acceleration_.as_vector()});

                }


                /**
                 * @brief Get the momentum vector
                 * 
                 * @return constexpr vector<DIM> 
                 */
                constexpr vector<DIM> momentum() const { 
                    
                    return mass_ * linear_velocity_.as_vector();
                    
                }


                /**
                 * @brief Get the angular momentum vector
                 * 
                 * @return constexpr vector<DIM> 
                 */
                constexpr vector<DIM> angular_momentum() const { 
                    
                    return mass_ * cross(position_.as_vector(), linear_velocity_.as_vector()); 
                    
                }


                /**
                 * @brief Get the kinetic energy
                 * 
                 * @return constexpr measurement 
                 */
                constexpr measurement kinetic_energy() const { 
                    
                    return 0.5 * mass_ * linear_velocity_.norm2(); 
                    
                }


                /**
                 * @brief Get the gravitational field flag
                 * 
                 * @return bool 
                 */
                constexpr bool gravitational_field() const noexcept { 
                    
                    return gravitational_field_; 
                    
                }


                /**
                 * @brief Get the gravitational field flag
                 * 
                 * @return bool&
                 */
                constexpr bool& gravitational_field() noexcept { 
                    
                    return gravitational_field_; 
                    
                }


                /**
                 * @brief Get the id
                 * 
                 * @return constexpr uint32_t 
                 */
                constexpr uint32_t id() const noexcept { 
                    
                    return id_;
                    
                }


                /**
                 * @brief Get the id
                 * 
                 * @return constexpr uint32_t& 
                 */
                constexpr uint32_t& id() noexcept { 
                    
                    return id_;
                    
                }

                
                /**
                 * @brief Get the type
                 * 
                 * @return constexpr std::string
                 */
                constexpr std::string type() const noexcept { 
                    
                    return type_;
                    
                }

                
                /**
                 * @brief Get the type
                 * 
                 * @return constexpr std::string& 
                 */
                constexpr std::string& type() noexcept { 
                    
                    return type_;
                    
                }


                /**
                 * @brief Get the name
                 * 
                 * @return constexpr std::string 
                 */
                constexpr std::string name() const noexcept { 
                    
                    return name_;
                    
                }


                /**
                 * @brief Get the name
                 * 
                 * @return constexpr std::string& 
                 */
                constexpr std::string& name() noexcept { 
                    
                    return name_;
                    
                }


/*
                inline shape& shape() const noexcept { 
                    
                    return *shape_;
                
                }
*/


                /**
                 * @brief Get the dimention of the space
                 * 
                 * @return constexpr size_t 
                 */
                constexpr size_t dim() const noexcept {

                    return DIM; 

                }


                /**
                 * @brief Get the mass object
                 * 
                 * @return constexpr object<DIM> 
                 */
                constexpr mass<DIM> as_mass_object() const noexcept { 
                    
                    return *this; 
                    
                }


                /**
                 * @brief Get the mass object
                 * 
                 * @return constexpr object<DIM>& 
                 */
                constexpr mass<DIM>& as_mass_object() noexcept { 
                    
                    return *this; 
                    
                }


                /// @brief Print the mass object to the standard output
                constexpr void print() const noexcept {

                    std::cout << "\nobject:\n";
                    if (id_ != 0) 
                        std::cout << "id = " << id_ << "\n";
                    if (!type_.empty())
                        std::cout << "type = " << type_ << "\n";
                    if (!name_.empty()) 
                        std::cout << "name = " << name_ << "\n"; 
                    // shape_->print(); 
                    position_.print(); 
                    linear_velocity_.print(); 
                    linear_acceleration_.print();

                }               


            protected:
                    
            // =============================================
            // class members
            // =============================================
                
                measurement mass_; ///< object mass measurement

                position<DIM> position_; ///< object position vector

                linear_velocity<DIM> linear_velocity_; ///< object linear velocity vector

                linear_acceleration<DIM> linear_acceleration_; ///< object linear acceleration vector

                // shape* shape_; ///< object shape

                uint32_t id_; ///< object id

                std::string type_; ///< object type

                std::string name_; ///< object name

                bool gravitational_field_; ///< if true, the object is affected by the gravitational field


        }; // class mass        

    
    } // namespace objects


} // namespace physics


namespace math {
    

    // namespace defining some tools for math
    namespace tools {


        using namespace physics::tools;


        // namespace defining some descriptive statistic functions
        namespace descriptive_statistics {
                

            // using measurement = physics::measurements::measurement;


            double mean(const std::vector<double>& v) {

                if (v.size() == 0) 
                    throw std::invalid_argument("Can't operate a descriptive statistic funtion on an empty vector"); 
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
                for (auto x : v) accu += pow(x - average, 2); 
                return accu / v.size();

            }


            inline double sd(const std::vector<double>& v) { 
                
                return std::sqrt(variance(v)); 
                
            }


            // inline double sdom(const std::vector<double>& v) { 
                
            //     return std::sqrt(variance(v) / v.size()); 
                
            // }


            /**
             * @brief Compute the mean value of a vector of measurements
             * @note The uncertainty is computed as the standard deviation of mean (sdom)
             * 
             * @param vec: vector of measurements
             * 
             * @return uncertain_measurement
             */
            uncertain_measurement mean(const std::vector<measurement>& vec) {

                if (vec.size() == 0) 
                    throw std::invalid_argument("Can't operate a descriptive statistic funtion on an empty vector"); 
                
                measurement average = std::accumulate(vec.begin(), vec.end(), measurement(0., vec[0].units())) / vec.size();
                measurement sigma_sq = measurement(0., vec[0].units().square()); 
                for (auto x : vec) 
                    sigma_sq += square(x - average); 

                return uncertain_measurement(average, sqrt(sigma_sq / (vec.size() * (vec.size() - 1))));                 

            }


            /**
             * @brief Compute the mean value of a vector of uncertain_measurements.
             * @note The uncertainty is computed as the standard deviation of mean (sdom)
             * 
             * @param vec: vector of uncertain_measurements
             * 
             * @return uncertain_measurement
             */
            uncertain_measurement mean(const std::vector<uncertain_measurement>& vec) {

                if (vec.size() == 0) 
                    throw std::invalid_argument("Can't operate a descriptive statistic funtion on an empty vector"); 
                
                measurement average = (std::accumulate(vec.begin(), vec.end(), uncertain_measurement(0., 0., vec[0].units())) / vec.size()).as_measurement();
                measurement sigma_sq = measurement(0., vec[0].units().square()); 
                for (auto x : vec) 
                    sigma_sq += square(x - average).as_measurement(); 

                return uncertain_measurement(average, sqrt(sigma_sq / (vec.size() * (vec.size() - 1))));                 

            }


            /**
             * @brief Compute the median of a vector of measurements
             * 
             * @param vec: vector of measurements
             * 
             * @return measurement 
             */
            measurement median(const std::vector<measurement>& v) {

                if (v.size() == 0) 
                    throw std::invalid_argument("Can't operate a descriptive statistic funtion on an empty vector"); 

                std::vector<measurement> copy{v}; 
                if (std::is_sorted(copy.begin(), copy.end()) == false) 
                    std::sort(copy.begin(), copy.end());

                if (v.size() % 2 != 0) 
                    return copy[v.size() / 2];
                else 
                    return (copy[v.size() / 2] + copy[(v.size() / 2) - 1]) / 2; 

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
                for (auto x : vec) 
                    sigma_sq += square(x - average); 

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
                
                return sqrt(variance(vec)); 
            
            }


            /**
             * @brief Compute the standard error of the mean of a vector of measurements
             * 
             * @param vec: vector of measurements
             * 
             * @return measurement 
             */
            inline measurement sdom(const std::vector<measurement>& vec) { 
                
                return sqrt((variance(vec) / (vec.size() - 1))); 
            
            }

            
            /**
             * @brief Compute the weighted mean of a vector of measurements
             * 
             * @param vec 
             * 
             * @return uncertain_measurement 
             */
            uncertain_measurement wmean(const std::vector<uncertain_measurement>& vec) {

                if (vec.size() == 0) 
                    throw std::invalid_argument("Can't operate a descriptive statistic funtion on an empty vector"); 
                
                measurement weighted = measurement(0., vec[0].units().inv());
                measurement weights = measurement(0., vec[0].units().inv().square());
                for (auto& x : vec) {
                    weighted += x.as_measurement() * x.weight(); 
                    weights += x.weight();
                }

                return uncertain_measurement(weighted / weights, sqrt(weights.inv()));

            }


            /**
             * @brief Compute the weighted variance of a vector of measurements
             * 
             * @param vec 
             * 
             * @return measurement 
             */
            measurement wvariance(const std::vector<uncertain_measurement>& vec) {

                if (vec.size() == 0) 
                    throw std::invalid_argument("Can't operate a descriptive statistic funtion on an empty vector"); 
                
                measurement weights = square(measurement(0., vec[0].units().inv()));
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

                return sqrt(wvariance(vec));
            
            }


            measurement chi_sq(const std::vector<measurement>& v, 
                               const std::vector<measurement>& expected) {

                if (v.size() != expected.size()) 
                    throw std::invalid_argument("Can't operate a chi square funtion on vectors of different size"); 
                
                measurement accu = measurement(0., v[0].units()); 
                for (size_t i{}; i < v.size(); ++i) 
                    accu += square((v[i] - expected[i])) / expected[i]; 
                
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
                            sum_XX += square(xData[i]); 
                            sum_Y += yData[i];
                            sum_XY += xData[i] * yData[i]; 

                        }

                        delta = (N * sum_XX - square(sum_X));
                        slope_ = (N * sum_XY - sum_X * sum_Y) / delta;
                        intercept_ = (sum_XX * sum_Y - sum_X * sum_XY) / delta;

                        for (size_t i{}; i < N; ++i) 
                            sigma_y += square(yData[i] - intercept_.as_measurement() - slope_.as_measurement() * xData[i]); 

                        sigma_y /= (N - 2);

                        slope_.uncertainty() = sqrt(N * sigma_y / delta).value();
                        intercept_.uncertainty() = sqrt(sigma_y * sum_XX / delta).value();

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

                            weight = square(yData[i].uncertainty_as_measurement().inv());
                            wsum += weight; 
                            wsum_X += xData[i].as_measurement() * weight;
                            wsum_Y += yData[i].as_measurement() * weight;
                            wsum_XX += square(xData[i].as_measurement()) * weight; 
                            wsum_XY += xData[i].as_measurement() * yData[i].as_measurement() * weight; 

                        }

                        delta = (wsum * wsum_XX - square(wsum_X));
                        slope_ = uncertain_measurement((wsum * wsum_XY - wsum_X * wsum_Y) / delta, sqrt(wsum / delta).value());
                        intercept_ = uncertain_measurement((wsum_XX * wsum_Y - wsum_X * wsum_XY) / delta, sqrt(wsum_XX / delta).value());

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

                            weight = square((square((xData[i].uncertainty_as_measurement() * sigma_y_from_x)) + square(yData[i].uncertainty_as_measurement())).inv());
                            wsum += weight; 
                            wsum_X += xData[i].as_measurement() * weight;
                            wsum_Y += yData[i].as_measurement() * weight;
                            wsum_XX += square(xData[i].as_measurement()) * weight; 
                            wsum_XY += xData[i].as_measurement() * yData[i].as_measurement() * weight; 

                        }

                        delta = (wsum * wsum_XX - square(wsum_X));
                        slope_ = uncertain_measurement((wsum * wsum_XY - wsum_X * wsum_Y) / delta, sqrt(wsum / delta).value());
                        intercept_ = uncertain_measurement((wsum_XX * wsum_Y - wsum_X * wsum_XY) / delta, sqrt(wsum_XX / delta).value());

                    }

 
                    constexpr uncertain_measurement predict(const measurement x) const noexcept {

                        return intercept_ + slope_ * x;

                    }


                    constexpr uncertain_measurement predict(const uncertain_measurement x) const noexcept {

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

                    constexpr double twom12{0.000244140625};
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


                constexpr double gauss(const double& mean, 
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
