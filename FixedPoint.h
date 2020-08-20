/*
 * WiseMansFixedPoint header only C++ fixed point module. This fixed point
 * implementation supports numbers of varying integer and fractional lengths
 * and it employs the most basic arithmetic functions +,-,*,/ with proper
 * rounding when desiered. It is created with a hardware implementation point of
 * view and strives to behave as close to the VHDL IEEE signed and unsigned data
 * types.
 *
 * For the penalty of some greater run-time, the user can enable over-/underflow
 * checks by compiling the header with preprocessor macro
 * '_DEBUG_SHOW_OVERFLOW_INFO' defined (commandline option
 * '-D_DEBUG_SHOW_OVERFLOW_INFO' for GCC or CLANG) which will display some
 * over-/underflow info during execution.
 *
 * CAVIATS:
 *
 *
 * Author: Mikael Henriksson [www.github.com/miklhh]
 */

#ifndef _WISE_MANS_FIXED_POINT_H
#define _WISE_MANS_FIXED_POINT_H

#include "ttmath/ttmath.h"

#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <cstdint>


/*
 * Wide integer types. The 128 bit integers are used as underlying data type
 * for the fixed point numbers and the 256 bit integers are used for multiplying
 * and dividing the fixed point numbers.
 */
using int128_t = ttmath::Int<2>;
using int256_t = ttmath::Int<4>;
using uint128_t = ttmath::UInt<2>;
using uint256_t = ttmath::UInt<4>;
using float64_t = ttmath::Big<11,52>; // TTMath IEEE 754 double floating point.


/*
 * Compile time template structures for retrieving the extended integer size of
 * the underlying signed and unsigned integer types.
 */
template<typename T> struct extend_int {};
template<> struct extend_int<int128_t> { using type = int256_t; };
template<> struct extend_int<uint128_t> { using type = uint256_t; };


/*
 * Simple TTMath to_string function that generates a string of the underlying 
 * 2's complement data as a hex string. Not padded in any way.
 */
template <unsigned long _N> 
static inline std::string to_string_hex(const ttmath::Int<_N> &a) 
{
    const char hex_alphabet[] = "0123456789ABCDEF";

    ttmath::UInt<_N> n{ a };
    std::stringstream ss{};
    std::string result{};
    while (n != 0)
    {   
        ss << hex_alphabet[(n%0x10).ToUInt()];
        n >>= 4;
    }   
    result = ss.str();
    std::reverse(result.begin(), result.end());
    return result.empty() ? std::string("0") : result;
}


/*
 * Fixed point base type for common operations between signed and unsigned fixed
 * point types.
 */
template <int INT_BITS, int FRAC_BITS, typename _128_INT_TYPE>
class BaseFixedPoint
{
public:

    /*
     * The length of the integer part of the fixed point number should be less
     * than or equal to 64 bits due to the underlying 128 bit data type. For the
     * fractional part the same condition holds, but one extra bit is required 
     * to guaranteeing correct rounding when needed.
     */
    static_assert(INT_BITS <= 64,
            "Integer bits need to be less than or equal to 64 bits.");
    static_assert(FRAC_BITS < 64,
            "Fractional bits need to be strictly less than 64 bits.");
    static_assert(INT_BITS + FRAC_BITS > 0,
            "Need at least one bit of representation.");

    /*
     * Friend declaration of addition and subtraction arithmetic operators on fixed 
     * point numbers.
     */
    template<
        int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
        int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
    LHS<std::max(LHS_INT_BITS,RHS_INT_BITS)+1,std::max(LHS_FRAC_BITS,RHS_FRAC_BITS)>
    friend operator+(const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs, 
              const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs);

    template<
        int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
        int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
    LHS<std::max(LHS_INT_BITS,RHS_INT_BITS)+1,std::max(LHS_FRAC_BITS,RHS_FRAC_BITS)>
    friend operator-(const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs, 
              const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs);

    template<int _INT_BITS, int _FRAC_BITS, template<int,int> class RHS>
    friend RHS<_INT_BITS,_FRAC_BITS> operator-(const RHS<_INT_BITS,_FRAC_BITS> &rhs);


    /*
     * Friend declaration of multiplication and division arithmetic operators on 
     * fixed point numbers.
     */
    template<
        int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
        int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
    friend LHS<LHS_INT_BITS+RHS_INT_BITS,LHS_FRAC_BITS+RHS_FRAC_BITS>
    operator*(const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs, 
              const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs);

    template<
        int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
        int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
    friend LHS<LHS_INT_BITS,LHS_FRAC_BITS>
    operator/(const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs, 
              const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs);


    /*
     * Functions for getting the underlying data type.
     */
    _128_INT_TYPE get_num() const noexcept { return num; }

    virtual _128_INT_TYPE get_num_sign_extended() const noexcept = 0;


    /*
     * Friend declaration of rounding of fixed point numbers.
     */
    template<
        int LHS_INT_BITS, int LHS_FRAC_BITS,
        int RHS_INT_BITS, int RHS_FRAC_BITS, template<int,int> class RHS>
    friend RHS<LHS_INT_BITS,LHS_FRAC_BITS> rnd(const RHS<RHS_INT_BITS, RHS_FRAC_BITS> &rhs);


protected:

    /*
     * Helper function for retrieving a string of the fixed point fractional 
     * part on quotient form.
     */
    std::string get_frac_quotient() const noexcept
    {
        uint128_t num_ext = get_num_sign_extended();
        uint64_t numerator{ uint64_t(num_ext.ToUInt()) >> (64 - FRAC_BITS) };
        uint64_t denominator{ 1ull << FRAC_BITS };
        return std::to_string(numerator) + "/" + std::to_string(denominator);
    }


    /*
     * Friend declarations for accessing protected members between different
     * types of fixed point numbers, i.e, between template instances with
     * different wordlenths.
     */
    template <int _INT_BITS, int _FRAC_BITS, typename __128_INT_TYPE>
    friend class BaseFixedPoint;
    template <int _INT_BITS, int _FRAC_BITS>
    friend class SignedFixedPoint;
    template <int _INT_BITS, int _FRAC_BITS>
    friend class UnsignedFixedPoint;


    /*
     * Apply a bit mask to the internal representation to erase any bits outside
     * the intended number range.
     */
    void apply_bit_mask() noexcept
    {
        constexpr int WIDTH = INT_BITS+FRAC_BITS;
        this->num &= ((int128_t(1) << WIDTH)-1) << (64-FRAC_BITS);
    }

    void apply_bit_mask_int() noexcept
    {
        this->num &= (int128_t(1) << (64+INT_BITS)) - 1;
    }

    void apply_bit_mask_frac() noexcept
    {
        this->num &= ~((int128_t(1) << (64-FRAC_BITS)) - 1);
    }


    /*
     * Method for rounding the result of some operation to the closest fixed 
     * point number in the current representation.
     */
    void round() noexcept
    {
        num += _128_INT_TYPE(1) << (63-FRAC_BITS);
    }


    /*
     * Underlying data type. It will either be a signed or unsigned 128-bit
     * integer.
     */
    _128_INT_TYPE num{};
    using int_type = _128_INT_TYPE;


    /*
     * Protected constructor since base type should not be uninstantiatable.
     */
    BaseFixedPoint() = default;


    /*
     * Protected destructor.
     */
    virtual ~BaseFixedPoint() = default;
};


/*
 * Signed fixed point data type.
 */
template <int INT_BITS, int FRAC_BITS>
class SignedFixedPoint : public BaseFixedPoint<INT_BITS,FRAC_BITS,int128_t>
{
public:

    SignedFixedPoint() = default;


    /*
     * Conversion to and from floating point numbers.
     */
    explicit SignedFixedPoint(double a)
    {
        long n = lround(std::ceil(std::log2( std::abs(a) + 1 ) + 1));
        this->num = std::lround(a * double(1ul << (64-n)));
        this->num <<= n;
        this->round();
        this->apply_bit_mask_frac();
        this->num = this->get_num_sign_extended();
    }

    explicit operator double() const
    {
        return float64_t(this->get_num_sign_extended()).ToDouble() / 
               std::pow(2.0, 64);
    }


    /*
     * Specialized to_string function for signed fixed point numbers.
     */
    std::string to_string() const noexcept
    {
        // Extract integer part.
        int128_t sign_num = this->get_num_sign_extended();
        long int integer = static_cast<long int>((sign_num >> 64).ToUInt());
        std::string integer_str( std::to_string(integer) );

        // Append fractional part.
        return integer_str + " + " + this->get_frac_quotient();
    }


    /*
     * Display the state of the fixed point number through the retuned string.
     * The string contains formated output for debuging purposes.
     */
    std::string get_state() const noexcept
    {
        std::stringstream ss{};
        std::string internal = to_string_hex(this->num);
        internal = std::string(32-internal.length(), '0') + internal;
        internal.insert(16, 1, '|');
        return internal;
    }


    /*
     * Assignment operator and assignment constructor for fixed point numbers.
     */
    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_128_INT_TYPE>
    SignedFixedPoint<INT_BITS, FRAC_BITS> &
       operator=(const BaseFixedPoint<RHS_INT_BITS,
                                      RHS_FRAC_BITS, 
                                      RHS_128_INT_TYPE> &rhs) noexcept
    {
        this->num = rhs.get_num_sign_extended();
        this->apply_bit_mask_frac();
        return *this;
    }

    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_128_INT_TYPE>
    SignedFixedPoint(const BaseFixedPoint<RHS_INT_BITS,
                                          RHS_FRAC_BITS, 
                                          RHS_128_INT_TYPE> &rhs) noexcept
    {
        this->num = rhs.get_num_sign_extended();
        this->apply_bit_mask_frac();
    }


    /*
     * Assignment from other fixed point number with proper rounind.
     */
    template <int RHS_INT_BITS,int RHS_FRAC_BITS, typename RHS_128_INT_TYPE>
    SignedFixedPoint<INT_BITS, FRAC_BITS> 
        &rnd(const BaseFixedPoint<RHS_INT_BITS, 
                                  RHS_FRAC_BITS, 
                                  RHS_128_INT_TYPE> &rhs)
    {
        this->num = rhs.num;
        this->round();
        this->apply_bit_mask_frac();
        this->num = this->get_num_sign_extended();
        return *this;
    }

    /*
     * Friend declaration for different template instansiations of one self.
     */
    template <int _INT_BITS, int _FRAC_BITS>
    friend class SignedFixedPoint;


    /*
     * Friend functions for rounding and saturation.
     */
    template<
        int LHS_INT_BITS, int LHS_FRAC_BITS,
        int RHS_INT_BITS, int RHS_FRAC_BITS, template<int,int> class RHS>
    friend RHS<LHS_INT_BITS,LHS_FRAC_BITS> rnd(const RHS<RHS_INT_BITS, RHS_FRAC_BITS> &rhs);

    template <int LHS_INT,int LHS_FRAC,int RHS_INT,int RHS_FRAC>
    friend SignedFixedPoint<LHS_INT, LHS_FRAC> sat(
                const SignedFixedPoint<RHS_INT, RHS_FRAC> &rhs);

    /*
     * Friend arithmetic.
     */
    template<
        int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
        int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
    friend LHS<LHS_INT_BITS,LHS_FRAC_BITS>
    operator/(const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs, 
              const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs);

    template<int _INT_BITS, int _FRAC_BITS, template<int,int> class RHS>
    friend RHS<_INT_BITS,_FRAC_BITS> operator-(const RHS<_INT_BITS,_FRAC_BITS> &rhs);


    /*
     * Returns the internal num representation sign extended, that is, num with
     * all bits more significant than the sign bit set to the value of the sign
     * bit.
     */
    int128_t get_num_sign_extended() const noexcept override
    {
        if ( sign() )
            return this->num | ~((int128_t(1) << (64+INT_BITS)) - 1);
        else
            return this->num &  ((int128_t(1) << (64+INT_BITS)) - 1);
    }


private:
    /*
     * Get the sign of the number.
     */
    bool sign() const noexcept
    {
        return int128_t(0) != ( this->num & (int128_t(1) << (64+INT_BITS-1)) );
    }
};


/*
 * Unsigned fixed point data type.
 */
template <int INT_BITS, int FRAC_BITS>
class UnsignedFixedPoint : public BaseFixedPoint<INT_BITS,FRAC_BITS,uint128_t>
{
    /*
     * get_num_sign_extended() should always return a fully masked num in this
     * class.
     */
};


/*
 * Addition operator for fixed point numbers.
 */
template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
LHS<std::max(LHS_INT_BITS,RHS_INT_BITS)+1,std::max(LHS_FRAC_BITS,RHS_FRAC_BITS)>
operator+(const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs, 
          const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs)
{
    // No sign extension or masking needed due to correct word length of result.
    constexpr int RES_INT_BITS = std::max(LHS_INT_BITS,RHS_INT_BITS)+1;
    constexpr int RES_FRAC_BITS = std::max(LHS_FRAC_BITS,RHS_FRAC_BITS);
    LHS<RES_INT_BITS,RES_FRAC_BITS> res{};
    res.num = lhs.num + rhs.num;
    return res;
}

template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
LHS<LHS_INT_BITS,LHS_FRAC_BITS> &
operator+=(LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs, 
           const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs)
{
    // Sign extension and masking is performed in assigment operator.
    return lhs = lhs + rhs;
}


/*
 * Subtraction operator for fixed point numbers.
 */
template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
LHS<std::max(LHS_INT_BITS,RHS_INT_BITS)+1,std::max(LHS_FRAC_BITS,RHS_FRAC_BITS)>
operator-(const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs, 
          const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs)
{
    // No sign extension or masking needed due to correct word length of result.
    constexpr int RES_INT_BITS = std::max(LHS_INT_BITS,RHS_INT_BITS)+1;
    constexpr int RES_FRAC_BITS = std::max(LHS_FRAC_BITS,RHS_FRAC_BITS);
    LHS<RES_INT_BITS,RES_FRAC_BITS> res{};
    res.num = lhs.num - rhs.num;
    return res;
}

template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
LHS<LHS_INT_BITS,LHS_FRAC_BITS> &
operator-=(LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs, 
           const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs)
{
    // Sign extension and masking is performed in assigment operator.
    return lhs = lhs - rhs;
}


/*
 * Multiplication operator for fixed point numbers.
 */
template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
LHS<LHS_INT_BITS+RHS_INT_BITS,LHS_FRAC_BITS+RHS_FRAC_BITS>
operator*(const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs, 
          const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs)
{
    // No sign extension or masking needed due to correct word length of result.
    LHS<LHS_INT_BITS+RHS_INT_BITS, LHS_FRAC_BITS+RHS_FRAC_BITS> res{};
    typename extend_int<typename LHS<1,0>::int_type>::type long_lhs = lhs.num;
    typename extend_int<typename LHS<1,0>::int_type>::type long_rhs = rhs.num;
    res.num = (long_lhs * long_rhs) >> 64;
    return res;
}

template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
LHS<LHS_INT_BITS,LHS_FRAC_BITS> &
operator*=(LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs, 
           const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs)
{
    // Sign extension and masking is performed in assigment operator.
    return lhs = lhs * rhs;
}


/*
 * Division operator for fixed point numbers.
 */
template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
LHS<LHS_INT_BITS,LHS_FRAC_BITS>
operator/(const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs, 
          const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs)
{
    // Sign extension and masking needed due to uncorrect result word length.
    using long_int = typename extend_int<typename LHS<1,0>::int_type>::type;
    LHS<LHS_INT_BITS, LHS_FRAC_BITS> res{};
    long_int long_lhs = lhs.num;
    long_int long_rhs = rhs.num;
    res.num = (long_int(long_lhs << 128) / long_rhs) >> 64;
    res.num = res.get_num_sign_extended();
    res.apply_bit_mask_frac();
    return res;
}

template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
LHS<LHS_INT_BITS,LHS_FRAC_BITS> &
operator/=(LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs, 
           const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs)
{
    // Sign extension and masking is performed in assigment operator.
    return lhs = lhs / rhs;
}


/*
 * Comparison operators for fixed point numbers.
 */
template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
bool operator==(
        const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs,
        const BaseFixedPoint<RHS_INT_BITS, RHS_FRAC_BITS, RHS_INT_TYPE> &rhs)
{
    return lhs.get_num() == rhs.get_num();
}

template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
bool operator!=(
        const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs,
        const BaseFixedPoint<RHS_INT_BITS, RHS_FRAC_BITS, RHS_INT_TYPE> &rhs)
{
    return lhs.get_num() != rhs.get_num();
}

template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
bool operator<(
        const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs,
        const BaseFixedPoint<RHS_INT_BITS, RHS_FRAC_BITS, RHS_INT_TYPE> &rhs)
{
    return lhs.get_num() < rhs.get_num();
}

template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
bool operator<=(
        const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs,
        const BaseFixedPoint<RHS_INT_BITS, RHS_FRAC_BITS, RHS_INT_TYPE> &rhs)
{
    return lhs.get_num() <= rhs.get_num();
}

template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
bool operator>(
        const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs,
        const BaseFixedPoint<RHS_INT_BITS, RHS_FRAC_BITS, RHS_INT_TYPE> &rhs)
{
    return lhs.get_num() > rhs.get_num();
}

template<
    int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
bool operator>=(
        const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs,
        const BaseFixedPoint<RHS_INT_BITS, RHS_FRAC_BITS, RHS_INT_TYPE> &rhs)
{
    return lhs.get_num() >= rhs.get_num();
}


/*
 * Unary negation of fixed point numbers.
 */
template<int INT_BITS, int FRAC_BITS, template<int,int> class RHS>
RHS<INT_BITS,FRAC_BITS> operator-(const RHS<INT_BITS,FRAC_BITS> &rhs)
{
    RHS<INT_BITS,FRAC_BITS> res{};
    res.num = -rhs.num;
    res.num = res.get_num_sign_extended();
    res.apply_bit_mask_frac();
    return res;
}


/*
 * Rounding for fixed point numbers.
 */
template<
    int LHS_INT_BITS, int LHS_FRAC_BITS,
    int RHS_INT_BITS, int RHS_FRAC_BITS, template<int,int> class RHS>
RHS<LHS_INT_BITS,LHS_FRAC_BITS> rnd(const RHS<RHS_INT_BITS, RHS_FRAC_BITS> &rhs)
{
    RHS<LHS_INT_BITS, LHS_FRAC_BITS> res{};
    res.num = rhs.num;
    res.round();
    res.apply_bit_mask_frac();
    res.num = res.get_num_sign_extended();
    return res;
}


/*
 * Saturation for signed fixed point numbers.
 */
template <int LHS_INT_BITS,int LHS_FRAC_BITS,int RHS_INT_BITS,int RHS_FRAC_BITS>
SignedFixedPoint<LHS_INT_BITS, LHS_FRAC_BITS> sat(
        const SignedFixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs)
{
    SignedFixedPoint<LHS_INT_BITS, LHS_FRAC_BITS> res{};
    int128_t num = rhs.get_num_sign_extended();
    int128_t overflow_mask = ~((int128_t(1) << (64+LHS_INT_BITS-1))-1);
    if ( (num ^ overflow_mask) != 0 )
    {
        if ( num < 0 )
        {
            // Min value.
            res.num = ~((int128_t(1) << (64+LHS_INT_BITS-1)) - 1);
        }
        else
        {
            // Max value.
            res.num = (int128_t(1) << (64+LHS_INT_BITS-1)) - 1;
        }
        res.apply_bit_mask_frac();
        return res;
    }
    else
    {
        res.num = num;
        res.apply_bit_mask_frac();
        return res;
    }
}


/*
 * Print-out to C++ stream object on the form '<int> + <frac>/<2^<frac_bits>'.
 * Good for debuging'n'stuff.
 */
template <int INT_BITS, int FRAC_BITS>
std::ostream &operator<<(
        std::ostream &os, const SignedFixedPoint<INT_BITS, FRAC_BITS> &rhs)
{
    return os << rhs.to_string();
}


#endif // Header guard ending.

