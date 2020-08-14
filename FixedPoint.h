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
#include <iostream> // DEBUGING!
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
 * Fixed point base type for common upperation between signed and unsigned
 * types.
 */
template <int INT_BITS, int FRAC_BITS, typename _128_INT_TYPE>
class BaseFixedPoint
{
public:
    /*
     * The length of the integer part of the fixed point number should be less
     * than or equal to 64 bits due to the underlying 128 bit data type. For the
     * fractional part the same thingn holds, but one extra bit is required for
     * guaranteeing correct rounding.
     */
    static_assert(INT_BITS <= 64,
            "Integer bits need to be less than or equal to 64 bits.");
    static_assert(FRAC_BITS < 64,
            "Fractional bits need to be strictly less than 64 bits.");
    static_assert(INT_BITS + FRAC_BITS > 0,
            "Need at least one bit of representation.");

protected:
    /*
     * Helper function for retrieving a string of the fractional part on
     * quotient form.
     */
    std::string get_frac_quotient() const noexcept
    {
        uint64_t numerator{ uint64_t(num.ToInt()) >> (64 - FRAC_BITS) };
        uint64_t denominator{ 1ull << FRAC_BITS };
        return std::to_string(numerator) + "/" + std::to_string(denominator);
    }

    /*
     * Friend declaration for accessing protected members between different 
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
     * Apply the bit mask to the internal representation to erase any bits
     * outside the intended number range.
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
     * This method will round the result of some operation to the closesd fixed
     * point number in the current representation.
     */
    void round() noexcept
    {
        num += _128_INT_TYPE(1) << (63-FRAC_BITS);
    }
    void round_and_mask() noexcept
    {
        round();
        apply_bit_mask();
    }

    /*
     * The internal representation of the fixed point number is considered dirty
     * if any bit outside of the desiered range <a,b> is set.
     */
    bool dirty() const noexcept
    {
        constexpr int WIDTH = INT_BITS+FRAC_BITS;
        return this->num & ~(((int128_t(1) << WIDTH)-1) << (64-FRAC_BITS));
    }
    bool dirty_int() const noexcept
    {
        return this->num & ~((int128_t(1) << (64+INT_BITS)) - 1);
    }
    bool dirty_frac() const noexcept
    {
        return this->num & ~((int128_t(1) << (64-INT_BITS)) - 1);
    }

    /*
     * Returns the internal number representation 'num', sign extended. For the
     * unsigned derived class, this method just return num.
     */
    virtual _128_INT_TYPE get_num_sign_extended() const noexcept;

    /*
     * Underlying data type. It will either be a signed or unsigned 128-bit
     * integer.
     */
    _128_INT_TYPE num{};

    /*
     * Protected constructor since base type should not be instantiatable.
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
        // Extract integer part. Result of the right shift operation is the 
        // sign extended number as an unsigned type.
        int128_t sign_num = this->get_num_sign_extended();
        long int integer = static_cast<long int>((sign_num >> 64).ToUInt());
        std::string integer_str( std::to_string(integer) );

        // Append fractional part.
        return integer_str + " + " + this->get_frac_quotient();
    }


    /*
     * Addition and subtraction arithmetic.
     */
    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename _RHS_128_INT_TYPE>
    SignedFixedPoint<INT_BITS, FRAC_BITS>
        operator+(const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,_RHS_128_INT_TYPE> &rhs) const noexcept
    {
        SignedFixedPoint<INT_BITS, FRAC_BITS> res{};
        res.num = this->num + rhs.num;
        res.apply_bit_mask_frac();
        return res;
    }

    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename _RHS_128_INT_TYPE>
    SignedFixedPoint<INT_BITS, FRAC_BITS>
        operator-(const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,_RHS_128_INT_TYPE> &rhs) const noexcept
    {
        SignedFixedPoint<INT_BITS, FRAC_BITS> res{};
        res.num = this->num - rhs.num;
        res.apply_bit_mask_frac();
        return res;
    }

    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename _RHS_128_INT_TYPE>
    SignedFixedPoint<INT_BITS, FRAC_BITS>
        &operator+=(const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,_RHS_128_INT_TYPE> &rhs) noexcept
    {
        this->num += rhs.num;
        this->apply_bit_mask_frac();
        return *this;
    }

    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename _RHS_128_INT_TYPE>
    SignedFixedPoint<INT_BITS, FRAC_BITS>
        &operator-=(const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,_RHS_128_INT_TYPE> &rhs) noexcept
    {
        this->num -= rhs.num;
        this->apply_bit_mask_frac();
        return *this;
    }


    /*
     * Multiplication and division airthmetic. Note that, just like in VHDL and
     * Verilog, the result of <a,b>*<c,d> = <a+c,b+d>. Using operator*= will 
     * of course result in a fixed point number of word length equal to the left
     * hand side. Like wise, to conform with the VHDL and Verilog way of doing
     * division, the result of <a,b>/<c,d> = <a,b>.
     */
    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename _RHS_128_INT_TYPE>
    SignedFixedPoint<INT_BITS+RHS_INT_BITS, FRAC_BITS+RHS_FRAC_BITS>
        operator*(const BaseFixedPoint<RHS_INT_BITS,
                                       RHS_FRAC_BITS,
                                       _RHS_128_INT_TYPE> &rhs) const noexcept
    {
        SignedFixedPoint<INT_BITS+RHS_INT_BITS, FRAC_BITS+RHS_FRAC_BITS> res{};
        int256_t long_lhs = this->get_num_sign_extended();
        int256_t long_rhs = rhs.get_num_sign_extended();
        int256_t long_res = long_lhs * long_rhs;
        long_res >>= 64;
        res.num = long_res;
        res.apply_bit_mask_frac();
        return res;
    }

    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename _RHS_128_INT_TYPE>
    SignedFixedPoint<INT_BITS, FRAC_BITS>
        &operator*=(const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,_RHS_128_INT_TYPE> &rhs) noexcept
    {
        SignedFixedPoint<INT_BITS+RHS_INT_BITS, FRAC_BITS+RHS_FRAC_BITS> res{};
        res = *this * rhs;
        return *this = res;
    }

    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename _RHS_128_INT_TYPE>
    SignedFixedPoint<INT_BITS, FRAC_BITS>
        operator/(const BaseFixedPoint<RHS_INT_BITS,
                                       RHS_FRAC_BITS,
                                       _RHS_128_INT_TYPE> &rhs) const
    {
        SignedFixedPoint<INT_BITS+RHS_INT_BITS, FRAC_BITS+RHS_FRAC_BITS> res{};
        int256_t long_lhs = this->get_num_sign_extended();
        int256_t long_rhs = rhs.get_num_sign_extended();
        long_lhs <<= 128;
        int256_t long_res = long_lhs / long_rhs;
        long_res >>= 64;
        res.num = long_res;
        res.apply_bit_mask_frac();
        return res;
    }

    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename _RHS_128_INT_TYPE>
    SignedFixedPoint<INT_BITS, FRAC_BITS>
        &operator/=(const BaseFixedPoint<RHS_INT_BITS,
                                         RHS_FRAC_BITS,
                                         _RHS_128_INT_TYPE> &rhs)
    {
        SignedFixedPoint<INT_BITS, FRAC_BITS> res{};
        res = *this / rhs;
        return *this = res;
    }


    /*
     * Unary negation.
     */
    SignedFixedPoint<INT_BITS, FRAC_BITS> operator-() const noexcept
    {
        SignedFixedPoint<INT_BITS, FRAC_BITS> res{};
        res.num = -( this->get_num_sign_extended() );
        res.apply_bit_mask_frac();
        return res;
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
     * Comparison operators.
     */
    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_128_INT_TYPE>
    bool operator==(const BaseFixedPoint<RHS_INT_BITS, 
                                         RHS_FRAC_BITS, 
                                         RHS_128_INT_TYPE> &rhs) const noexcept
    {
        return this->get_num_sign_extended() == rhs.get_num_sign_extended();
    }


private:
    /*
     * Get the sign of the number.
     */
    bool sign() const noexcept
    {
        return int128_t(0) != ( this->num & (int128_t(1) << (64+INT_BITS-1)) );
    }

    /*
     * Returns the internal num representation sign extended, that is, num with
     * all bits more significant than the sign bit set to the value of the sign
     * bit.
     */
    int128_t get_num_sign_extended() const noexcept
    {
        if ( sign() )
            return this->num | ~((int128_t(1) << (64+INT_BITS)) - 1);
        else
            return this->num &  ((int128_t(1) << (64+INT_BITS)) - 1);
    }
};

template <int INT_BITS, int FRAC_BITS>
class UnsignedFixedPoint : public BaseFixedPoint<INT_BITS,FRAC_BITS,uint128_t>
{

};

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


/*
 * ----------------------------------------------------------------------------
 *  Underneath this line lies code from the old PoorMansFixedPoint for ref.
 * ----------------------------------------------------------------------------
 */



///*
// * Debuging stuff for being able to display over-/underflows.
// */
//#ifdef _DEBUG_SHOW_OVERFLOW_INFO
//    #include <sstream>
//
//    /*
//     * If using this debuging functionality, define in this function how an
//     * output string should be passed to the user.
//     */
//    #include <iostream>
//    void _DEBUG_PRINT_FUNC(const char *str)
//    {
//        // Example debug feed forward.
//        std::cerr << str << std::endl;
//    }
//#endif
//
//
///*
// * Type FixedPoint begin.
// */
//template <int INT_BITS, int FRAC_BITS>
//class FixedPoint
//{
//    /*
//     * Constraints put on the FixedPoint numbers by this implementaiton.
//     */
//    static_assert(INT_BITS <= 32,
//            "Integer bits need to be less than or equal to 32 bits.");
//    static_assert(FRAC_BITS <= 32,
//            "Fractional bits need to be less than or equal to 32 bits.");
//    static_assert(INT_BITS + FRAC_BITS > 0,
//            "Need at least one bit of representation.");
//
//    /*
//     * We rely on right shift of signed long long integers to be equivilant with
//     * arithmetic right shifts. When writting this comment, this is proposed to
//     * be the standard behaviour of right shifts in the C++20 standard draft,
//     * but for current standards (<= C++17) right shifting of signed integers is
//     * implementation defined.
//     */
//    static_assert( (-2ll >> 1) == -1ll,
//            "We rely on signed right shifts to be arithmetic." );
//
//protected:
//    /*
//     * Long long is guaranteed to be atleast 64-bits wide. We use the 32 most
//     * significant bits to store the integer part and the 32 least significant
//     * bits to store the fraction.
//     */
//    long long num{};
//
//    /*
//     * Friend declaration for accessing 'num' between different types, i.e,
//     * between template instances with different wordlenth.
//     */
//    template <int _INT_BITS, int _FRAC_BITS>
//    friend class FixedPoint;
//
//    /*
//     * Private rounding method. This method will round the result of some
//     * operation to the closed fixed point number in the current representation.
//     * It also contains support for displaying over-/underflows, see code.
//     */
//    void round() noexcept
//    {
//        /*
//         * Perform rounding by adding 2^(-FRAC_BITS)/2.
//         */
//        if (FRAC_BITS < 32)
//            this->num += 1ll << (31-FRAC_BITS);
//
//    #ifdef _DEBUG_SHOW_OVERFLOW_INFO
//        /*
//         * If debug overflow info mode is enebaled, test for over-/underflow
//         * in the result and present user with a warning.
//         */
//        if ( test_over_or_underflow() )
//        {
//            // Print operation result.
//            std::stringstream ss{};
//            ss << ( ((1ull << 63) & this->num) ? "Underflow " : "Overflow " );
//            ss << "in node <" << INT_BITS << "," << FRAC_BITS << "> ";
//            ss << "'of value: " << (this->num >> 32) << " + ";
//            ss << this->get_frac_quotient() << ", ";
//
//            // Apply bitmask and print truncated result.
//            if (INT_BITS+FRAC_BITS < 64)
//                this->num &= ((1ull<<(INT_BITS+FRAC_BITS))-1) << (32-FRAC_BITS);
//            ss << "truncated to: " << (this->get_num_sign_extended() >> 32);
//            ss << " + " << this->get_frac_quotient();
//            _DEBUG_PRINT_FUNC(ss.str().c_str());
//        }
//        else
//        {
//            if (INT_BITS+FRAC_BITS < 64)
//                this->num &= ((1ull<<(INT_BITS+FRAC_BITS))-1) << (32-FRAC_BITS);
//        }
//    #else
//        /*
//         * Debugmode disabled. Just apply the bitmask.
//         */
//        if (INT_BITS+FRAC_BITS < 64)
//            this->num &= ((1ull<<(INT_BITS+FRAC_BITS))-1) << (32-FRAC_BITS);
//    #endif
//    }
//
//    /*
//     * Get the current number sign extended to Q(32, 32) format.
//     */
//    long long get_num_sign_extended() const noexcept
//    {
//        /*
//         * Instead of of testing whether we need to sign extend the internal
//         * number num, we left shift it (unsigned, logically) all the way to the
//         * MSb, and then shift it (signed, arithmeticaly) back to its original
//         * position. This seem to have performance benifits.
//         */
//        using uns_ll = unsigned long long;
//        uns_ll l = static_cast<uns_ll>(num) << (32-INT_BITS);
//        return static_cast<long long>(l) >> (32-INT_BITS);
//    }
//
//    /*
//     * Private method for testing over/underflow.
//     */
//    bool test_over_or_underflow() const noexcept
//    {
//        unsigned long long msb_extended = (this->num >> (31+INT_BITS));
//        return !( (msb_extended == -1ull) || (msb_extended == 0ull) );
//    }
//
//public:
//    FixedPoint() = default;
//    virtual ~FixedPoint() = default;
//
//    /*
//     * Constructor for initialization from other fixed point number. Note that
//     * if the number cannot fit into the FixedPoint type, it will be truncated.
//     */
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    FixedPoint(const FixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs) noexcept
//    {
//        this->num = rhs.get_num_sign_extended();
//        this->round();
//    }
//    FixedPoint(const FixedPoint<INT_BITS, FRAC_BITS> &rhs) noexcept
//    {
//        /*
//         * Initialization from FixedPoint numbers with same length dont need
//         * rounding.
//         */
//        this->num = rhs.num;
//    }
//
//    /*
//     * Constructor for floating point number inputs.
//     */
//    explicit FixedPoint(double a)
//    {
//        this->num = std::llround(a * static_cast<double>(1ll << 32));
//        this->round();
//    }
//
//    /*
//     * Constructor from integers.
//     */
//    explicit FixedPoint(int n) noexcept
//    {
//        this->num = static_cast<long long>(n) << 32;
//        this->round();
//    }
//
//    /*
//     * Constructor for setting the bit pattern of a FixedPoint number.
//     */
//    explicit FixedPoint(int i, unsigned f) noexcept
//    {
//        this->num = (static_cast<long long>(i) << 32) | (0xFFFFFFFFll & num);
//        this->num &= 0xFFFFFFFF00000000ll;
//        this->num |= 0xFFFFFFFFll & (static_cast<long long>(f)<<(32-FRAC_BITS));
//        this->round();
//    }
//
//    /*
//     * Get template arguments from FixedPoint.
//     */
//    constexpr int get_int_bits() const noexcept { return INT_BITS; }
//    constexpr int get_frac_bits() const noexcept { return FRAC_BITS; }
//
//    /*
//     * Retrieve a string of the fractional part of the FixedPoint number, useful
//     * for displaying the FixedPoint number content. The string will be on the
//     * format <frac>/2^{FRAC_BITS}, e.g, "5/32" or "13/2048".
//     */
//    std::string get_frac_quotient() const noexcept
//    {
//        using std::string; using std::to_string;
//        string numerator = to_string((this->num & 0xFFFFFFFF)>>(32-FRAC_BITS));
//        string denominator = to_string(1ll << FRAC_BITS);
//        return numerator + "/" + denominator;
//    }
//
//    /*
//     * To string function, good for FixedPoint printouts.
//     */
//    std::string to_string() const noexcept
//    {
//        long long num{ this->get_num_sign_extended() };
//        std::string res{ std::to_string(num >> 32) };
//        return res += std::string(" + ") += this->get_frac_quotient();
//    }
//
//    /*
//     * Assigment operators of FixedPoint numbers.
//     */
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    FixedPoint<INT_BITS, FRAC_BITS> &
//        operator=(const FixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs) noexcept
//    {
//        this->num = rhs.get_num_sign_extended();
//        this->round();
//        return *this;
//    }
//    FixedPoint<INT_BITS, FRAC_BITS> &
//        operator=(const FixedPoint<INT_BITS, FRAC_BITS> &rhs) noexcept
//    {
//        this->num = rhs.num;
//        return *this;
//    }
//
//    /*
//     * (explicit) Conversion to floating point number.
//     */
//    explicit operator double() const noexcept
//    {
//        // Test if sign extension is needed.
//        return static_cast<double>(this->get_num_sign_extended()) /
//               static_cast<double>(1ll << 32);
//    }
//
//    /*
//     * Unary negation operator.
//     */
//    FixedPoint<INT_BITS, FRAC_BITS> operator-() const noexcept
//    {
//        FixedPoint<INT_BITS, FRAC_BITS> res{};
//        res.num = -( this->get_num_sign_extended() );
//        res.round();
//        return res;
//    }
//
//    /*
//     * Addition/subtraction of FixedPoint numbers. Result will have word length
//     * equal to that of the left hand side operand.
//     */
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    FixedPoint<INT_BITS, FRAC_BITS>
//        operator+(const FixedPoint<RHS_INT_BITS,RHS_FRAC_BITS> &rhs)
//        const noexcept
//    {
//        FixedPoint<INT_BITS, FRAC_BITS> res;
//        res.num = this->get_num_sign_extended() + rhs.get_num_sign_extended();
//        res.round();
//        return res;
//    }
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    FixedPoint<INT_BITS, FRAC_BITS> &
//        operator+=(const FixedPoint<RHS_INT_BITS,RHS_FRAC_BITS> &rhs) noexcept
//    {
//        this->num = this->get_num_sign_extended() + rhs.get_num_sign_extended();
//        this->round();
//        return *this;
//    }
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    FixedPoint<INT_BITS, FRAC_BITS>
//        operator-(const FixedPoint<RHS_INT_BITS,RHS_FRAC_BITS> &rhs)
//        const noexcept
//    {
//        FixedPoint<INT_BITS, FRAC_BITS> res;
//        res.num = this->get_num_sign_extended() - rhs.get_num_sign_extended();
//        res.round();
//        return res;
//    }
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    FixedPoint<INT_BITS, FRAC_BITS> &
//        operator-=(const FixedPoint<RHS_INT_BITS,RHS_FRAC_BITS> &rhs) noexcept
//    {
//        this->num = this->get_num_sign_extended() - rhs.get_num_sign_extended();
//        this->round();
//        return *this;
//    }
//
//    /*
//     * Multilication of FixedPoint numbers. Result will have an integer and
//     * fractional wordlength equal to that of the sum of left hand side and
//     * right hand side operand integer and fractional wordlengths, but no
//     * longer than <32,32>.
//     */
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    FixedPoint<std::min(INT_BITS+RHS_INT_BITS, 32),
//               std::min(FRAC_BITS+RHS_FRAC_BITS, 32)>
//        operator*(const FixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs)
//        const noexcept
//    {
//        FixedPoint<std::min(INT_BITS+RHS_INT_BITS, 32),
//                   std::min(FRAC_BITS+RHS_FRAC_BITS, 32)> res{};
//
//        /*
//         * Scenario 1:
//         * The entire result of the multiplication can fit into one 64-bit
//         * integer. This code produces faster result when applicable.
//         */
//        if (INT_BITS+FRAC_BITS <= 32 && RHS_INT_BITS+RHS_FRAC_BITS <= 32)
//        {
//            long long op_a{ this->get_num_sign_extended() >> INT_BITS     };
//            long long op_b{   rhs.get_num_sign_extended() >> RHS_INT_BITS };
//            res.num = op_a * op_b;
//
//            // Shift result to the correct place.
//            if (INT_BITS + RHS_INT_BITS > 32)
//            {
//                res.num <<= INT_BITS + RHS_INT_BITS - 32;
//            }
//            else
//            {
//                res.num >>= 32 - INT_BITS - RHS_INT_BITS;
//            }
//            res.round();
//            return res;
//        }
//        /*
//         * Scenario 2:
//         * The entire result of the multiplication can fit into one 128-bit
//         * integer. Running this code takes a little longer time than running
//         * the code of scenario 1, probably due to the fact that this code won't
//         * be accelerated by any integer vectorization. However, this piece of
//         * code seems to work for all sizes of FixedPoints.
//         */
//        else
//        {
//            // Utilize the compiler extension of 128-bit wide integers to be
//            // able to store the exact result.
//            __extension__ __int128 op_a{ this->get_num_sign_extended() };
//            __extension__ __int128 op_b{   rhs.get_num_sign_extended() };
//            __extension__ __int128 res_128{ op_a * op_b };
//
//            // Shift result back to the form Q(32,32) and return result.
//            res.num = static_cast<long long>(res_128 >> 32);
//            res.round();
//            return res;
//        }
//    }
//    /*
//     * NOTE: Result of this operator will not change the wordlength of the the
//     * number.
//     */
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    FixedPoint<INT_BITS, FRAC_BITS> &
//        operator*=(const FixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs) noexcept
//    {
//        // Rounding is performed in operator*().
//        FixedPoint<INT_BITS, FRAC_BITS> res{ *this * rhs };
//        this->num = res.num;
//        return *this;
//    }
//
//    /*
//     * Division of FixedPoint numbers with integers. Result will have word
//     * length equal to that of the left hand side of the operator, but the
//     * precision of the result will not necessary represent such a wide number.
//     */
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    FixedPoint<INT_BITS, FRAC_BITS>
//        operator/(const FixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs)
//        const
//    {
//        // Note that Q(a,64) / Q(b,32) == Q(a-b,32).
//        FixedPoint<INT_BITS, FRAC_BITS> res{};
//        __extension__ __int128 dividend{ this->get_num_sign_extended() };
//        __extension__ __int128 divisor {   rhs.get_num_sign_extended() };
//        dividend <<= 32;
//
//        // Create and return result. Note that Q(a,64) / Q(b,32) == Q(c,32).
//        res.num = static_cast<long long>( dividend/divisor );
//        res.round();
//        return res;
//    }
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    FixedPoint<INT_BITS, FRAC_BITS> &
//        operator/=(const FixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs)
//    {
//        // Rounding is performed in operator/().
//        FixedPoint<INT_BITS, FRAC_BITS> res{ *this / rhs };
//        this->num = res.num;
//        return *this;
//    }
//
//    /*
//     * Comparison operators.
//     */
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    bool operator==(const FixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs)
//        const noexcept
//    {
//        return this->get_num_sign_extended() == rhs.get_num_sign_extended();
//    }
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    bool operator!=(const FixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs)
//        const noexcept
//    {
//        return this->get_num_sign_extended() != rhs.get_num_sign_extended();
//    }
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    bool operator<(const FixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs)
//        const noexcept
//    {
//        return this->get_num_sign_extended() < rhs.get_num_sign_extended();
//    }
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    bool operator<=(const FixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs)
//        const noexcept
//    {
//        return this->get_num_sign_extended() <= rhs.get_num_sign_extended();
//    }
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    bool operator>(const FixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs)
//        const noexcept
//    {
//        return this->get_num_sign_extended() > rhs.get_num_sign_extended();
//    }
//    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
//    bool operator>=(const FixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs)
//        const noexcept
//    {
//        return this->get_num_sign_extended() >= rhs.get_num_sign_extended();
//    }
//
//};
//
//
///*
// * Print-out to C++ stream object on the form '<int> + <frac>/<2^<frac_bits>'.
// * Good for debuging'n'stuff.
// */
//template <int INT_BITS, int FRAC_BITS>
//std::ostream &operator<<(
//        std::ostream &os, const FixedPoint<INT_BITS, FRAC_BITS> &rhs)
//{
//    return os << rhs.to_string();
//}

/*
 * Include guard end.
 */
