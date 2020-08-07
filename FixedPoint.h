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
 *   * For rounding to work properly, the fractional part of the FixedPoint
 *     number has to be strictly less than 32 bits. For numbers with 32
 *     fractional bits, the result will round towards -INF.
 *
 *   * The implementation requires support for the compiler extension '__int128'
 *     which is used for longer fixed point multiplications and divisions. As
 *     far as I know, this feature is not supported by the Visual C++ compiler.
 *
 *   * The implementation makes use of some implementation defined behaviour, in
 *     that is depends on right shifts of signed integers to realize arithmetic
 *     shifts. This is the default behaviour of many modern compilers, and if
 *     your compiler does not support this functionality, the code will generate
 *     a compile time error through a static assertion.
 *
 * Author: Mikael Henriksson [www.github.com/miklhh]
 */

#ifndef _WISE_MANS_FIXED_POINT_H
#define _WISE_MANS_FIXED_POINT_H

#include <ostream>
#include <string>
#include <cmath>
#include <iostream> // DEBUGING!
#include <algorithm>
#include <type_traits>
#include <cstdint>

/*
 * Test for compiler support of the __int128 integer data type. As of yet there
 * is no fallback for the lack of this functionality.
 */
#if ( !defined( __SIZEOF_INT128__) )
    #error "WiseMansFixedPoint: no support for extension __int128 detected."
#else
    static_assert( sizeof(void *) >= 8,
           "WiseMansFixedPoint: no support for extension __int128 detected.");
#endif

/*
 * We rely on signed right shift to be arithmetic. Signed right shift of
 * negative values are normally implementation defined behaviour and therefore
 * needs testing to assure correct functionality.
 */
static_assert( (-2ll >> 1) == -1ll,
        "WiseMansFixedPoint rely on signed right shifts to be arithmetic." );
static_assert( __int128(-2) >> 1 == __int128(-1),
        "WiseMansFixedPoint rely on signed right shifts to be arithmetic." );

/*
 * Helper functions for converting the buildin __int128 data types to
 * std::string objects.
 */
//static inline std::string int128_to_string_base10(__int128 num)
//{
//    // Result can always fit in 42 characters including negative sign and null
//    // terminating character. Start applying characters from the back of the
//    // buffer.
//    char buffer[42] = { 0 };
//    std::size_t idx{ 40 };
//    bool negative{ num < 0 };
//    while (num != 0)
//    {
//        buffer[idx--] = negative ? '0' - num%10 : '0' + num%10;
//        num /= 10;
//    }
//    if (negative)
//    {
//        buffer[idx--] = '-';
//    }
//    return std::string( &buffer[idx+1] );
//}
//
//static inline std::string uint128_to_string_base10(unsigned __int128 num)
//{
//    // Result can always fit in 41 characters.
//    char buffer[41] = { 0 };
//    std::size_t idx{ 39 };
//    while (num != 0)
//    {
//        buffer[idx--] = '0' + num%10;
//        num /= 10;
//    }
//    return std::string( &buffer[idx+1] );
//}

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

    /*
     * Helper function for retrieving a string of the fractional part on
     * quotient form.
     */
    std::string get_frac_quotient() const noexcept
    {
        uint64_t denominator{ 1ull << FRAC_BITS };
        uint64_t numerator{static_cast<uint64_t>(num & 0xFFFFFFFFFFFFFFFFull)};
        numerator >>= 64 - FRAC_BITS;
        return std::to_string(numerator) + "/" + std::to_string(denominator);
    }

    _128_INT_TYPE get_num_sign_extended() const noexcept
    {
        return (num << (64-INT_BITS)) >> 64-INT_BITS;
    }


protected:

    /*
     * Apply the bit mask to the internal representation to erase any bits
     * outside the intended numbers range.
     */
    void apply_bit_mask() noexcept
    {
        constexpr _128_INT_TYPE upper_mask = 1;
        constexpr _128_INT_TYPE lower_mask = 2;
    }

    /*
     * This method will round the result of some operation to the closesd fixed
     * point number in the current representation.
     */
    void round() noexcept
    {
        // Perform the rounding.
        num += _128_INT_TYPE(1) << (63-FRAC_BITS);
    }

    /*
     * Underlying data type. It will either be a signed or unsigned 128-bit
     * integer.
     */
    _128_INT_TYPE num{};

    /*
     * Protected constructor since base type should not be instantiatable.
     */
    BaseFixedPoint() = default;
};

/*
 * Signed fixed point data type.
 */
template <int INT_BITS, int FRAC_BITS>
class SignedFixedPoint : public BaseFixedPoint<
                         INT_BITS,FRAC_BITS,__int128>
{
public:
    explicit SignedFixedPoint(double a)
    {
        int64_t b = std::llround(a * static_cast<double>(1ll << 32));
        this->num = b; this->num <<= 32;
        this->num = this->get_num_sign_extended();
        this->round();
    }

    /*
     * Specialized to_string function for signed fixed point numbers.
     */
    std::string to_string() const noexcept
    {
        // Extract integer part.
        std::string integer{ 
            std::to_string(int64_t(this->get_num_sign_extended() >> 64)) 
        };

        // Append fractional part.
        return integer + " + " + this->get_frac_quotient();
    }
};

template <int INT_BITS, int FRAC_BITS>
class UnsignedFixedPoint : public BaseFixedPoint<
                           INT_BITS,FRAC_BITS,unsigned __int128>
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
