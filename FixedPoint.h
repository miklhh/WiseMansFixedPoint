/*
 * WiseMansFixedPoint header only C++ fixed point module. This fixed point
 * implementation supports numbers of varying integer and fractional word
 * lengths, and it employs the most basic arithmetic functions +,-,*,/ with
 * proper rounding when desiered. It is created with a hardware implementation
 * of fixed point numbers as leader, and it strives to behave as close to the
 * VHDL IEEE signed and unsigned data types as possible.
 *
 * For the penalty of some greater run-time, the user can enable overflow checks
 * by compiling the header with preprocessor macro '_DEBUG_SHOW_OVERFLOW_INFO'
 * defined (commandline option '-D_DEBUG_SHOW_OVERFLOW_INFO' for GCC or CLANG)
 * which will display some over-/underflow information during code execution.
 *
 *
 * Author: Mikael Henriksson [www.github.com/miklhh]
 * Repo: https://github.com/miklhh/WiseMansFixedPoint
 */

#ifndef _WISE_MANS_FIXED_POINT_H
#define _WISE_MANS_FIXED_POINT_H

#include "ttmath/ttmath.h"

#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <type_traits>


/*
 * Enabling compiler flag _DEBUG_SHOW_OVERFLOW_INFO will help the user find
 * fixed point assignments that overflows. Overflow info will be printed through
 * the debug routine defined below.
 */
#ifdef _DEBUG_SHOW_OVERFLOW_INFO
    #include <sstream>

    #ifdef MX_API_VER
        /*
         * Running in MATLAB/Mex environment. Use mexPrint for printing
         * overflows.
         */
        #include "mex.h"
        void _DEBUG_PRINT_FUNC(const char *str)
        {
            mexPrintf("%s\n", str);
        }
    #else
        /*
         * Running in basic environment. Use stderr for printing overflows.
         */
        void _DEBUG_PRINT_FUNC(const char *str)
        {
            std::cerr << str << std::endl;
        }
    #endif // MX_API_VER

#endif // _DEBUG_SHOW_OVERFLOW_INFO


/*
 * Test for if constexpr support.
 */
#ifdef __cpp_if_constexpr
    #define CONSTEXPR constexpr
#else
    #define CONSTEXPR /* no if constexpr */
#endif


/*
 * Wide integer types. The 128 bit integers are used as underlying data type for
 * fixed point numbers and the 256 bit integers are used in the multiplication
 * and division operators.
 */
using int128_t = ttmath::Int<2>;
using int256_t = ttmath::Int<4>;
using uint128_t = ttmath::UInt<2>;
using uint256_t = ttmath::UInt<4>;


/*
 * Compile time template structures for retrieving the extended or narrowed
 * integer type of an underlying signed and unsigned integer type. Extension of
 * 64-bit numbers result in the signed or unsigned __(u)int128_t.
 */
template<typename T> struct extend_int {};
template<typename T> struct narrow_int {};
template<> struct extend_int<int128_t> { using type = int256_t; };
template<> struct extend_int<uint128_t> { using type = uint256_t; };
template<> struct narrow_int<int128_t> { using type = int64_t; };
template<> struct narrow_int<uint128_t> { using type = uint64_t; };
template<> struct extend_int<int64_t> { using type = __int128_t; };
template<> struct extend_int<uint64_t> { using type = __uint128_t; };


namespace detail
{
    /*
     * Fast integer log2(double) function.
     */
    static inline int ilog2_fast(double d)
    {
        int result;
        std::frexp(d, &result);
        return result-1;
    }

    /*
     * Constexpr function for generating a ttmath 128 bit data type with the
     * value 1 << N, where: 0 <= N < 128. N outside of that range is undefined
     * behaviour.
     */
    template<typename TTMATH_INT>
    constexpr TTMATH_INT ONE_SHL(int N)
    {
        TTMATH_INT res{};
        if (N >= 64)
        {
            res.table[0] = 0ull;
            res.table[1] = 1ull << (N-64);
        }
        else
        {
            res.table[0] = 1ull << N;
            res.table[1] = 0ull;
        }
        return ttmath::UInt<2>(res);
    }

    /*
     * Constexpr function for generating a ttmath 128 bit data type with the
     * value (1 << N) - 1, where: 0 <= N < 128. N outside of that range is
     * undefined behaviour.
     */
    template<typename TTMATH_INT>
    constexpr TTMATH_INT ONE_SHL_M1(int N)
    {
        TTMATH_INT res{};
        if (N >= 64)
        {
            res.table[0] = ~0ull;
            res.table[1] = (1ull << (N-64)) - 1;
        }
        else
        {
            res.table[0] = (1ull << N) - 1;
            res.table[1] = 0ull;
        }
        return res;
    }

    /*
     * Constexpr function for generating a ttmath 128 bit data type with the
     * value ~((1 << N) - 1), where: 0 <= N < 128. N outside of that range
     * causes undefined behaviour.
     */
    template<typename TTMATH_INT>
    constexpr TTMATH_INT ONE_SHL_M1_INV(int N)
    {
        TTMATH_INT res{};
        if (N >= 64)
        {
            res.table[0] = 0ull;
            res.table[1] = ~( (1ull << (N-64)) - 1 );
        }
        else
        {
            res.table[0] = ~( (1ull << N) - 1 );
            res.table[1] = ~0ull;
        }
        return res;
    }

    /*
     * Multiplication of type 64bit x 64bit -> 128 bit.
     */
    template <typename T>
    static inline typename extend_int<T>::type mul_64_to_128(T a, T b)
    {
        return static_cast<typename extend_int<T>::type>(a) * b;
    }

    /*
     * Constexpr bounded shift functions. Needed due to negative wordlengths
     * causing alot of left shift wider thatn 64 or smaller than zero. These
     * constexpr function will shift the (u)int64_t argument a, left/right, n
     * times when possible, otherwise return the propriate result.
     */
    template <typename T>
    constexpr T BOUND_SHL(T a, int n)
    {
        if (n >= 64)
            return 0;
        else if (n <= 0)
            return a;
        else
            return a << n;
    }
    template <typename T>
    constexpr T BOUND_SHR(T a, int n)
    {
        if (n >= 64)
            return 0;
        else if (n <= 0)
            return a;
        else
            return a >> n;
    }
}


/*
 * Simple TTMath to_string function that generates a hex string of the underlying
 * ttmaths 2's complement data.
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
template <unsigned long _N>
static inline std::string to_string_hex(const ttmath::UInt<_N> &a)
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
     * Friend declaration of addition and subtraction arithmetic operators on
     * fixed point numbers.
     */
    template<
        int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
        int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
    LHS<std::max(LHS_INT_BITS,RHS_INT_BITS)+1,
        std::max(LHS_FRAC_BITS,RHS_FRAC_BITS) >
    friend operator+(
        const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs,
        const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs);

    template<
        int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
        int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
    LHS<std::max(LHS_INT_BITS,RHS_INT_BITS)+1,
        std::max(LHS_FRAC_BITS,RHS_FRAC_BITS) >
    friend operator-(
        const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs,
        const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs);

    template<int _INT_BITS, int _FRAC_BITS, template<int,int> class RHS>
    friend RHS<_INT_BITS,_FRAC_BITS> operator-(
        const RHS<_INT_BITS,_FRAC_BITS> &rhs);


    /*
     * Friend declaration of multiplication and division arithmetic operators on
     * fixed point numbers.
     */
    template<
        int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
        int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
    friend LHS<LHS_INT_BITS+RHS_INT_BITS,LHS_FRAC_BITS+RHS_FRAC_BITS>
    operator*(
        const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs,
        const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs);

    template<
        int LHS_INT_BITS, int LHS_FRAC_BITS, template<int,int> class LHS,
        int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_INT_TYPE >
    friend LHS<LHS_INT_BITS+RHS_FRAC_BITS,LHS_FRAC_BITS-RHS_FRAC_BITS>
    operator/(
        const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs,
        const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs);


    /*
     * Friend declaration for rounding of fixed point numbers.
     */
    template<
        int LHS_INT_BITS, int LHS_FRAC_BITS,
        int RHS_INT_BITS, int RHS_FRAC_BITS, template<int,int> class RHS>
    friend RHS<LHS_INT_BITS,LHS_FRAC_BITS> rnd(
        const RHS<RHS_INT_BITS, RHS_FRAC_BITS> &rhs);


    /*
     * Functions for getting the underlying data type.
     */
    _128_INT_TYPE get_num() const noexcept { return num; }

    virtual _128_INT_TYPE get_num_sign_extended() const noexcept = 0;


    /*
     * Function for setting the underlying data type to its sign extended
     * representation. This could have performance benefits over using
     * 'this->num = this->get_num_sign_extended()', if INT_BITS > 0. Profiling
     * the Cooridnate Descent project shows an approx. 20% performance increase
     * using this optimization.
     */
    virtual void set_num_sign_extended() noexcept = 0;



    /*
     * Common masking and sign extension of the signed and unsigned fixed point
     * assignment operator. When compiler flag _DEBUG_SHOW_OVERFLOW_INFO is
     * enabled, the assignment operator will forward an overflow message string
     * to _DEBUG_PRINT_FUNC(char *).
     */
    template <int RHS_INT_BITS, int RHS_FRAC_BITS>
    void assignment_common()
    {
        if CONSTEXPR (INT_BITS < RHS_INT_BITS)
        {
            #ifdef _DEBUG_SHOW_OVERFLOW_INFO
            {
                // Debug overflow checks
                if (this->test_overflow())
                {
                    std::stringstream ss{};
                    ss << "Overflow in assignment ";
                    ss << "<" << RHS_INT_BITS << "," << RHS_FRAC_BITS << "> ";
                    ss << "--> " << "<" << INT_BITS << "," << FRAC_BITS << "> ";
                    ss << "of value: " << this->to_string() << " ";
                    this->set_num_sign_extended();
                    ss << "truncated to: " << this->to_string();
                    _DEBUG_PRINT_FUNC(ss.str().c_str());
                }
                else
                {
                    // Sign extend (possibly truncate) MSB side.
                    this->set_num_sign_extended();
                }
            }
            #else
            {
                // Sign extend (possibly truncate) MSB side.
                this->set_num_sign_extended();
            }
            #endif
        }

        // Truncate fractional bits if necessary.
        if CONSTEXPR (FRAC_BITS < RHS_FRAC_BITS)
        {
            this->apply_bit_mask_frac();
        }
    }



    /*
     * Test for overflow in the underlying data type. The function is overloaded
     * in
     */
    virtual bool test_overflow() const noexcept;


    /*
     * To string method for signed and unsigned fixed point numbers. The
     * resulting string will be on the form: "a + b/2^FRAC_BITS" where a is the
     * integer part and b is the fractional part.
     */
    std::string to_string() const noexcept
    {
        // Extract integer part.
        using narrow_int_type = typename narrow_int<int_type>::type;
        narrow_int_type integer = narrow_int_type((this->num >> 64).ToUInt());
        std::string integer_str( std::to_string(integer) );

        // Append fractional part if it exists.
        if CONSTEXPR (FRAC_BITS > 0)
        {
            return integer_str + " + " + this->get_frac_quotient();
        }
        else
        {
            return integer_str;
        }
    }


    /*
     * Explilcit conversion to double data type.
     */
    explicit operator double() const
    {
        /*
         * Truncate num to 64 bits, with as many fractional bits remaining as
         * possible without truncating any integer bits.
         */
        using std::min; using std::max;
        using narrow_int_type = typename narrow_int<int_type>::type;
        constexpr int SHIFT_WIDTH = max( min(64, 64-FRAC_BITS), INT_BITS );
        narrow_int_type num_small = (this->num >> SHIFT_WIDTH).ToUInt();
        return double(num_small) / double(1ull << (64-SHIFT_WIDTH));
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


protected:
    /*
     * Construct a fixed point number from a floating point number.
     */
    void construct_from_double(double a)
    {
        /*
         * NOTE: This magic number is (exactly) the greatest IEEE 754
         * Double-Precision Floating-Point number smaller than 1.0. It is
         * used to create a fast std::ceil(std::log2(std::abs(x)+1)) from
         * ilog2_fast(std::abs(x)+magic)+1.
         */
        constexpr double MAGIC_CEIL = 0.9999999999999999;
        long n = detail::ilog2_fast(std::abs(a) + MAGIC_CEIL) + 2;
        int64_t num = std::llround(a * double(1ull << (64-n)));
        this->num.table[0] = num << n;
        this->num.table[1] = num >> (64-n);
        this->round();
        this->apply_bit_mask_frac();
        this->set_num_sign_extended();
    }


    /*
     * Helper function for retrieving a string of the fixed point fractional
     * part on quotient form.
     */
    std::string get_frac_quotient() const noexcept
    {
        using std::to_string;
        int_type num_masked = num;
        num_masked.table[1] = 0;
        uint64_t numerator{ (num_masked >> (64-FRAC_BITS)).ToUInt() };
        if CONSTEXPR (FRAC_BITS > 0)
        {
            uint64_t denominator{ 1ull << FRAC_BITS };
            return to_string(numerator) + "/" + to_string(denominator);
        }
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
    void apply_bit_mask_frac() noexcept
    {
        this->num &= detail::ONE_SHL_M1_INV<int_type>(64-FRAC_BITS);
    }


    /*
     * Method for rounding the result of some operation to the closest fixed
     * point number in the current representation.
     */
    void round() noexcept
    {
        num += detail::ONE_SHL<int_type>(63-FRAC_BITS);
    }


    /*
     * Underlying data type. It will either be a signed or unsigned 128-bit
     * integer.
     */
    using int_type = _128_INT_TYPE;
    int_type num{};


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
     * Conversion from floating point numbers.
     */
    explicit SignedFixedPoint(double a) { this->construct_from_double(a); }


    /*
     * Copy assignment operator for signed fixed point numbers.
     */
    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_128_INT_TYPE>
    SignedFixedPoint<INT_BITS, FRAC_BITS> &
    operator=(const BaseFixedPoint<
        RHS_INT_BITS,RHS_FRAC_BITS, RHS_128_INT_TYPE > &rhs) noexcept
    {
        this->num = rhs.num;
        this->template assignment_common<RHS_INT_BITS,RHS_FRAC_BITS>();
        return *this;
    }


    /*
     * Copy constructor for signed fixed point numbers.
     */
    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_128_INT_TYPE>
    SignedFixedPoint(const BaseFixedPoint<
            RHS_INT_BITS,RHS_FRAC_BITS,RHS_128_INT_TYPE > &rhs) noexcept
    {
        // Reuse copy assignment.
        *this = rhs;
    }


    /*
     * Assignment from other fixed point number with proper rounding.
     */
    template <int RHS_INT_BITS,int RHS_FRAC_BITS, typename RHS_128_INT_TYPE>
    SignedFixedPoint<INT_BITS, FRAC_BITS> &
        rnd(const BaseFixedPoint<
            RHS_INT_BITS,RHS_FRAC_BITS,RHS_128_INT_TYPE > &rhs)
    {
        this->num = rhs.num;
        this->round();
        this->apply_bit_mask_frac();
        this->set_num_sign_extended();
        return *this;
    }


    /*
     * Friend declaration for saturation function.
     */
    template <
        int LHS_INT_BITS,int LHS_FRAC_BITS,int RHS_INT_BITS,int RHS_FRAC_BITS>
    friend SignedFixedPoint<LHS_INT_BITS, LHS_FRAC_BITS> sat(
            const SignedFixedPoint<RHS_INT_BITS, RHS_FRAC_BITS> &rhs);


    /*
     * Returns the internal num representation sign extended, that is, num with
     * all bits more significant than the sign bit set to the value of the sign
     * bit.
     */
    int128_t get_num_sign_extended() const noexcept override
    {
        if ( sign() )
            return this->num | detail::ONE_SHL_M1_INV<int128_t>(64+INT_BITS);
        else
            return this->num & detail::ONE_SHL_M1<int128_t>(64+INT_BITS);
    }


    /*
     * Set the internal num representation sign extended, that is, num with all
     * bits more significant than the sign bit set to the value of the sign bit.
     */
    void set_num_sign_extended() noexcept override
    {
        if CONSTEXPR (INT_BITS <= 0)
        {
            if ( sign() )
                this->num |= detail::ONE_SHL_M1_INV<int128_t>(64+INT_BITS);
            else
                this->num &= detail::ONE_SHL_M1<int128_t>(64+INT_BITS);
        }
        else
        {
            if (sign())
                this->num.table[1] |= ~((1ull << (INT_BITS-1)) - 1);
            else
                this->num.table[1] &= (1ull << (INT_BITS-1)) - 1;
        }
    }

    /*
     * Test if signed overflow has occured before possible sign extension. This
     * is used when _DEBUG_SHOW_OVERFLOW_INFO is enabled and in the saturion
     * function.
     */
    bool test_overflow() const noexcept override
    {
        using detail::ONE_SHL_M1_INV;
        constexpr uint128_t MASK = ONE_SHL_M1_INV<uint128_t>(64+INT_BITS-1);
        return !( (this->num & MASK) == 0 || (this->num & MASK) == MASK );
    }


private:
    /*
     * Get the sign of the number.
     */
    bool sign() const noexcept
    {
        if CONSTEXPR (INT_BITS <= 0)
            return this->num.table[0] & (1ull << (64+INT_BITS-1));
        else
            return this->num.table[1] & (1ull << (INT_BITS-1));
    }
};


/*
 * Unsigned fixed point data type.
 */
template <int INT_BITS, int FRAC_BITS>
class UnsignedFixedPoint : public BaseFixedPoint<INT_BITS,FRAC_BITS,uint128_t>
{
public:
    UnsignedFixedPoint() = default;


    /*
     * Conversion from floating point numbers.
     */
    explicit UnsignedFixedPoint(double a) { this->construct_from_double(a); }


    /*
     * Copy assignment operator for unsigned fixed point numbers.
     */
    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_128_INT_TYPE>
    UnsignedFixedPoint<INT_BITS, FRAC_BITS> &
    operator=(const BaseFixedPoint<
        RHS_INT_BITS,RHS_FRAC_BITS, RHS_128_INT_TYPE > &rhs) noexcept
    {
        this->num = rhs.num;
        this->template assignment_common<RHS_INT_BITS,RHS_FRAC_BITS>();
        return *this;
    }


    /*
     * Copy constructor for unsigned fixed point numbers.
     */
    template <int RHS_INT_BITS, int RHS_FRAC_BITS, typename RHS_128_INT_TYPE>
    UnsignedFixedPoint(const BaseFixedPoint<
            RHS_INT_BITS,RHS_FRAC_BITS,RHS_128_INT_TYPE > &rhs) noexcept
    {
        // Reuse copy assignment.
        *this = rhs;
    }


    /*
     * Returns the internal num representation sign extended. For unsigned
     * numbers this means returning num with all bits greater than the most
     * significat bit zeroed.
     */
    uint128_t get_num_sign_extended() const noexcept override
    {
        return this->num & detail::ONE_SHL_M1<int128_t>(64+INT_BITS);
    }


    /*
     * Set the internal num representation sign extended. For unsigned numbers
     * this means zeroing all bits greater than the most significant bit.
     */
    void set_num_sign_extended() noexcept override
    {
        if CONSTEXPR (INT_BITS <= 0)
        {
            this->num &= detail::ONE_SHL_M1<int128_t>(64+INT_BITS);
        }
        else
        {
            this->num.table[1] &= (1ull << INT_BITS) - 1;
        }
    }


    /*
     * Test if unsigned overflow has occured before possible sign extension.
     * This is used when _DEBUG_SHOW_OVERFLOW_INFO is enabled and in the
     * saturion function.
     */
    bool test_overflow() const noexcept override
    {
        using detail::ONE_SHL_M1_INV;
        constexpr uint128_t MASK = ONE_SHL_M1_INV<uint128_t>(64+INT_BITS);
        return (this->num & MASK) != 0;
    }
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
    // Disallow arithmetic between signed and unsigned numbers.
    static_assert(
        std::is_same<RHS_INT_TYPE, typename LHS<1,0>::int_type >::value,
        "Arithmetic between signed and unsigned fixed point types disallowed. "
        "Use explicit type conversion and convert LHS or RHS to a common type."
    );

    constexpr int RES_INT_BITS = std::max(LHS_INT_BITS,RHS_INT_BITS)+1;
    constexpr int RES_FRAC_BITS = std::max(LHS_FRAC_BITS,RHS_FRAC_BITS);
    LHS<RES_INT_BITS,RES_FRAC_BITS> res{};

    // No sign extension or masking needed due to correct word length. The
    // following code seems to be the most consistent way of generating addition
    // with carry (x86 instruction 'adc') throughout the tests.
    res.num.table[1] = lhs.num.table[1] + rhs.num.table[1];
    res.num.table[0] = lhs.num.table[0] + rhs.num.table[0];
    res.num.table[1] += res.num.table[0] < lhs.num.table[0];
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
    // Disallow arithmetic between signed and unsigned numbers.
    static_assert(
        std::is_same<RHS_INT_TYPE, typename LHS<1,0>::int_type >::value,
        "Arithmetic between signed and unsigned fixed point types disallowed. "
        "Use explicit type conversion and convert LHS or RHS to a common type."
    );

    constexpr int RES_INT_BITS = std::max(LHS_INT_BITS,RHS_INT_BITS)+1;
    constexpr int RES_FRAC_BITS = std::max(LHS_FRAC_BITS,RHS_FRAC_BITS);
    LHS<RES_INT_BITS,RES_FRAC_BITS> res{};

    // No sign extension or masking needed due to correct word length. The
    // following code seems to be the most consistent way of generating
    // subtraction with borrow (x86 instruction 'sbb') throughtout the tests.
    res.num.table[1] = lhs.num.table[1] - rhs.num.table[1];
    res.num.table[0] = lhs.num.table[0] - rhs.num.table[0];
    res.num.table[1] -= res.num.table[0] > lhs.num.table[0];
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
    // Disallow arithmetic between signed and unsigned numbers.
    static_assert(
        std::is_same<RHS_INT_TYPE, typename LHS<1,0>::int_type >::value,
        "Arithmetic between signed and unsigned fixed point types disallowed. "
        "Use explicit type conversion and convert LHS or RHS to a common type."
    );

    LHS<LHS_INT_BITS+RHS_INT_BITS, LHS_FRAC_BITS+RHS_FRAC_BITS> res{};

    /*
     * Specialized multiplication operator for when the resulting word length
     * is smaller than or equal to 64 bits. This specialized version is faster
     * to execute since it does not need to perform the wide multiplication.
     */
    constexpr int LHS_TOTAL_BITS = LHS_INT_BITS+LHS_FRAC_BITS;
    constexpr int RHS_TOTAL_BITS = RHS_INT_BITS+RHS_FRAC_BITS;
    if CONSTEXPR (LHS_TOTAL_BITS+RHS_TOTAL_BITS <= 64)
    {
        uint64_t lhs_int, lhs_frac, rhs_int, rhs_frac;
        lhs_int = detail::BOUND_SHL(lhs.num.table[1], LHS_FRAC_BITS);
        rhs_int = detail::BOUND_SHL(rhs.num.table[1], RHS_FRAC_BITS);
        lhs_frac = detail::BOUND_SHR(lhs.num.table[0], 64-LHS_FRAC_BITS);
        rhs_frac = detail::BOUND_SHR(rhs.num.table[0], 64-RHS_FRAC_BITS);

        using short_int = typename narrow_int<typename LHS<1,0>::int_type>::type;
        short_int rhs_short = rhs_int | rhs_frac;
        short_int lhs_short = lhs_int | lhs_frac;
        short_int res_short = lhs_short * rhs_short;
        if CONSTEXPR (LHS_FRAC_BITS <= 0 && RHS_FRAC_BITS > 0)
        {
            res.num.table[0] = res_short << (64-RHS_FRAC_BITS);
            res.num.table[1] = res_short >> RHS_FRAC_BITS;
        }
        else if CONSTEXPR (LHS_FRAC_BITS > 0 && RHS_FRAC_BITS <= 0)
        {
            res.num.table[0] = res_short << (64-LHS_FRAC_BITS);
            res.num.table[1] = res_short >> LHS_FRAC_BITS;
        }
        else
        {
            constexpr int SHIFT_WIDTH = LHS_FRAC_BITS + RHS_FRAC_BITS;
            res.num.table[0] = detail::BOUND_SHL(res_short, 64-SHIFT_WIDTH);
            res.num.table[1] = detail::BOUND_SHR(res_short, SHIFT_WIDTH);
        }
    }
    /*
     * Yet another specialized multiplication operator. This utilizes the
     * specialized 64x64->128 bit multiplication that most computers can perform
     * to get slightly high performance than the fully 128x128 bit
     * multiplication yields.
     */
    else if CONSTEXPR (LHS_TOTAL_BITS <= 64 && RHS_TOTAL_BITS <= 64)
    {
        int64_t lhs_short =
            detail::BOUND_SHR(lhs.num.table[0], 64-LHS_FRAC_BITS) |
            detail::BOUND_SHL(lhs.num.table[1], LHS_FRAC_BITS);
        int64_t rhs_short =
            detail::BOUND_SHR(rhs.num.table[0], 64-RHS_FRAC_BITS) |
            detail::BOUND_SHL(rhs.num.table[1], RHS_FRAC_BITS);
        __int128_t res_long = detail::mul_64_to_128<int64_t>(lhs_short, rhs_short);
        if CONSTEXPR (LHS_FRAC_BITS <= 0 && RHS_FRAC_BITS > 0)
        {
            res_long <<= 64 - RHS_FRAC_BITS;
        }
        else if CONSTEXPR (LHS_FRAC_BITS > 0 && RHS_FRAC_BITS <= 0)
        {
            res_long <<= 64 - LHS_FRAC_BITS;
        }
        else
        {
            res_long <<= 64 - RHS_FRAC_BITS - LHS_FRAC_BITS;
        }
        res.num.table[1] = res_long >> 64;
        res.num.table[0] = res_long;
    }
    /*
     * General base case multiplication. This is the slowest to execute, but it
     * is able to perform all multiplications with all different word lengths.
     */
    else
    {
        using int_type = typename extend_int<
            typename narrow_int<typename LHS<1,0>::int_type>::type >::type;
        int_type lhs_frac = lhs.num.table[0];
        int_type lhs_int = lhs.num.table[1];
        int_type lhs_long = lhs_int << 64 | lhs_frac;
        int_type rhs_frac = rhs.num.table[0];
        int_type rhs_int = rhs.num.table[1];
        int_type rhs_long = rhs_int << 64 | rhs_frac;
        lhs_long >>= (64 - LHS_FRAC_BITS);
        rhs_long >>= (64 - RHS_FRAC_BITS);
        int_type res_long = lhs_long * rhs_long;
        res_long <<= (64 - RHS_FRAC_BITS - LHS_FRAC_BITS);
        res.num.table[1] = res_long >> 64;
        res.num.table[0] = res_long;
    }
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
LHS<LHS_INT_BITS+RHS_FRAC_BITS,LHS_FRAC_BITS-RHS_FRAC_BITS>
operator/(const LHS<LHS_INT_BITS,LHS_FRAC_BITS> &lhs,
          const BaseFixedPoint<RHS_INT_BITS,RHS_FRAC_BITS,RHS_INT_TYPE> &rhs)
{
    // Disallow arithmetic between signed and unsigned numbers.
    static_assert(
        std::is_same<RHS_INT_TYPE, typename LHS<1,0>::int_type >::value,
        "Arithmetic between signed and unsigned fixed point types disallowed. "
        "Use explicit type conversion and convert LHS or RHS to a common type."
    );

    // Get the 128-bit underlying type.
    using int_type = typename extend_int<
        typename narrow_int<typename LHS<1,0>::int_type>::type >::type;
    LHS<LHS_INT_BITS+RHS_FRAC_BITS,LHS_FRAC_BITS-RHS_FRAC_BITS> res{};

    // Extract LHS num to 128-bit integer.
    int_type lhs_frac = lhs.num.table[0];
    int_type lhs_int = lhs.num.table[1];
    int_type lhs_long = lhs_int << 64 | lhs_frac;
    lhs_long <<= 64-LHS_INT_BITS;

    // Extract RHS num to 128-bit integer.
    int_type rhs_frac = rhs.num.table[0];
    int_type rhs_int = rhs.num.table[1];
    int_type rhs_long = rhs_int << 64 | rhs_frac;
    rhs_long >>= 64-RHS_FRAC_BITS;

    // Perform division and move result in place.
    int_type res_long = lhs_long / rhs_long;
    if CONSTEXPR ( 64-LHS_INT_BITS-RHS_FRAC_BITS > 0 )
    {
        res_long >>= 64-LHS_INT_BITS-RHS_FRAC_BITS;
    }
    else
    {
        res_long <<= 64-LHS_INT_BITS-RHS_FRAC_BITS;
    }

    // Back to ttmath integer.
    res.num.table[1] = res_long >> 64;
    res.num.table[0] = res_long;

    // Truncate fractional side and return.
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
    res.set_num_sign_extended();
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
    res.set_num_sign_extended();
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
    res.num = rhs.get_num_sign_extended();
    if (res.test_overflow())
    {
        if ( res.num < 0 )
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
        res.apply_bit_mask_frac();
        return res;
    }
}


/*
 * Print-out to C++ stream object on the form '<int> + <frac>/<2^<frac_bits>'.
 * Good for debuging'n'stuff.
 */
template <int INT_BITS, int FRAC_BITS, template<int,int> class RHS>
std::ostream &operator<<(
        std::ostream &os, const RHS<INT_BITS, FRAC_BITS> &rhs)
{
    return os << rhs.to_string();
}


#endif // Header guard ending.

