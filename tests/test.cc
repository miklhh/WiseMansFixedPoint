#include "catch.hpp"
#include "FixedPoint.h"
#include <sstream>
#include <iostream>
#include <cmath>
#include <chrono>

template <int INT_BITS, int FRAC_BITS>
using FixedPoint = SignedFixedPoint<INT_BITS, FRAC_BITS>;

TEST_CASE("Floating-point constructor")
{
    /*
     * Two simple tests.
     */
    {
        std::stringstream result{};
        SignedFixedPoint<10,10> fix_a{ 3.25 }, fix_b{ -19.125 };
        result << fix_a.to_string() << "|" << fix_b.to_string();
        REQUIRE(result.str() == std::string("3 + 256/1024|-20 + 896/1024"));
    }

    /*
     * More tests.
     */
    {
        std::stringstream result{};
        FixedPoint<8,8> fix_a{ -1.55555555 }, fix_b{ -0.555555555 };
        result << fix_a << "|" << fix_b;
        REQUIRE(result.str() == std::string("-2 + 114/256|-1 + 114/256"));
    }

    /*
     * Test Zero.
     */
    {
        std::stringstream result{};
        FixedPoint<12,12> fix_a{ 0.0 };
        result << fix_a;
        REQUIRE(result.str() == std::string("0 + 0/4096"));
    }

    /*
     * Test correct rounding close to zero.
     */
    {
        std::stringstream result{};
        FixedPoint<12,12> fix_a{ -0.0001 }, fix_b{ -0.0002 };
        result << fix_a << "|" << fix_b;
        REQUIRE(result.str() == std::string("0 + 0/4096|-1 + 4095/4096"));
    }

}

//TEST_CASE("Instance going out of scope should reset value.")
//{
//    for (int i=0; i<5; ++i)
//    {
//        FixedPoint<10,12> fix{ 0.0 };
//        fix += FixedPoint<10,12>(i);
//        fix /= FixedPoint<3,0>(2);
//        REQUIRE( static_cast<double>(fix) == static_cast<double>(i) / 2 );
//    }
//}
//
//TEST_CASE("Assigment of FixedPoint values (just need to compile).")
//{
//    FixedPoint<10,10> fix_a{}, fix_b{};
//    FixedPoint<10,10> fix_c = fix_a * fix_b;
//    fix_c = fix_a * fix_b;
//}
//
//TEST_CASE("Unary negation.")
//{
//    /*
//     * Test all 'sides'.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<12,12> fix_a{ 12.25 }, fix_b{ -17.05 }, fix_c{ 0.0 };
//        result << -fix_a << "|" << -fix_b << "|" << -fix_c;
//        REQUIRE(result.str() == std::string("-13 + 3072/4096|17 + 205/4096|0 + 0/4096"));
//    }
//
//    /*
//     * -FIX_MIN = FIX_MIN due to overflow.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<9,9> fix{ -256.0 };
//        result << -fix;
//        REQUIRE(result.str() == std::string("-256 + 0/512"));
//    }
//}
//
//
//TEST_CASE("Floating-point constructor")
//{
//    /*
//     * Two simple tests.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<10,10> fix_a{ 3.25 }, fix_b{ -19.125 };
//        result << fix_a << "|" << fix_b;
//        REQUIRE(result.str() == std::string("3 + 256/1024|-20 + 896/1024"));
//    }
//
//    /*
//     * More tests.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<8,8> fix_a{ -1.55555555 }, fix_b{ -0.555555555 };
//        result << fix_a << "|" << fix_b;
//        REQUIRE(result.str() == std::string("-2 + 114/256|-1 + 114/256"));
//    }
//
//    /*
//     * Test Zero.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<12,12> fix_a{ 0.0 };
//        result << fix_a;
//        REQUIRE(result.str() == std::string("0 + 0/4096"));
//    }
//
//    /*
//     * Test correct rounding close to zero.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<12,12> fix_a{ -0.0001 }, fix_b{ -0.0002 };
//        result << fix_a << "|" << fix_b;
//        REQUIRE(result.str() == std::string("0 + 0/4096|-1 + 4095/4096"));
//    }
//
//}
//
//TEST_CASE("Fixed point to floating point conversion introductory test.")
//{
//    FixedPoint<6,10> fix_a{ -5.25 };
//    FixedPoint<9,16> fix_b{ 2.33 };
//    REQUIRE(static_cast<double>(fix_a) == -5.25);
//    REQUIRE(std::abs(static_cast<double>(fix_b) -2.33) < 0.0001);
//}
//
//
//TEST_CASE("Addition tests")
//{
//    /*
//     * Introductory tests.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<10,10> fix_a{3.25};
//        FixedPoint<11,11> fix_b{7.50};
//        result << fix_a + fix_b;
//        REQUIRE(result.str() == std::string("10 + 768/1024"));
//    }
//    {
//        std::stringstream result{};
//        FixedPoint<10,10> fix_a{3.3333333}; // (3.33301) when rounded.
//        FixedPoint<10,10> fix_b{7.4444444}; // (4.44433) when rounded.
//        result << fix_a + fix_b;
//        REQUIRE(result.str() == std::string("10 + 796/1024"));
//    }
//
//    /*
//     * Addition with negative operands of different word lengths.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<9,10> fix_a{ -3.25 };
//        FixedPoint<6,12> fix_b{ 7.75 };
//        result <<  fix_a +  fix_b << "|";
//        result << (fix_a += fix_b);
//        REQUIRE(result.str() == std::string("4 + 512/1024|4 + 512/1024"));
//    }
//    {
//        std::stringstream result{};
//        FixedPoint<9,10> fix_a{ 3.25 };
//        FixedPoint<6,12> fix_b{ -7.75 };
//        result <<  fix_a +  fix_b << "|";
//        result << (fix_a += fix_b);
//        REQUIRE(result.str() == std::string("-5 + 512/1024|-5 + 512/1024"));
//    }
//
//    /*
//     * Subtraction with negative operands of different word lengths.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<9,10> fix_a{ -3.50 };
//        FixedPoint<6,12> fix_b{ 7.75 };
//        result <<  fix_a -  fix_b << "|";
//        result << (fix_a -= fix_b);
//        REQUIRE(result.str() == std::string("-12 + 768/1024|-12 + 768/1024"));
//    }
//    {
//        std::stringstream result{};
//        FixedPoint<9,10> fix_a{ 5.50 };
//        FixedPoint<6,12> fix_b{ -9.75 };
//        result <<  fix_a -  fix_b << "|";
//        result << (fix_a -= fix_b);
//        REQUIRE(result.str() == std::string("15 + 256/1024|15 + 256/1024"));
//    }
//}
//
//
//TEST_CASE("Multiplication of Fixed Point Numbers")
//{
//    /*
//     * Basic multiplication with (pos, pos), (pos, neg), (neg, pos) and (neg,
//     * neg)
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<8,4> fix_a{ 3.0 };
//        FixedPoint<7,8> fix_b{ 2.0 };
//        result << fix_b * fix_a;
//        REQUIRE(result.str() == std::string("6 + 0/4096"));
//    }
//    {
//        std::stringstream result{};
//        FixedPoint<10,10> fix_a{ 3.25 };
//        FixedPoint<12,12> fix_b{ 1.925 };
//        result << fix_b * fix_a;
//        REQUIRE(result.str() == std::string("6 + 1075456/4194304"));
//    }
//    {
//        std::stringstream result{};
//        FixedPoint<10,10> fix_a{ -7.02 }; // (-7.01953) when rounded.
//        FixedPoint<14,10> fix_b{ 1.925 }; // ( 1.92480) when rounded.
//        result << fix_b * fix_a;
//        REQUIRE(result.str() == std::string("-14 + 512516/1048576"));
//    }
//    {
//        std::stringstream result{};
//        FixedPoint<9,10> fix_a{ 3.25 };
//        FixedPoint<11,12> fix_b{ -1.925 };
//        result << fix_b * fix_a;
//        REQUIRE(result.str() == std::string("-7 + 3118848/4194304"));
//    }
//    {
//        std::stringstream result{};
//        FixedPoint<10,10> fix_a{ -3.25 };
//        FixedPoint<19,12> fix_b{ -1.925 };
//        result << fix_b * fix_a;
//        REQUIRE(result.str() == std::string("6 + 1075456/4194304"));
//    }
//
//    /*
//     * Mulitplication with zero should equal zero.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<10,8> fix_a{ -3.25 };
//        FixedPoint<12,4> fix_b{ 0 };
//        result << fix_b * fix_a;
//        REQUIRE(result.str() == std::string("0 + 0/4096"));
//    }
//
//    /*
//     * Multiplication when INT_BITS+FRAC_BITS > 32.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<25,21> fix_a{ 1050.239 };
//        FixedPoint<20,21> fix_b{  238.052 };
//        result << fix_b * fix_a;
//        REQUIRE(result.str() == std::string("250011 + 2123598665/4294967296"));
//    }
//    {
//        std::stringstream result{};
//        FixedPoint<27,21> fix_a{ 15.25 };
//        FixedPoint<18,17> fix_b{ -4.75 };
//        result << fix_b * fix_a;
//        REQUIRE(result.str() == std::string("-73 + 2415919104/4294967296"));
//    }
//    {
//        std::stringstream result{};
//        FixedPoint<27,21> fix_a{ -15.25 };
//        FixedPoint<18,17> fix_b{   4.75 };
//        result << fix_b * fix_a;
//        REQUIRE(result.str() == std::string("-73 + 2415919104/4294967296"));
//    }
//
//    /*
//     * Test of longest fixed point numbers.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<32,32> fix_a{ -78.55 };
//        FixedPoint<32,29> fix_b{ -99.75 };
//        result << fix_b * fix_a;
//        REQUIRE(result.str() == std::string("7835 + 1556925664/4294967296"));
//    }
//
//
//    /*
//     * Some more tests.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<22,5> fix_a{ 2 };
//        FixedPoint<27,9> fix_b{ 524288.0 };
//        result << fix_b * fix_a;
//        REQUIRE(result.str() == std::string("1048576 + 0/16384"));
//    }
//
//}
//
//TEST_CASE("Fixed point division")
//{
//    /*
//     * Simple introductory test.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<13,22> fix_a{ 7.60 };
//        FixedPoint<14,17> fix_b{ 3.40 };
//        result << fix_a/fix_b;
//        REQUIRE(result.str() == std::string("2 + 986891/4194304"));
//    }
//
//    /*
//     * Negative operands.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<6,23> fix_a{ -7.60 };
//        FixedPoint<5,20> fix_b{ 3.40 };
//        result << fix_a/fix_b;
//        REQUIRE(result.str() == std::string("-3 + 6414816/8388608"));
//    }
//    {
//        std::stringstream result{};
//        FixedPoint<6,23> fix_a{ 7.60 };
//        FixedPoint<5,20> fix_b{ -3.40 };
//        result << fix_a/fix_b;
//        REQUIRE(result.str() == std::string("-3 + 6414816/8388608"));
//    }
//    {
//        std::stringstream result{};
//        FixedPoint<10,23> fix_a{ -7.60 };
//        FixedPoint<5,25> fix_b{ -3.40 };
//        result << fix_a/fix_b;
//        REQUIRE(result.str() == std::string("2 + 1973790/8388608"));
//    }
//}
//
//TEST_CASE("Approximate pi using Leibniz formula")
//{
//    /*
//     * 10 000 000 iterations of Leibniz formula should result in a number close
//     * to pi, correct up to 7 decimals (including the 3) when correctly rounded.
//     */
//    const double pi = 3.1415926535;
//    const int ITERATIONS=10000000;
//
//    FixedPoint<4,32> pi_fixed{ 4.0 };
//    FixedPoint<32,0> divisor{ 3.0 };
//    for (int i=0; i<ITERATIONS; ++i)
//    {
//        if (i % 2)
//        {
//            // Odd iteration.
//            pi_fixed += FixedPoint<4,32>{4.0}/divisor;
//        }
//        else
//        {
//            // Even iteration.
//            pi_fixed -= FixedPoint<4,32>{4.0}/divisor;
//        }
//        divisor += FixedPoint<3,0>{2};
//    }
//
//    std::cout << std::endl;
//    std::cout << "Result from Leibniz formula of " << ITERATIONS << " iterations:" << std::endl;
//    std::cout.precision(9);
//    std::cout << "    Fixed     (fixed form)   : " << pi_fixed << std::endl;
//    std::cout << "    Fixed     (decimal form) : " << static_cast<double>(pi_fixed) << std::endl;
//    std::cout << "    Reference (decimal form) : " << pi << std::endl << std::endl;
//
//    // We can acquire around 7 significant digits using this method.
//    REQUIRE(std::abs(static_cast<double>(pi_fixed) - pi) < 0.000001);
//}
//
//
//TEST_CASE("Approximate e with Bernoulli limit.")
//{
//    /*
//     * A note on this approch. Since the sum of fractional bits of the left hand
//     * side and the right hand side of the multiplication is greater than 32, we
//     * are going to lose some precision at the for the performed product. There
//     * is no point of taking this algorithm further then it is already taken, it
//     * won't yield more significant digits.
//     */
//    const int ITERATIONS=99500;
//    const double e = 2.71828183;
//    FixedPoint<3,32> product_fixed{ 1.0 + 1.0/static_cast<double>(ITERATIONS) };
//    FixedPoint<3,32> e_fixed{ 1.0 };
//    for (int i=0; i<ITERATIONS; ++i)
//    {
//        e_fixed *= product_fixed;
//    }
//    std::cout << "Result from Bernoulli limit: n=" << ITERATIONS << std::endl;
//    std::cout << "    Fixed     (fixed form)   : " << e_fixed << std::endl;
//    std::cout << "    Fixed     (decimal form) : " << static_cast<double>(e_fixed) << std::endl;
//    std::cout << "    Reference (decimal form) : " << e << std::endl << std::endl;
//
//    // We can acquire around 6 significant digits using this method.
//    REQUIRE(std::abs(static_cast<double>(e_fixed) - e) < 0.00001);
//
//}
//
//TEST_CASE("Conversion from BIG fixed point numbers to floating point conversion.")
//{
//    FixedPoint<29,29> a{ 178956970, 357913941 };
//    double a_ref{ 178956970.66666666604 };
//    double a_error{ std::abs(a_ref - static_cast<double>(a)) };
//
//    FixedPoint<30,30> b{ 536870911, 178956970 };
//    double b_ref{ 536870911.16666666604 };
//    double b_error{ std::abs(b_ref - static_cast<double>(b)) };
//
//    FixedPoint<31,31> c{ -1073741823, 195225801 };
//    double c_ref{ -1073741822.9090909018 };
//    double c_error{ std::abs(c_ref - static_cast<double>(c)) };
//
//    std::cout << "Result from BIG fixed point to floating point conversion:" << std::endl;
//    std::cout.width(41); std::cout << "Reference";
//    std::cout.width(36); std::cout << "Fixed->Float      | Error" << std::endl;
//    std::cout.precision(20);
//    std::cout << "    FixedPoint<29,29>: "; std::cout.width(23);
//    std::cout << a_ref; std::cout.width(23); std::cout << static_cast<double>(a);
//    std::cout << " | "; std::cout.width(3); std::cout << a_error << std::endl;
//
//    std::cout << "    FixedPoint<30,30>: "; std::cout.width(23);
//    std::cout << b_ref; std::cout.width(23); std::cout << static_cast<double>(b);
//    std::cout << " | "; std::cout.width(3); std::cout << b_error << std::endl;
//
//    std::cout << "    FixedPoint<31,31>: "; std::cout.width(23);
//    std::cout << c_ref; std::cout.width(23); std::cout << static_cast<double>(c);
//    std::cout << " | "; std::cout.width(3); std::cout << c_error << std::endl;
//    std::cout << std::endl;
//
//    // Error should be very small.
//    double err_tol = 0.000001;
//    REQUIRE((a_error < err_tol && b_error < err_tol && c_error < err_tol));
//}
//
//TEST_CASE("Multiplication performance.")
//{
//    /*
//     * Multiplication with 'short' will result in multiplication scenario 1
//     * (described in FixedPoint.h) and multiplication with 'long' will result
//     * in multiplication scenario 2.
//     */
//    using namespace std::chrono;
//    const int ITERATIONS=10000000;
//    const FixedPoint<1,30> fix_factor{ 0.9999995 }; // Rounded to 0.99999949988.
//    std::cout << "Results from multiplication performance test:" << std::endl;
//    std::cout.precision(7);
//
//    /*
//     * Double-precision floating-point multiplication test.
//     */
//    {
//        const double float_factor{ 0.9999995 };
//        double float_res{ 0.9999995 };
//        auto t1 = high_resolution_clock::now();
//        for (int i=0; i<ITERATIONS; ++i)
//        {
//            float_res *= float_factor;
//        }
//        auto t2 = high_resolution_clock::now();
//        auto time = duration_cast<microseconds>(t2 - t1);
//        std::cout << "    Float res:       ";
//        std::cout << static_cast<double>(float_res) << " @ ";
//        std::cout << time.count() << "us" << std::endl;
//    }
//
//    /*
//     * Short multiplication, scenario 1.
//     */
//    {
//        FixedPoint<1,30> fix_short{ 0.9999995 };
//        auto t1 = high_resolution_clock::now();
//        for (int i=0; i<ITERATIONS; ++i)
//        {
//            fix_short *= fix_factor;
//        }
//        auto t2 = high_resolution_clock::now();
//        auto time = duration_cast<microseconds>(t2 - t1);
//        std::cout << "    Short fixed res: ";
//        std::cout << static_cast<double>(fix_short) << " @ ";
//        std::cout << time.count() << "us" << std::endl;
//    }
//
//    /*
//     * Long multiplication, scenario 2.
//     */
//    {
//        FixedPoint<3,30> fix_long{ 0.9999995 };
//        auto t1 = high_resolution_clock::now();
//        for (int i=0; i<ITERATIONS; ++i)
//        {
//            fix_long *= fix_factor;
//        }
//        auto t2 = high_resolution_clock::now();
//        auto time = duration_cast<microseconds>(t2 - t1);
//        std::cout << "    Long fixed res:  ";
//        std::cout << static_cast<double>(fix_long) << " @ ";
//        std::cout << time.count() << "us" << std::endl;
//    }
//
//}
//
//TEST_CASE("Simple comparison test.")
//{
//   FixedPoint<10,10> a { 5.125 };
//   FixedPoint<20,23> b { 5.125 };
//   FixedPoint<4,20>  c { 2.0095 };
//   FixedPoint<2,3>   d { 0.0 };
//   FixedPoint<9,7>   e { 0.0 };
//   REQUIRE(a == b);
//   REQUIRE( (c < a && c <= a && c < b && c <= b) );
//   REQUIRE( (!(d < e) && !(d > e)) );
//   REQUIRE( (d <= e && d >= e) );
//   REQUIRE( (a >= d && a > e) );
//}
//
//TEST_CASE("Assignment where the wordlength changes.")
//{
//    /*
//     * Test constructors.
//     */
//    {
//        FixedPoint<10,10> a{ 4.75 }, b{ -13.0625 };
//        FixedPoint<14,14> a_longer{a}, b_longer{b};
//        REQUIRE( (a == a_longer && b == b_longer) );
//    }
//    {
//        FixedPoint<14,14> a{ 4.75 }, b{ -13.0625 };
//        FixedPoint<10,10> a_longer{a}, b_longer{b};
//        REQUIRE( (a == a_longer && b == b_longer) );
//    }
//
//    /*
//     * Test assignment operators.
//     */
//    {
//        FixedPoint<10,10> a{ 4.75 }, b{ -13.0625 };
//        FixedPoint<14,14> a_longer{}, b_longer{};
//        a_longer = a;
//        b_longer = b;
//        REQUIRE( (a == a_longer && b == b_longer) );
//        a = FixedPoint<1,0>(0); a = a_longer;
//        b = FixedPoint<1,0>(0); b = b_longer;
//        REQUIRE( (a == a_longer && b == b_longer) );
//    }
//}
//
//TEST_CASE("Overflow tests")
//{
//    {
//        std::stringstream result{};
//        FixedPoint<5,5> fix_a{ 10.0 };
//        FixedPoint<1,2> fix_b{ 0.25 };
//
//        // Result = 40, trucanted to 8.
//        result << fix_a / fix_b;
//        REQUIRE( result.str() == std::string("8 + 0/32") );
//    }
//    {
//        std::stringstream result{};
//        FixedPoint<10,10> fix_a{ -511.0 };
//        FixedPoint<10,10> fix_b{    1.0 };
//        FixedPoint<10,10> fix_c{    2.0 };
//
//        result << (fix_a - fix_b) << "|" << (fix_a - fix_c);
//        REQUIRE( result.str() == std::string("-512 + 0/1024|511 + 0/1024") );
//    }
//}

