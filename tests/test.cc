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

    /*
     * Test correctness of construction from very small floating point numbers.
     */
    {
        std::stringstream result_a{}, result_b{}, result_c{};
        FixedPoint<10,63> fix_a{ 3.25261e-19 }; // approx.   3 * 2^(-63).
        FixedPoint<10,63> fix_b{ 2.17925e-17 }; // approx. 201 * 2^(-63).
        FixedPoint<10,56> fix_c{ 1.28231e-14 }; // approx. 924 * 2^(-56).

        result_a << fix_a; result_b << fix_b; result_c << fix_c;
        REQUIRE(result_a.str() == std::string("0 + 3/9223372036854775808"));
        REQUIRE(result_b.str() == std::string("0 + 201/9223372036854775808"));
        REQUIRE(result_c.str() == std::string("0 + 924/72057594037927936"));
    }

    /*
     * Test correctness of construction from bigger numbers. Note especially
     * that with the size of these floating point numbers, the precision is less
     * than one.
     */
    {
        std::stringstream result_a{}, result_b{}, result_c{};
        FixedPoint<64, 4> fix_a{ 9223372036854775807.0 };
        FixedPoint<59,10> fix_b{ 288230376151711743.97 };
        FixedPoint<42,11> fix_c{ -2199023255552.0 };

        result_a << fix_a; result_b << fix_b; result_c << fix_c;
        //REQUIRE(result_a.str() == std::string("9223372036854775807 + 0/16"));
        //REQUIRE(result_b.str() == std::string("288230376151711743 + 992/1024"));
        //REQUIRE(result_c.str() == std::string("-2199023255552 + 0/2048"));
    }

}


TEST_CASE("Addition tests")
{
    /*
     * Introductory tests.
     */
    {
        std::stringstream result{};
        FixedPoint<10,10> fix_a{3.25};
        FixedPoint<11,11> fix_b{7.50};
        result << fix_a + fix_b;
        REQUIRE(result.str() == std::string("10 + 1536/2048"));
    }
    {
        std::stringstream result{};
        FixedPoint<10,10> fix_a{3.3333333}; // (3.33301) when rounded.
        FixedPoint<10,10> fix_b{7.4444444}; // (4.44433) when rounded.
        result << fix_a + fix_b;
        REQUIRE(result.str() == std::string("10 + 796/1024"));
    }

    /*
     * Addition with negative operands of different word lengths.
     */
    {
        std::stringstream result{};
        FixedPoint<9,10> fix_a{ -3.25 };
        FixedPoint<6,12> fix_b{ 7.75 };
        result <<  fix_a + fix_b << "|";
        result << (fix_a += fix_b);
        REQUIRE(result.str() == std::string("4 + 2048/4096|4 + 512/1024"));
    }
    {
        std::stringstream result{};
        FixedPoint<9,10> fix_a{ 3.25 };
        FixedPoint<6,12> fix_b{ -7.75 };
        result <<  fix_a + fix_b << "|";
        result << (fix_a += fix_b);
        REQUIRE(result.str() == std::string("-5 + 2048/4096|-5 + 512/1024"));
    }

    /*
     * Subtraction with negative operands of different word lengths.
     */
    {
        std::stringstream result{};
        FixedPoint<9,10> fix_a{ -3.50 };
        FixedPoint<6,12> fix_b{ 7.75 };
        result <<  fix_a -  fix_b << "|";
        result << (fix_a -= fix_b);
        REQUIRE(result.str() == std::string("-12 + 3072/4096|-12 + 768/1024"));
    }
    {
        std::stringstream result{};
        FixedPoint<9,10> fix_a{ 5.50 };
        FixedPoint<6,12> fix_b{ -9.75 };
        result <<  fix_a -  fix_b << "|";
        result << (fix_a -= fix_b);
        REQUIRE(result.str() == std::string("15 + 1024/4096|15 + 256/1024"));
    }

    /*
     * Immediate result of addition and subtraction should not overflow.
     */
    {
        std::stringstream result{};
        FixedPoint<4,2> fix_a{ 7.25 };
        FixedPoint<5,1> fix_b{ 15.50 };
        result <<  fix_a +  fix_b;;
        REQUIRE(result.str() == std::string("22 + 3/4"));
    }
}


TEST_CASE("Overflow tests")
{
    {
        std::stringstream result{};
        FixedPoint<5,5> fix_a{ 10.0 };
        FixedPoint<1,2> fix_b{ 0.25 };

        // Result = 40, trucanted to 8.
        result << fix_a / fix_b;
        REQUIRE( result.str() == std::string("8 + 0/32") );
    }
    {
        std::stringstream result{};
        FixedPoint<10,10> fix_a{ -511.0 }, fix_b{ 1.0 }, fix_c{ 2.0 };
        FixedPoint<10,10> res_a{ fix_a - fix_b };
        FixedPoint<10,10> res_b{ fix_a - fix_c };
        result << res_a << "|" << res_b;
        REQUIRE( result.str() == std::string("-512 + 0/1024|511 + 0/1024") );
    }
}


TEST_CASE("Instance going out of scope should reset value.")
{
    for (int i=0; i<5; ++i)
    {
        FixedPoint<10,12> fix{ 0.0 };
        fix += FixedPoint<10,12>(i);
        fix /= FixedPoint<3,0>(2);
        REQUIRE(double(fix) == double(i) / 2.0) ;
    }
}


TEST_CASE("Assignment of FixedPoint that truncates.")
{
    {
        std::stringstream result{};
        FixedPoint<5,10> fix_a{ 15.0 };
        FixedPoint<3,10> fix_b{};
        fix_b = fix_a;
        result << fix_b;
        REQUIRE(result.str() == std::string("-1 + 0/1024"));
    }
}


TEST_CASE("Assigment of FixedPoint values (just need to compile).")
{
    FixedPoint<10,10> fix_a{}, fix_b{};
    FixedPoint<10,10> fix_c = fix_a * fix_b;
    fix_c = fix_a * fix_b;
}


TEST_CASE("Unary negation.")
{
    /*
     * Test all 'sides'.
     */
    {
        std::stringstream result{};
        FixedPoint<12,12> fix_a{ 12.25 }, fix_b{ -17.05 }, fix_c{ 0.0 };
        result << -fix_a << "|" << -fix_b << "|" << -fix_c;
        REQUIRE(result.str() == std::string("-13 + 3072/4096|17 + 205/4096|0 + 0/4096"));
    }

    /*
     * -FIX_MIN = FIX_MIN due to overflow.
     */
    {
        std::stringstream result{};
        FixedPoint<9,9> fix{ -256.0 };
        result << -fix;
        REQUIRE(result.str() == std::string("-256 + 0/512"));
    }
}


TEST_CASE("Fixed point to floating point conversion introductory test.")
{
    FixedPoint<6,10> fix_a{ -5.25 };
    FixedPoint<9,16> fix_b{ 2.33 };
    REQUIRE(static_cast<double>(fix_a) == -5.25);
    REQUIRE(std::abs(static_cast<double>(fix_b) -2.33) < 0.0001);
}


TEST_CASE("Multiplication of Fixed Point Numbers")
{
    /*
     * Basic multiplication with (pos, pos), (pos, neg), (neg, pos) and (neg,
     * neg)
     */
    {
        std::stringstream result{};
        FixedPoint<8,4> fix_a{ 3.0 };
        FixedPoint<7,8> fix_b{ 2.0 };
        result << fix_b * fix_a;
        REQUIRE(result.str() == std::string("6 + 0/4096"));
    }
    {
        std::stringstream result{};
        FixedPoint<10,10> fix_a{ 3.25 };
        FixedPoint<12,12> fix_b{ 1.925 };
        result << fix_b * fix_a;
        REQUIRE(result.str() == std::string("6 + 1075456/4194304"));
    }
    {
        std::stringstream result{};
        FixedPoint<10,10> fix_a{ -7.02 }; // (-7.01953) when rounded.
        FixedPoint<14,10> fix_b{ 1.925 }; // ( 1.92480) when rounded.
        result << fix_b * fix_a;
        REQUIRE(result.str() == std::string("-14 + 512516/1048576"));
    }
    {
        std::stringstream result{};
        FixedPoint<9,10> fix_a{ 3.25 };
        FixedPoint<11,12> fix_b{ -1.925 };
        result << fix_b * fix_a;
        REQUIRE(result.str() == std::string("-7 + 3118848/4194304"));
    }
    {
        std::stringstream result{};
        FixedPoint<10,10> fix_a{ -3.25 };
        FixedPoint<19,12> fix_b{ -1.925 };
        result << fix_b * fix_a;
        REQUIRE(result.str() == std::string("6 + 1075456/4194304"));
    }

    /*
     * Mulitplication with zero should equal zero.
     */
    {
        std::stringstream result{};
        FixedPoint<10,8> fix_a{ -3.25 };
        FixedPoint<12,4> fix_b{ 0 };
        result << fix_b * fix_a;
        REQUIRE(result.str() == std::string("0 + 0/4096"));
    }

    /*
     * Multiplication when INT_BITS+FRAC_BITS > 32.
     */
    {
        std::stringstream result{};
        FixedPoint<25,21> fix_a{ 1050.239 };
        FixedPoint<20,21> fix_b{  238.052 };
        result << fix_b * fix_a;
        REQUIRE(result.str() == std::string("250011 + 2174565033588/4398046511104"));
    }
    {
        std::stringstream result{};
        FixedPoint<27,21> fix_a{ 15.25 };
        FixedPoint<18,17> fix_b{ -4.75 };
        result << fix_b * fix_a;
        REQUIRE(result.str() == std::string("-73 + 154618822656/274877906944"));
    }
    {
        std::stringstream result{};
        FixedPoint<27,21> fix_a{ -15.25 };
        FixedPoint<18,17> fix_b{   4.75 };
        result << fix_b * fix_a;
        REQUIRE(result.str() == std::string("-73 + 154618822656/274877906944"));
    }

    /*
     * Test of longest fixed point numbers.
     */
    {
        std::stringstream result{};
        FixedPoint<32,32> fix_a{ -78.55 };
        FixedPoint<32,29> fix_b{ -99.75 };
        result << fix_b * fix_a;
        REQUIRE(result.str() == std::string("7835 + 835868101550538752/2305843009213693952"));
    }


    /*
     * Some more tests.
     */
    {
        std::stringstream result{};
        FixedPoint<22,5> fix_a{ 2 };
        FixedPoint<27,9> fix_b{ 524288.0 };
        result << fix_b * fix_a;
        REQUIRE(result.str() == std::string("1048576 + 0/16384"));
    }
}


TEST_CASE("Fixed point division")
{
    /*
     * Simple introductory test.
     */
    {
        std::stringstream result{};
        FixedPoint<13,22> fix_a{ 7.60 };
        FixedPoint<14,17> fix_b{ 3.40 };
        result << fix_a/fix_b;
        REQUIRE(result.str() == std::string("2 + 986890/4194304"));
    }

    /*
     * Negative operands.
     */
    {
        std::stringstream result{};
        FixedPoint<6,23> fix_a{ -7.60 };
        FixedPoint<5,20> fix_b{ 3.40 };
        result << fix_a/fix_b;
        REQUIRE(result.str() == std::string("-3 + 6414815/8388608"));
    }
    {
        std::stringstream result{};
        FixedPoint<6,23> fix_a{ 7.60 };
        FixedPoint<5,20> fix_b{ -3.40 };
        result << fix_a/fix_b;
        REQUIRE(result.str() == std::string("-3 + 6414815/8388608"));
    }
    {
        std::stringstream result{};
        FixedPoint<10,23> fix_a{ -7.60 };
        FixedPoint<5,25> fix_b{ -3.40 };
        result << fix_a/fix_b;
        REQUIRE(result.str() == std::string("2 + 1973790/8388608"));
    }
}


TEST_CASE("Rounding test.")
{
    {
        std::stringstream result{};
        FixedPoint<5,2> fix_a{ 10.5 };
        FixedPoint<5,5> fix_b{ 0.1875 };
        FixedPoint<5,2> res{ fix_a + fix_b };
        result << res;
        REQUIRE(result.str() == std::string("10 + 2/4"));
    }
    {
        std::stringstream result{};
        FixedPoint<5,2> fix_a{ 10.5 };
        FixedPoint<5,5> fix_b{ 0.1875 };
        FixedPoint<5,2> res{ rnd<5,2>(fix_a + fix_b) };
        result << res;
        REQUIRE(result.str() == std::string("10 + 3/4"));
    }
    {
        std::stringstream result{};
        FixedPoint<5,2> fix_a{ 10.5 };
        FixedPoint<5,5> fix_b{ 0.1875 };
        FixedPoint<5,2> res{};
        res.rnd( fix_a + fix_b );
        result << res;
        REQUIRE(result.str() == std::string("10 + 3/4"));
    }
}


TEST_CASE("Test of negative wordlengths.")
{
//    /*
//     * These fixed point numbers should, of course, act as any other fixed point 
//     * numbers and therefore should be covered by the other tests. Here we just
//     * test the different kinds of operator to see if the compiler generates any
//     * wierd compiler warnings due to negative shifts or likewise.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<-4,10> fix_a{ 0.125 / 4.0 };
//        fix_a + fix_a;
//        fix_a - fix_a;
//        fix_a * fix_a;
//        fix_a / fix_a;
//        fix_a += fix_a;
//        fix_a -= fix_a;
//        fix_a *= fix_a;
//        fix_a /= fix_a;
//        result << fix_a;
//    }
//    {
//        std::stringstream result{};
//        FixedPoint<10,-4> fix_a{ 32 };
//        fix_a + fix_a;
//        fix_a - fix_a;
//        fix_a * fix_a;
//        fix_a / fix_a;
//        fix_a += fix_a;
//        fix_a -= fix_a;
//        fix_a *= fix_a;
//        fix_a /= fix_a;
//        result << fix_a;
//    }
//
//
//    /*
//     * Bonus tests just for fun.
//     */
//    {
//        std::stringstream result{};
//        FixedPoint<6,-2> fix_a{ 16 };
//        FixedPoint<6,0> fix_b{ 3 };
//        FixedPoint<6,-2> fix_res{};
//        fix_res = fix_a + fix_b;
//        result << fix_res;
//        fix_res = rnd<6,-2>(fix_a + fix_b);
//        result << "|" << fix_res;
//        REQUIRE(result.str() == std::string("16|20"));
//    }
}


TEST_CASE("Saturation testing.")
{
    {
        std::stringstream result{};
        FixedPoint<10,5> fix_a{ 204.96875 };
        FixedPoint<10,5> fix_b{ -93.96875 };
        result << sat<4,4>(fix_a) << "|" << sat<4,5>(fix_a);
        REQUIRE(result.str() == std::string("7 + 15/16|7 + 31/32"));
        result.str(std::string()); result.clear();
        result << sat<5,4>(fix_a) << "|" << sat<5,5>(fix_a);
        REQUIRE(result.str() == std::string("15 + 15/16|15 + 31/32"));
        result.str(std::string()); result.clear();
        result << sat<4,4>(fix_b) << "|" << sat<4,5>(fix_b);
        REQUIRE(result.str() == std::string("-8 + 0/16|-8 + 0/32"));
        result.str(std::string()); result.clear();
        result << sat<5,4>(fix_b) << "|" << sat<5,5>(fix_b);
        REQUIRE(result.str() == std::string("-16 + 0/16|-16 + 0/32"));
    }
}


TEST_CASE("Breakage of non correctly sign extended numbers.")
{
    /*
     * This test only exists in due to a bug found in early code which happend
     * due to some optimizations that was 'to good' to work.
     */
    std::stringstream result{};
    FixedPoint<3,3> fix_a{ 24.00625 };  // Truncated to 0.0
    FixedPoint<10,4> fix_b{ 13.5 };
    result << (fix_b += fix_a);
    REQUIRE(result.str() == std::string("13 + 8/16"));
}


TEST_CASE("Approximate pi using Leibniz formula")
{
    /*
     * Time this test.
     */
    auto t1 = std::chrono::high_resolution_clock::now();

    /*
     * 950 000 iterations of Leibniz formula should result in a number close
     * to pi, correct up to 7 decimals (including the 3).
     */
    const double pi = 3.1415926535;
    const int ITERATIONS=950000;


    FixedPoint<4,63> pi_fixed{ 4.0 };
    FixedPoint<32,0> divisor{ 3.0 };
    for (int i=0; i<ITERATIONS; ++i)
    {
        if (i % 2)
        {
            // Odd iteration.
            pi_fixed += FixedPoint<4,32>{4.0}/divisor;
        }
        else
        {
            // Even iteration.
            pi_fixed -= FixedPoint<4,32>{4.0}/divisor;
        }
        divisor += FixedPoint<3,0>{2};
    }

    std::cout << std::endl;
    std::cout << "Result from Leibniz formula of " << ITERATIONS << " iterations:" << std::endl;
    std::cout.precision(9);
    std::cout << "    Fixed     (fixed form)   : " << pi_fixed << std::endl;
    std::cout << "    Fixed     (decimal form) : " << static_cast<double>(pi_fixed) << std::endl;
    std::cout << "    Reference (decimal form) : " << pi << std::endl;

    // We can acquire around 7 significant digits using this method.
    REQUIRE(std::abs(static_cast<double>(pi_fixed) - pi) < 0.000001);

    /*
     * Time this test.
     */
    auto t2 = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "    Time to run this test    : " << time.count();
    std::cout << " ms" << std::endl << std::endl;
}


TEST_CASE("Approximate e with Bernoulli limit.")
{
    /*
     * Time this test.
     */
    auto t1 = std::chrono::high_resolution_clock::now();

    /*
     * Note especially that iterations is choosen to be a power of two such that
     * iterations^(-1) will fit exactly into a fixed point number.
     */
    const double e = 2.71828183;
    const int ITERATIONS = 1048574;
    FixedPoint<3,20> product_fixed{ 1.0 + 1.0/static_cast<double>(ITERATIONS) };
    FixedPoint<3,43> e_fixed{ 1.0 };
    for (int i=0; i<ITERATIONS; ++i)
    {
        e_fixed *= product_fixed;
    }
    std::cout << "Result from Bernoulli limit: n=" << ITERATIONS << std::endl;
    std::cout << "    Fixed     (fixed form)   : " << e_fixed << std::endl;
    std::cout << "    Fixed     (decimal form) : " << static_cast<double>(e_fixed) << std::endl;
    std::cout << "    Reference (decimal form) : " << e << std::endl;

    // We can acquire around 6 significant digits using this method.
    REQUIRE(std::abs(static_cast<double>(e_fixed) - e) < 0.00001);

    /*
     * Time this test.
     */
    auto t2 = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << "    Time to run this test    : " << time.count();
    std::cout << " ms" << std::endl << std::endl;

}


TEST_CASE("Addition performance.")
{
    /*
     * Add up the fibbonacci sequence.
     */
    using namespace std::chrono;
    const int ITERATIONS=1000000;
    //double factor = 0.999995;
    double factor = 1.0 - 5.0/std::pow(2, 20);
    const FixedPoint<1,20> fix_factor{ factor };
    std::cout << "Results from addition performance test:" << std::endl;
    std::cout.precision(7);

    /*
     * Double-precision floating-point addition test.
     */
    std::chrono::microseconds time_float{};
    {
        double float_res{ 0.0 };

        auto t1 = high_resolution_clock::now();
        for (int i=0; i<ITERATIONS; ++i)
        {
            float_res += 12.345678;
            float_res -= 12.345677;
        }
        auto t2 = high_resolution_clock::now();

        time_float = duration_cast<microseconds>(t2 - t1);
        std::cout << "    Float res:       ";
        std::cout << static_cast<double>(float_res) << " @ ";
        std::cout << time_float.count() << "us" << std::endl;
    }

    /*
     * Addition test using fixed point.
     */
    std::chrono::microseconds time_fix{};
    {
        FixedPoint<63,63> fix_res{ 0.0 };

        auto t1 = high_resolution_clock::now();
        for (int i=0; i<ITERATIONS; ++i)
        {
            fix_res += FixedPoint<63,63>{ 12.345678 };
            fix_res -= FixedPoint<63,63>{ 12.345677 };
        }
        auto t2 = high_resolution_clock::now();

        time_fix = duration_cast<microseconds>(t2 - t1);
        std::cout << "    Fixed point res: ";
        std::cout << static_cast<double>(fix_res) << " @ ";
        std::cout << time_fix.count() << "us" << std::endl;
    }

    /*
     * Display difference in addition timings.
     */
    std::cout << "    Fix point addition slower with a factor: ";
    std::cout << double(time_fix.count()) / double(time_float.count());
    std::cout << std::endl << std::endl;
}


TEST_CASE("Multiplication performance.")
{
    /*
     * Multiplication with 'short' will result in multiplication scenario 1
     * (described in FixedPoint.h) and multiplication with 'long' will result
     * in multiplication scenario 2.
     */
    using namespace std::chrono;
    const int ITERATIONS=1000000;
    double factor = 1.0 - 5.0/std::pow(2, 20);
    const FixedPoint<1,20> fix_factor{ factor };
    std::cout << "Results from multiplication performance test:" << std::endl;
    std::cout.precision(7);

    /*
     * Double-precision floating-point multiplication test.
     */
    std::chrono::microseconds time_float{};
    {
        const double float_factor{ factor };
        double float_res{ factor };
        auto t1 = high_resolution_clock::now();
        for (int i=0; i<ITERATIONS; ++i)
        {
            float_res *= float_factor;
        }
        auto t2 = high_resolution_clock::now();
        time_float = duration_cast<microseconds>(t2 - t1);
        std::cout << "    Float res:       ";
        std::cout << static_cast<double>(float_res) << " @ ";
        std::cout << time_float.count() << "us" << std::endl;
    }

    /*
     * Multiplication test using fixed point.
     */
    std::chrono::microseconds time_fix{};
    {
        FixedPoint<3,43> fix_long{ factor };
        auto t1 = high_resolution_clock::now();
        for (int i=0; i<ITERATIONS; ++i)
        {
            fix_long *= fix_factor;
        }
        auto t2 = high_resolution_clock::now();
        time_fix = duration_cast<microseconds>(t2 - t1);
        std::cout << "    Fixed point res: ";
        std::cout << static_cast<double>(fix_long) << " @ ";
        std::cout << time_fix.count() << "us" << std::endl;
    }

    /*
     * Display difference in multiplication timings.
     */
    std::cout << "    Fix point multiplication slower with a factor: ";
    std::cout << double(time_fix.count()) / double(time_float.count()) << std::endl;
}


TEST_CASE("Simple comparison test.")
{
   FixedPoint<10,10> a { 5.125 };
   FixedPoint<20,23> b { 5.125 };
   FixedPoint<4,20>  c { 2.0095 };
   FixedPoint<2,3>   d { 0.0 };
   FixedPoint<9,7>   e { 0.0 };
   REQUIRE(a == b);
   REQUIRE( (c < a && c <= a && c < b && c <= b) );
   REQUIRE( (!(d < e) && !(d > e)) );
   REQUIRE( (d <= e && d >= e) );
   REQUIRE( (a >= d && a > e) );
}


TEST_CASE("Assignment where the wordlength changes.")
{
    /*
     * Test constructors.
     */
    {
        FixedPoint<10,10> a{ 4.75 }, b{ -13.0625 };
        FixedPoint<14,14> a_longer{a}, b_longer{b};
        REQUIRE( (a == a_longer && b == b_longer) );
    }
    {
        FixedPoint<14,14> a{ 4.75 }, b{ -13.0625 };
        FixedPoint<10,10> a_longer{a}, b_longer{b};
        REQUIRE( (a == a_longer && b == b_longer) );
    }

    /*
     * Test assignment operators.
     */
    {
        FixedPoint<10,10> a{ 4.75 }, b{ -13.0625 };
        FixedPoint<14,14> a_longer{}, b_longer{};
        a_longer = a;
        b_longer = b;
        REQUIRE( (a == a_longer && b == b_longer) );
        a = FixedPoint<1,0>{0}; a = a_longer;
        b = FixedPoint<1,0>{0}; b = b_longer;
        REQUIRE( (a == a_longer && b == b_longer) );
    }
}

