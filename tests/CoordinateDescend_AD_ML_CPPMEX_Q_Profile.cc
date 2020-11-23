/*
 * A profiling test program used to measure the performance of 
 * WiseMansFixedPoint. It is a model of an ASIC implementation that evaluates
 * coordinate descent to minimize the maximum liklihood error estimate of 
 * channel strengths in an block-fading channel model for some large numbers of
 * anthennas (a massive MIMO scenario). The program generates its own test data
 * and then performes 10 iterations of the algorithm.
 *
 * ----------------------------------------------------------------------------
 *
 * Compile with MATLAB-Mex by invoking: 
 * mex -R2018a CXXFLAGS='$CXXFLAGS -std=c++14' CoordinateDescend_AD_ML_CPPMEX_Q.cc
 * from MATLAB command line. This unfortunately seems to be as slow/slower than
 * the current MATLAB implementation of the ML Coordinate Descent algorithm. 
 * 
 * Update 2020-02-03: After marking the complex arithmetic functions static 
 * inline, the C++ implementation now performes almost twice as good as the
 * MATLAB implementation.
 * 
 * Author: Mikael Henriksson
 */

#include "FixedPoint.h"

#include <vector>
#include <numeric>
#include <algorithm>
#include <random>
#include <sstream>
#include <cstdlib>
#include <random>
#include <complex>
#include <iostream>
#include <chrono>

#define COORDINATE_DESCEND_ITERATIONS 10

template <int INT_BITS,int FRAC_BITS>
using FixedPoint = SignedFixedPoint<INT_BITS,FRAC_BITS>;

/*
 * Some word length settings.
 */

// <6,25>
constexpr int S2_INT = 6;
constexpr int S2_FRAC = 25;

// <5,31>
constexpr int Sigma_inv_INT = 5;
constexpr int Sigma_inv_FRAC = 31;

// <10,16>
constexpr int GAMMA_HAT_INT = 10;
constexpr int GAMMA_HAT_FRAC = 16;


/*
 * Coordinate descent routine.
 */
static void CoordinateDescent(

    const std::vector<std::vector<int>> &A_exp,
    const std::vector<std::vector<FixedPoint<10,40>>> &Sigma_Y_real,
    const std::vector<std::vector<FixedPoint<10,40>>> &Sigma_Y_imag,
    std::vector<FixedPoint<S2_INT,S2_FRAC>> &s2_real,
    std::vector<FixedPoint<S2_INT,S2_FRAC>> &s2_imag,
    std::vector<FixedPoint<GAMMA_HAT_INT,GAMMA_HAT_FRAC>> &gamma_hat,
    std::vector<std::vector<FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC>>> &Sigma_inv_real,
    std::vector<std::vector<FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC>>> &Sigma_inv_imag,
    const int Kc, const int L, const int i)

{
    // Iteration profiling.
    using namespace std::chrono;
    auto t1 = high_resolution_clock::now();


    for (int k_acc = 0; k_acc < Kc; ++k_acc)
    {
        // Psudo randomize a user 'k' with the module random algorithm.
        int k = (2*i+1)*k_acc;
        k &= (1 << 11) - 1;

        // Column vector s2.
        for (int row=0; row<L; ++row)
        {
            s2_real[row] = FixedPoint<S2_INT,S2_FRAC>(0);
            s2_imag[row] = FixedPoint<S2_INT,S2_FRAC>(0);
            for (int j=0; j<L; ++j)
            {
                FixedPoint<5,25> Sigma_inv_real_short{ Sigma_inv_real[row][j] };
                FixedPoint<5,25> Sigma_inv_imag_short{ Sigma_inv_imag[row][j] };
                switch (A_exp[k][j])
                {
                    case 0: /* exp(i*pi*0) */
                        s2_real[row] += Sigma_inv_real_short;
                        s2_imag[row] += Sigma_inv_imag_short; 
                        break;
                    case 1:  /* exp(i*pi*1) */
                        s2_real[row] -= Sigma_inv_imag_short;
                        s2_imag[row] += Sigma_inv_real_short;
                        break;
                    case 2: /* exp(i*pi*2) */
                        s2_real[row] -= Sigma_inv_real_short;
                        s2_imag[row] -= Sigma_inv_imag_short;
                        break;
                    case 3: /* exp(i*pi*3) */
                        s2_real[row] += Sigma_inv_imag_short;
                        s2_imag[row] -= Sigma_inv_real_short;
                        break;
                }
            }
        }  

        std::vector<FixedPoint<6,13>> s2_real_tmp(L);
        std::vector<FixedPoint<6,13>> s2_imag_tmp(L);
        for (int i=0; i<L; i++)
        {
            s2_real_tmp[i] = s2_real[i];
            s2_imag_tmp[i] = s2_imag[i];
        }

        // Scalar N1. 
        FixedPoint<23,14>  N1{ 0 };
        for (int col=0; col<L; ++col)
        {
            // 6 15
            FixedPoint<19,15> col_acc_real{};
            FixedPoint<19,15> col_acc_imag{};
            for (int row=0; row<L; ++row)
            {
                // 6 13
                col_acc_real += s2_real_tmp[row] * Sigma_Y_real[col][row] +
                                s2_imag_tmp[row] * Sigma_Y_imag[col][row];
                col_acc_imag += s2_real_tmp[row] * Sigma_Y_imag[col][row] -
                                s2_imag_tmp[row] * Sigma_Y_real[col][row];
            }
            N1 += FixedPoint<19,12>(col_acc_real)*s2_real_tmp[col] -
                  FixedPoint<19,12>(col_acc_imag)*s2_imag_tmp[col];
        }


        // Scalar N2. 11.20 (11.25) 20 is ok.
        FixedPoint<11,20> N2{ 0 };
        for (int j=0; j<L; ++j)
        {
            switch (A_exp[k][j])
            {
                case 0: N2 += FixedPoint<11,20>(s2_real[j]); break;
                case 1: N2 += FixedPoint<11,20>(s2_imag[j]); break;
                case 2: N2 -= FixedPoint<11,20>(s2_real[j]); break;
                case 3: N2 -= FixedPoint<11,20>(s2_imag[j]); break;
            }
        }

        //
        // gamma_hat[k] = max( (N1-N2)/N2^2, -gamma_hat[k] ) 
        // 21,14, update: 21,5 OR LONGER!
        FixedPoint<21,6> N22{ FixedPoint<11,9>(N2)*FixedPoint<11,9>(N2) };
        FixedPoint<10,6> val{};
        if ( N22 == FixedPoint<2,0>{0} )
        {
            // From extensive simulation we know that N1-N2 is always less than 
            // zero, so val should be min{FixedPoint}.
            val = FixedPoint<10,6>{ -512.0 };
        }
        else
            val = FixedPoint<10,6>{ (N1-N2)/N22 };


        // FixedPoint<10,6> val{ (N1-N2)/N22 };
        FixedPoint<10,16> delta{ 0 };
        if ( val + gamma_hat[k] < FixedPoint<2,0>{0} )
        {
            delta = -gamma_hat[k];
        }
        else
        {
            delta = val;
        }
        gamma_hat[k] = gamma_hat[k] + delta;


        // Scalar N4 and q. (N4 > 0 by sim) (q < 0 by sim)
        FixedPoint<14,10> N4{ delta*N2 + FixedPoint<2,0>( 1.0 ) };
        FixedPoint<15,26> q{};
        if ( N4 == FixedPoint<2,0>(0) )
            q = FixedPoint<15,26>{ -16384.0 };
        else
            q = FixedPoint<15,26>{ FixedPoint<15,26>(delta)/N4 };


        /*
         * Update Sigma_inv which is a hermitian matrix. A good implementation
         * should only have to do half of the multiplications.
         */
        for (int row=0; row<L; ++row)
        {
            // 3,28
            FixedPoint<3,28> row_val_real{ s2_real[row]*q };
            FixedPoint<3,28> row_val_imag{ s2_imag[row]*q };
            for (int col=row; col<L; ++col)
            {
                FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC> value_real = 
                    Sigma_inv_real[row][col] - 
                    (
                      FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC>{
                        double(row_val_real) * double(s2_real[col]) +
                        double(row_val_imag) * double(s2_imag[col])
                      }
                    );
                FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC> value_imag = 
                    Sigma_inv_imag[row][col] - 
                    (
                      FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC>{ 
                        double(row_val_imag) * double(s2_real[col]) -
                        double(row_val_real) * double(s2_imag[col])
                      }  
                    );
                Sigma_inv_real[row][col] = value_real;
                Sigma_inv_imag[row][col] = value_imag;
                if (col != row)
                {
                    Sigma_inv_real[col][row] =  value_real;
                    Sigma_inv_imag[col][row] = -value_imag;
                }
            }
        }
    }

    // Iteration profiling.
    auto t2 = high_resolution_clock::now();
    auto iterationTime = duration_cast<milliseconds>(t2 - t1);
    std::cout << "[" << i << "]" << " iteration time: " << iterationTime.count() << " ms." << std::endl;
}

/*
 * C++/Mex function entry point.
 */
int main()
{
    // Read noise variance (diagonal) matrix data.
    const double sigma_sqr = 0.1;

    // Read anthena count.
    const int M = 196;

    // Read number of potential users.
    const int Kc = 2048;

    // Read number of pilot dimensions.
    const int L = 100;

    // Initialize row vector gamma hat.
    using std::vector;
    vector<FixedPoint<GAMMA_HAT_INT,GAMMA_HAT_FRAC>> gamma_hat(Kc, FixedPoint<30,31>(0.0)); // 10 16

    vector< vector<int> > A_exp( Kc, vector<int>(L, 0) );
    for (unsigned row=0; row<L; ++row)
    {
        for (unsigned col=0; col<Kc; ++col)
        {
            // Generate random number between 0 and 3.
            A_exp[col][row] = std::rand() % 4;
        }
    }

    /*
     * Initialze Sigma_Y (sample covariance matrix). Sigma_Y should equal: 1/M * Y*Y'
     * where Y' denotes that matrix hermitian. For further optimization, we know that
     * Sigma_Y has to be Sigma_Y'. Note that, unlike every other matrix in this
     * implementation, Sigma_Y is stored in column major format, this since iterative
     * reads from Sigma_Y only happens over rows.
     */
    vector< vector<FixedPoint<10,40>> > Sigma_Y_real(
        L, vector<FixedPoint<10,40>>(L, FixedPoint<10,40>(0) ));
    vector< vector<FixedPoint<10,40>> > Sigma_Y_imag(
        L, vector<FixedPoint<10,40>>(L, FixedPoint<10,40>(0) ));
    std::uniform_real_distribution<double> Y_dist(-40, 40);
    std::default_random_engine re;
    for (unsigned col=0; col<L; ++col)
    {
        for (unsigned row=0; row<L; ++row)
        {
            std::complex<double> elm{ 0.0, 0.0 };
            for (unsigned i=0; i<M; ++i)
            {
                double real = Y_dist(re);
                double imag = Y_dist(re);
                std::complex<double> element{real, imag};
                std::complex<double> hermit{real, -imag};
                elm += element * hermit;
            }
            elm.real(elm.real() / M);
            elm.imag(elm.imag() / M);
            Sigma_Y_real[col][row] = FixedPoint<10,40>(elm.real());
            Sigma_Y_imag[col][row] = FixedPoint<10,40>(elm.imag());
        }
    }


    /*
     * Initialize Sigma_inv.
     */
    vector< vector<FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC>> > Sigma_inv_real(
        L, vector<FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC>>(L, FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC>(0) ));
    vector< vector<FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC>> > Sigma_inv_imag(
        L, vector<FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC>>(L, FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC>(0) ));
    for (int i=0; i<L; ++i)
    {
        Sigma_inv_real[i][i] = FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC>(1.0)/FixedPoint<Sigma_inv_INT,Sigma_inv_FRAC>(sigma_sqr);
    }

    /*
     * Local vectors used inside inner loop.
     */
    std::vector< FixedPoint<S2_INT, S2_FRAC> > s2_real(L, FixedPoint<S2_INT, S2_FRAC>(0) );
    std::vector< FixedPoint<S2_INT, S2_FRAC> > s2_imag(L, FixedPoint<S2_INT, S2_FRAC>(0) );

    /*
     * Perform the Coordinate descent algorithm.
     */
    std::cout << "Running..." << std::endl;
    for (int i=0; i<COORDINATE_DESCEND_ITERATIONS; ++i)
    {
        CoordinateDescent(
            /*A_real, A_imag,*/ A_exp, Sigma_Y_real, Sigma_Y_imag, s2_real, s2_imag, 
            /*s1_real, s1_imag,*/ gamma_hat, Sigma_inv_real, Sigma_inv_imag, 
            Kc, L, i);
    }
    std::cout << "Complete." << std::endl;


    // Return gamma_hat to MATLAB.
    for (int i=0; i<Kc; ++i)
        (void) gamma_hat[i];
}
