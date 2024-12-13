#include <algorithm>
#include "c_assert.hpp"
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{

    random_device seed_gen;
    default_random_engine engine(seed_gen());

    uniform_int_distribution<uint32_t> Torus32dist(0, UINT32_MAX);

    cout << "Start LVL1 test." << endl;
    Polynomial<lvl1param> a;
    for (typename TFHEpp::lvl1param::T &i : a) i = Torus32dist(engine);
    PolynomialInFD<lvl1param> resfft;
    TFHEpp::TwistIFFT<lvl1param>(resfft, a);
    Polynomial<lvl1param> res;
    TFHEpp::TwistFFT<lvl1param>(res, resfft);
    for (int i = 0; i < lvl1param::n; i++) {
        auto b = abs(static_cast<int32_t>(a[i] - res[i]));
        cout << b;
        c_assert(a <= 1);
    }
    cout << "FFT Passed" << endl;


    return 0;
}