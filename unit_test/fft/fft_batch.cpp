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

    cout << "Start fft batch test." << endl;
    constexpr int batch = 2;
    

    Polynomialn<lvl1param, batch> a;
    for (int j = 0; j < batch; j++)
        for (typename TFHEpp::lvl1param::T &i : a[j]) i = Torus32dist(engine);
    PolynomialInFDn<lvl1param, batch> resfft;
    TFHEpp::TwistIFFTbatch<lvl1param, batch>(resfft, a);
    Polynomialn<lvl1param, batch> res;
    TFHEpp::TwistFFTbatch<lvl1param, batch>(res, resfft);
    for (int j = 0; j < batch; j++)
        for (int i = 0; i < lvl1param::n; i++) {
            auto b = abs(static_cast<int32_t>(a[j][i] - res[j][i]));
            cout << b << " ";
        //c_assert(b <= 1);
        }

    return 0;
}