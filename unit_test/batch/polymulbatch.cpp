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
    uniform_int_distribution<uint32_t> Bgdist(0, lvl1param::Bg);
    uniform_int_distribution<uint32_t> Torus32dist(0, UINT32_MAX);
    constexpr int batch = 2;

    alignas(64) Polynomialn<TFHEpp::lvl1param, batch> a;

    for (int j = 0; j < batch; j++)
        for (int i = 0; i < lvl1param::n; i++)
           a[j][i] = Bgdist(engine) - lvl1param::Bg / 2;

    for (int j = 0; j < batch; j++)
        for (typename TFHEpp::lvl1param::T &i : a[j])
            i = Bgdist(engine) - lvl1param::Bg / 2;

    alignas(64) Polynomialn<TFHEpp::lvl1param, batch> b;
    for (int j = 0; j < batch; j++)
        for (typename TFHEpp::lvl1param::T &i : b[j]) i = Torus32dist(engine);

    alignas(64) Polynomialn<TFHEpp::lvl1param, batch> polymul;

    TFHEpp::PolyMulbatch<TFHEpp::lvl1param,batch>(polymul, a, b);

    Polynomialn<TFHEpp::lvl1param, batch> naieve = {};

    for (int k = 0; k < batch; k++)
        for (int i = 0; i < lvl1param::n; i++) {
            for (int j = 0; j <= i; j++)
                 naieve[k][i] += static_cast<int32_t>(a[k][j]) * b[k][i - j];
            for (int j = i + 1; j < lvl1param::n; j++)
                  naieve[k][i] -= static_cast<int32_t>(a[k][j]) * b[k][lvl1param::n + i - j];
        }

    for (int j = 0; j < batch; j++)
       for (int i = 0; i < lvl1param::n; i++) {
            auto diff = abs(static_cast<int32_t>(naieve[j][i] - polymul[j][i]));
            std::cout <<naieve[j][i] << " " <<  polymul[j][i] << " " << diff << std::endl;
            c_assert(diff <= 900000);
         }

    cout << "PolyMul Passed" << endl;

    return 0;
}