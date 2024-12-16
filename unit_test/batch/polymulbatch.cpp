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

    alignas(64) array<typename TFHEpp::lvl1param::T, lvl1param::n> a;
    for (int i = 0; i < lvl1param::n; i++)
        a[i] = Bgdist(engine) - lvl1param::Bg / 2;
    for (typename TFHEpp::lvl1param::T &i : a)
        i = Bgdist(engine) - lvl1param::Bg / 2;
    alignas(64) array<typename TFHEpp::lvl1param::T, lvl1param::n> b;
    for (typename TFHEpp::lvl1param::T &i : b) i = Torus32dist(engine);

    alignas(64) Polynomial<lvl1param> polymul;
    TFHEpp::PolyMul<lvl1param>(polymul, a, b);
    Polynomial<lvl1param> naieve = {};
    for (int i = 0; i < lvl1param::n; i++) {
        for (int j = 0; j <= i; j++)
             naieve[i] += static_cast<int32_t>(a[j]) * b[i - j];
        for (int j = i + 1; j < lvl1param::n; j++)
              naieve[i] -= static_cast<int32_t>(a[j]) * b[lvl1param::n + i - j];
    }
    for (int i = 0; i < lvl1param::n; i++) {
        auto diff = abs(static_cast<int32_t>(naieve[i] - polymul[i]));
        std::cout <<naieve[i] << " " <<  polymul[i] << " " << diff << std::endl;
        c_assert(abs(static_cast<int32_t>(naieve[i] - polymul[i])) <= 90000);
     }

    cout << "PolyMul Passed" << endl;

    return 0;
}