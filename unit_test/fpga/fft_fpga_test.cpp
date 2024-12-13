#include "mult_fft_fpga.hpp"
#include <numeric>
#include <cmath>
#include "c_assert.hpp"
#include <tfhe++.hpp>
#include <random>
#include <iostream>


using namespace std;

template <int nbits>
void test_fft(const std::array<uint32_t, 1 << nbits>& p1)
{

    std::string string_msg = "Input p1";
    constexpr int N = 1 << nbits;

    print_results<uint32_t>(string_msg, p1.data(), p1.size());

    std::array<uint32_t, N> result{};
    std::fill(result.begin(), result.end(), 0);

    alignas(64) std::array<double, N> fft{};
    TwistFpgaIFFT<N>(fft, p1);


    string_msg = "TwistIFFT 32 bit";
    print_results<double>(string_msg,  fft.data(), fft.size());

    TwistFpgaFFT<N>(result, fft);
    string_msg = "TwistFFT 32 bit";
    print_results<int32_t>(string_msg,  reinterpret_cast<int32_t*>(result.data()), result.size());


    cout <<"\n Difference between input and result \n"
    for (int i = 0; i < N; i++) {
        auto b = abs(static_cast<int32_t>(p1[i] - result[i]));
        cout << b << " ";
        //c_assert(b <= 1);
    }
    cout <<"\n";
}



template <class P, int nbits>
void test_fft_p()
{
    constexpr int N = 1 << nbits;

    std::array<P,  N> p1{};
    using namespace TFHEpp;
    //std::iota(p1.begin(), p1.end(), 1);
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> Torus32dist(0, UINT32_MAX);
    for (P &i : p1) i = Torus32dist(engine);
    test_fft<nbits>(p1);
}

int main()
{
    constexpr int nbit = 10;
    test_fft_p<uint32_t,nbit>();
    return 0;
}