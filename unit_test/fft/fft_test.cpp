#include "mult_fft_fpga.hpp"
#include "c_assert.hpp"
#include <numeric>
#include <cmath>

constexpr int batch = 2;

template <int nbits>
void test_fft(const std::array<uint32_t, batch * (1 << nbits)>& p1)
{

    std::string string_msg = "Input p1";
    constexpr int N = 1 << nbits;

    print_results<uint32_t>(string_msg, p1.data(), p1.size());

    std::array<uint32_t, batch * N> result{};
    std::fill(result.begin(), result.end(), 0);

    alignas(64) std::array<double, batch *N> fft{};
    TwistFpgaIFFTbatch(fft.data(), p1.data(), batch);

    string_msg = "TwistIFFT 32 bit";
    print_results<double>(string_msg,  fft.data(), fft.size());

    TwistFpgaFFTbatch(result.data(), fft.data(), batch);
    string_msg = "TwistFFT 32 bit";
    print_results<int32_t>(string_msg,  reinterpret_cast<int32_t*>(result.data()), result.size());

}


template <class P, int nbits>
void test_fft_p()
{
    constexpr int N = 1 << nbits;

    std::array<P, batch*N> p1{};
    std::iota(p1.begin(), p1.end(), 1);
    test_fft<nbits>(p1);
}

int main()
{
    constexpr int nbit = 10;
    test_fft_p<uint32_t,nbit>();
    return 0;
}