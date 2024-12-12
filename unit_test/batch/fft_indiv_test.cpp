#include "mult_fft.hpp"
#include "mult_fft_fpga.hpp"
#include <numeric>
#include <cmath>
#include "c_assert.hpp"
#include <chrono>

template <int nbits>
void test_fft(const std::array<uint32_t, 1 << nbits>& p1, uint32_t num_test)
{
    constexpr int N = 1 << nbits;
    int i;
    auto start = std::chrono::high_resolution_clock::now();

    std::array<uint32_t, N> result{};
    std::fill(result.begin(), result.end(), 0);

    alignas(64) std::array<double, N> fft{};
    for (i=0; i< num_test; i++) TwistFpgaIFFT<N>(fft, p1);

    auto finish1 = std::chrono::high_resolution_clock::now();

    for (i=0; i< num_test; i++) TwistFpgaFFT<N>(result, fft);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = finish - start;
    std::chrono::duration<double, std::milli> elapsed1 = finish1 - start;
    std::chrono::duration<double, std::milli> elapsed2 = finish - finish1;

    std::cout << "Num tests: " << num_test << " s\n";
    std::cout << "Total Elapsed IFFT time: " << elapsed1.count() << " ms\n";
    std::cout << "Total Elapsed FFT time: " << elapsed2.count() << " ms\n";
    std::cout << "Total Elapsed time: " << elapsed.count() << " ms\n";
}



template <class P, int nbits>
void test_fft_p(uint32_t num_test)
{
    constexpr int N = 1 << nbits;

    std::array<P,  N> p1{};
    std::iota(p1.begin(), p1.end(), 1);
    test_fft<nbits>(p1, num_test);
}

int main(int argc, char* argv[])
{
    uint32_t num_test = 10;
    for(int i = 0; i < argc; i++) {
        std::cout << "Argument " << i << " is " << argv[i] << std::endl;
    }

    if(argc > 1) {
        std::stringstream str_stream(argv[1]);
        str_stream >> num_test;
    }

    std::cout << "num iterations: " << num_test << std::endl;
    constexpr int nbit = 10;
    test_fft_p<uint32_t,nbit>(num_test);
    return 0;
}