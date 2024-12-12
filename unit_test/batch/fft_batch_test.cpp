#include "mult_fft.hpp"
#include "mult_fft_fpga.hpp"
#include <numeric>
#include <cmath>
#include "c_assert.hpp"
#include <chrono>

template <int nbits,  int batch>
void test_fft(uint32_t num_test)
{
    constexpr int N = 1 << nbits;
    constexpr unsigned size = batch *N;
    int i;
    std::array<uint32_t, size> result{};
    std::fill(result.begin(), result.end(), 0);

    std::array<uint32_t,  size> p1{};
    std::iota(p1.begin(), p1.end(), 1);

    auto start = std::chrono::high_resolution_clock::now();

    alignas(64) std::array<double, size> fft{};
    for (i=0; i< num_test; i++) TwistFpgaIFFTbatch(fft.data(), p1.data(), batch);

    auto finish1 = std::chrono::high_resolution_clock::now();

    for (i=0; i< num_test; i++) TwistFpgaFFTbatch(result.data(), fft.data(), batch);

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = finish - start;
    std::chrono::duration<double, std::milli> elapsed1 = finish1 - start;
    std::chrono::duration<double, std::milli> elapsed2 = finish - finish1;

    std::cout << "batch: " << batch << "\n";
    std::cout << "num iterations: " << num_test << std::endl;
    std::cout << "total run: " << batch * num_test << std::endl;
    std::cout << "Total Elapsed IFFT time: " << elapsed1.count() << " ms\n";
    std::cout << "Total Elapsed FFT time: " << elapsed2.count() << " ms\n";
    std::cout << "Total Elapsed time: " << elapsed.count() << " ms\n";
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
    constexpr int batch = 500;
    test_fft<nbit, batch>(num_test);
    return 0;
}