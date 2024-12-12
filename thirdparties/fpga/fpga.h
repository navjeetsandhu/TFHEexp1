#pragma once
#include "fftfpga.h"

void fpga_initialize();
void fpga_close();
fpga_t fpga_fft(const unsigned num, const float2 *inp, float2 *out, const bool inv, const unsigned batch=1);

