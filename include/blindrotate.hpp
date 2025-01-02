#pragma once

#include <cmath>
#include <limits>
#include "cloudkey.hpp"
#include "detwfa.hpp"
#include "keyswitch.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "utils.hpp"

namespace TFHEpp {

template <class P, uint32_t num_out = 1>
void BlindRotate(TRLWE<typename P::targetP> &res,
                 const TLWE<typename P::domainP> &tlwe,
                 const BootstrappingKeyFFT<P> &bkfft,
                 const Polynomial<typename P::targetP> &testvector)
{
    //cout << "b " ;
    constexpr uint32_t bitwidth = bits_needed<num_out - 1>();
    const uint32_t bLong = 2 * P::targetP::n -
                       ((tlwe[P::domainP::k * P::domainP::n] >>
                         (std::numeric_limits<typename P::domainP::T>::digits -
                          1 - P::targetP::nbit + bitwidth))
                        << bitwidth);
    res = {};
    PolynomialMulByXai<typename P::targetP>(res[P::targetP::k], testvector, bLong);

    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        constexpr typename P::domainP::T roundoffset =
            1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                     P::targetP::nbit + bitwidth);
        const uint32_t aLong =
            (tlwe[i] + roundoffset) >>
            (std::numeric_limits<typename P::domainP::T>::digits - 1 -
             P::targetP::nbit + bitwidth)
                << bitwidth;
        //cout << " along " << aLong << " ";
        if (aLong == 0) continue;
        // Do not use CMUXFFT to avoid unnecessary copy.
        CMUXFFTwithPolynomialMulByXaiMinusOne<P>(res, bkfft[i], aLong);
    }
}


template <class P, int batch, uint32_t num_out = 1>
void BlindRotatebatch(TRLWEn<typename P::targetP, batch> &res,
                 const TLWEn<typename P::domainP, batch> &tlwe,
                 const BootstrappingKeyFFT<P> &bkfft,
                 const Polynomial<typename P::targetP> &testvector)
{
    res = {};
    constexpr uint32_t bitwidth = bits_needed<num_out - 1>();
    for (int j = 0; j < batch; j++) {
        const uint32_t bLong =
            2 * P::targetP::n -
            ((tlwe[j][P::domainP::k * P::domainP::n] >>
              (std::numeric_limits<typename P::domainP::T>::digits - 1 -
               P::targetP::nbit + bitwidth))
             << bitwidth);

        PolynomialMulByXai<typename P::targetP>(res[j][P::targetP::k], testvector,
                                                bLong);
    }
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        intArray<batch> aLongArray;
        for (int j = 0; j < batch; j++) {
            constexpr typename P::domainP::T roundoffset =
                1ULL << (std::numeric_limits<typename P::domainP::T>::digits -
                         2 - P::targetP::nbit + bitwidth);
            aLongArray[j] =
                (tlwe[j][i] + roundoffset) >>
                (std::numeric_limits<typename P::domainP::T>::digits - 1 -
                 P::targetP::nbit + bitwidth)
                    << bitwidth;
        }
        // Do not use CMUXFFT to avoid unnecessary copy.
        CMUXFFTwithPolynomialMulByXaiMinusOnebatch<P, batch>(res, bkfft[i], aLongArray);
    }
}

}  // namespace TFHEpp