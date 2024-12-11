#pragma once

#include <array>
#include <cstdint>
#include "mulfft.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "decomposition.hpp"


namespace TFHEpp {

template <class P>
void trgswfftExternalProduct(TRLWE<P> &res, const TRLWE<P> &trlwe,
                             const TRGSWFFT<P> &trgswfft)
{
    alignas(64) DecomposedPolynomial<P> decpoly;
    Decomposition<P>(decpoly, trlwe[0]);
    alignas(64) PolynomialInFD<P> decpolyfft;
    TwistIFFT<P>(decpolyfft, decpoly[0]);
    alignas(64) TRLWEInFD<P> restrlwefft;
    for (int m = 0; m < P::k + 1; m++)
        MulInFD<P::n>(restrlwefft[m], decpolyfft, trgswfft[0][m]);
    for (int i = 1; i < P::l; i++) {
        TwistIFFT<P>(decpolyfft, decpoly[i]);
        for (int m = 0; m < P::k + 1; m++)
            FMAInFD<P::n>(restrlwefft[m], decpolyfft, trgswfft[i][m]);
    }
    for (int k = 1; k < P::k + 1; k++) {
        Decomposition<P>(decpoly, trlwe[k]);
        for (int i = 0; i < P::l; i++) {
            TwistIFFT<P>(decpolyfft, decpoly[i]);
            for (int m = 0; m < P::k + 1; m++)
                FMAInFD<P::n>(restrlwefft[m], decpolyfft,
                              trgswfft[i + k * P::l][m]);
        }
    }
    for (int k = 0; k < P::k + 1; k++) TwistFFT<P>(res[k], restrlwefft[k]);
}

}  // namespace TFHEpp