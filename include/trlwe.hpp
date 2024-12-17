#pragma once

#include "mulfft.hpp"
#include "params.hpp"

namespace TFHEpp {
template <class P>
TRLWE<P> trlweSymEncryptZero(const double alpha, const Key<P> &key)
{
    constexpr auto numeric_limit = std::numeric_limits<typename P::T>::max(); // i.e. 0xFFFFFFFF
    constexpr auto dimension = P::n; // i.e. 1024
    constexpr auto k_max = P::k;  // i.e 1
    std::uniform_int_distribution<typename P::T> Torusdist(0, numeric_limit);
    TRLWE<P> c;
    for (typename P::T &i : c[k_max]) i = ModularGaussian<P>(0, alpha);
    for (int k = 0; k < k_max; k++) {
        for (typename P::T &i : c[k]) i = Torusdist(generator);
        std::array<typename P::T, dimension> partkey;
        for (int i = 0; i < dimension; i++) partkey[i] = key[k * dimension + i];
        Polynomial<P> temp;
        PolyMul<P>(temp, c[k], partkey);
        for (int i = 0; i < dimension; i++) c[k_max][i] += temp[i];
    }
    return c;
}

template <class P, int batch>
TRLWEn<P, batch> trlweSymEncryptZerobatch(const double alpha, const Key<P> &key)
{
    constexpr auto numeric_limit = std::numeric_limits<typename P::T>::max(); // i.e. 0xFFFFFFFF
    constexpr auto dimension = P::n; // i.e. 1024
    constexpr auto k_max = P::k;  // i.e 1
    std::uniform_int_distribution<typename P::T> Torusdist(0, numeric_limit);
    TRLWEn<P, batch> c;
    for (int j=0;j<batch;j++)
        for (typename P::T &i : c[k_max][j]) i = ModularGaussian<P>(0, alpha);

    for (int k = 0; k < k_max; k++) {
        for (int j=0;j<batch;j++)
            for (typename P::T &i : c[k][j]) i = Torusdist(generator);

        std::array<typename P::T, dimension> partkey;
        for (int i = 0; i < dimension; i++) partkey[i] = key[k * dimension + i];

        Polynomialn<P, batch> temp;
        for (int j=0;j<batch;j++) {
            PolyMul<P>(temp[j], c[k], partkey);
            for (int i = 0; i < dimension; i++) c[k_max][j][i] += temp[j][i];
        }
    }
    return c;
}

template <class P>
TRLWE<P> trlweSymEncryptZero(const uint eta, const Key<P> &key)
{
    constexpr auto numeric_limit = std::numeric_limits<typename P::T>::max(); // i.e. 0xFFFFFFFF
    constexpr auto dimension = P::n; // i.e. 1024
    constexpr auto k_max = P::k; // i.e 1
    std::uniform_int_distribution<typename P::T> Torusdist(0, numeric_limit);

    alignas(64) TRLWE<P> c;
    for (typename P::T &i : c[k_max])
        i = (CenteredBinomial<P>(eta) << std::numeric_limits<P>::digits) / P::q;


    for (int k = 0; k < k_max; k++) {
        for (typename P::T &i : c[k]) i = Torusdist(generator);
        alignas(64) std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < dimension; i++) partkey[i] = key[k * dimension + i];
        alignas(64) Polynomial<P> temp;
        PolyMul<P>(temp, c[k], partkey);
        for (int i = 0; i < dimension; i++) c[k_max][i] += temp[i];
    }
    return c;
}

template <class P, int batch>
TRLWEn<P, batch> trlweSymEncryptZerobatch(const uint eta, const Key<P> &key)
{
    constexpr auto numeric_limit = std::numeric_limits<typename P::T>::max(); // i.e. 0xFFFFFFFF
    constexpr auto dimension = P::n; // i.e. 1024
    constexpr auto k_max = P::k; // i.e 1
    std::uniform_int_distribution<typename P::T> Torusdist(0, numeric_limit);
    alignas(64) TRLWEn<P, batch> c;

    for (int j=0;j<batch;j++)
        for (typename P::T &i : c[k_max][j])
            i = (CenteredBinomial<P>(eta) << std::numeric_limits<P>::digits) / P::q;

    for (int k = 0; k < k_max; k++) {
        for (int j=0;j<batch;j++)
            for (typename P::T &i : c[k][j]) i = Torusdist(generator);

        alignas(64) std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < dimension; i++) partkey[i] = key[k * dimension + i];

        alignas(64) Polynomialn<P, batch> temp;
        for (int j=0;j<batch;j++) {
            PolyMul<P>(temp[j], c[k], partkey);
            for (int i = 0; i < dimension; i++) c[k_max][j][i] += temp[j][i];
        }
    }
    return c;
}

template <class P>
TRLWE<P> trlweSymEncryptZero(const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return trlweSymEncryptZero<P>(P::alpha, key);
    else
        return trlweSymEncryptZero<P>(P::eta, key);
}

template <class P, int batch>
TRLWEn<P,batch> trlweSymEncryptZerobatch(const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return trlweSymEncryptZerobatch<P,batch>(P::alpha, key);
    else
        return trlweSymEncryptZerobatch<P,batch>(P::eta, key);
}


template <class P>
TRLWE<P> trlweSymEncrypt(const std::array<typename P::T, P::n> &p,
                         const double alpha, const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(alpha, key);
    for (int i = 0; i < P::n; i++) c[P::k][i] += p[i];
    return c;
}


template <class P>
TRLWE<P> trlweSymEncrypt(const std::array<typename P::T, P::n> &p, const uint eta,
                         const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(eta, key);
    for (int i = 0; i < P::n; i++) c[P::k][i] += p[i];
    return c;
}

template <class P>
TRLWE<P> trlweSymEncrypt(const std::array<typename P::T, P::n> &p,
                         const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return trlweSymEncrypt<P>(p, P::alpha, key);
    else
        return trlweSymEncrypt<P>(p, P::eta, key);
}


template <class P>
TRLWE<P> trlweSymIntEncrypt(const std::array<typename P::T, P::n> &p,
                            const double alpha, const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(alpha, key);
    for (int i = 0; i < P::n; i++)
        c[P::k][i] += static_cast<typename P::T>(P::delta * p[i]);
    return c;
}

template <class P>
TRLWE<P> trlweSymIntEncrypt(const std::array<typename P::T, P::n> &p,
                            const uint eta, const Key<P> &key)
{
    TRLWE<P> c = trlweSymEncryptZero<P>(eta, key);
    for (int i = 0; i < P::n; i++)
        c[P::k][i] += static_cast<typename P::T>(P::delta * p[i]);
    return c;
}

template <class P>
TRLWE<P> trlweSymIntEncrypt(const std::array<typename P::T, P::n> &p,
                            const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return trlweSymIntEncrypt<P>(p, P::alpha, key);
    else
        return trlweSymIntEncrypt<P>(p, P::eta, key);
}

template <class P>
Polynomial<P> trlwePhase(const TRLWE<P> &c, const Key<P> &key)
{
    Polynomial<P> phase = c[P::k];
    for (int k = 0; k < P::k; k++) {
        Polynomial<P> mulres;
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        PolyMul<P>(mulres, c[k], partkey);
        for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];
    }
    return phase;
}

template <class P>
std::array<bool, P::n> trlweSymDecrypt(const TRLWE<P> &c, const Key<P> &key)
{
    Polynomial<P> phase = trlwePhase<P>(c, key);

    std::array<bool, P::n> p;
    if constexpr (hasq<P>::value) {
        for (int i = 0; i < P::n; i++) p[i] = (phase[i] % P::q) < P::q / 2;
    }
    else
        for (int i = 0; i < P::n; i++)
            p[i] = static_cast<typename std::make_signed<typename P::T>::type>(
                       phase[i]) > 0;
    return p;
}

template <class P>
Polynomial<P> trlweSymIntDecrypt(const TRLWE<P> &c, const Key<P> &key)
{
    Polynomial<P> phase = c[P::k];
    for (int k = 0; k < P::k; k++) {
        Polynomial<P> mulres;
        std::array<typename P::T, P::n> partkey;
        for (int i = 0; i < P::n; i++) partkey[i] = key[k * P::n + i];
        PolyMul<P>(mulres, c[k], partkey);
        for (int i = 0; i < P::n; i++) phase[i] -= mulres[i];
    }

    Polynomial<P> p;
    for (int i = 0; i < P::n; i++)
        p[i] = static_cast<typename P::T>(std::round(phase[i] / P::delta)) %
               P::plain_modulus;
    return p;
}

template <class P>
void SampleExtractIndex(TLWE<P> &tlwe, const TRLWE<P> &trlwe, const int index)
{
    for (int k = 0; k < P::k; k++) {
        for (int i = 0; i <= index; i++)
            tlwe[k * P::n + i] = trlwe[k][index - i];
        for (int i = index + 1; i < P::n; i++)
            tlwe[k * P::n + i] = -trlwe[k][P::n + index - i];
    }
    tlwe[P::k * P::n] = trlwe[P::k][index];
}

template <class P>
void InvSampleExtractIndex(TRLWE<P> &trlwe, const TLWE<P> &tlwe,
                           const int index)
{
    for (int k = 0; k < P::k; k++) {
        for (int i = 0; i <= index; i++)
            trlwe[k][index - i] = tlwe[k * P::n + i];
        for (int i = index + 1; i < P::n; i++)
            trlwe[k][P::n + index - i] = -tlwe[k * P::n + i];
    }
    trlwe[P::k] = {};
    trlwe[P::k][index] = tlwe[P::k * P::n];
}

}  // namespace TFHEpp