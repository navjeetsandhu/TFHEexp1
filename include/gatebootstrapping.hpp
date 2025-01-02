#pragma once

#include <cmath>
#include <limits>

#include "cloudkey.hpp"
#include "detwfa.hpp"
#include "keyswitch.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "utils.hpp"
#include "blindrotate.hpp"

namespace TFHEpp {

template <class P>
void GateBootstrappingTLWE2TLWEFFT(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector)
{
    alignas(64) TRLWE<typename P::targetP> acc;
    BlindRotate<P>(acc, tlwe, bkfft, testvector);
    SampleExtractIndex<typename P::targetP>(res, acc, 0);
}

template <class P, int batch>
void GateBootstrappingTLWE2TLWEFFTbatch(
    TLWEn<typename P::targetP, batch> &res, const TLWEn<typename P::domainP, batch> &tlwe,
    const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector)
{
    alignas(64) std::unique_ptr<TRLWEn<typename P::targetP, batch>> accPtr = std::make_unique<TRLWEn<typename P::targetP, batch>>();
    BlindRotatebatch<P, batch>(*accPtr, tlwe, bkfft, testvector);
    for (int j = 0; j < batch; j++)
        SampleExtractIndex<typename P::targetP>(res[j], (*accPtr)[j], 0);
}


template <class P, typename P::T mu>
constexpr Polynomial<P> mupolygen()
{
    Polynomial<P> poly;
    for (typename P::T &p : poly) p = mu;
    return poly;
}

template <class iksP, class bkP, typename bkP::targetP::T mu>
void GateBootstrapping(TLWE<typename iksP::domainP> &res,
                       const TLWE<typename iksP::domainP> &tlwe,
                       const EvalKey &ek)
{
    alignas(64) TLWE<typename iksP::targetP> tlwelvl0;
    IdentityKeySwitch<iksP>(tlwelvl0, tlwe, ek.getiksk<iksP>());
    GateBootstrappingTLWE2TLWEFFT<bkP>(res, tlwelvl0, ek.getbkfft<bkP>(),
                                       mupolygen<typename bkP::targetP, mu>());
}


template <class iksP, class bkP, typename bkP::targetP::T mu, int batch>
void GateBootstrappingbatch(TLWEn<typename iksP::domainP, batch> &res,
                       const TLWEn<typename iksP::domainP, batch> &tlwe,
                       const EvalKey &ek)
{
    alignas(64) std::unique_ptr<TLWEn<typename iksP::targetP, batch>> tlwelvl0Ptr = std::make_unique<TLWEn<typename iksP::targetP, batch>>();

    for (int j = 0; j < batch; j++)
        IdentityKeySwitch<iksP>((*tlwelvl0Ptr)[j], tlwe[j], ek.getiksk<iksP>());

    GateBootstrappingTLWE2TLWEFFTbatch<bkP, batch>(res, *tlwelvl0Ptr, ek.getbkfft<bkP>(),
                                       mupolygen<typename bkP::targetP, mu>());
}


}  // namespace TFHEpp