#pragma once
#include"../gate.hpp"

namespace TFHEpp{
#define INST(P) extern template void HomCONSTANTONE<P>(TLWE<P> & res)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P) extern template void HomCONSTANTZERO<P>(TLWE<P> & res)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P) extern template void HomNOT<P>(TLWE<P> & res, const TLWE<P> &ca)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P) extern template void HomCOPY<P>(TLWE<P> & res, const TLWE<P> &ca)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(iksP, brP, mu)                                                \
    extern template void HomNAND<iksP, brP, mu>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST

#define INST(iksP, brP, mu, batch)                                                \
    extern template void HomNANDbatch<iksP, brP, mu, batch>(TLWEn<typename brP::targetP, batch> &res, \
                                        const TLWEn<typename iksP::domainP, batch> &ca, \
                                        const TLWEn<typename iksP::domainP, batch> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_BATCH_IKSBR(INST)
#undef INST

#define INST(iksP, brP, mu)                                                \
    extern template void HomNOR<iksP, brP, mu>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST


#define INST(iksP, brP, mu)                                                \
    extern template void HomXNOR<iksP, brP, mu>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST


#define INST(iksP, brP, mu)                                                \
    extern template void HomAND<iksP, brP, mu>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST


#define INST(iksP, brP, mu)                                                \
    extern template void HomOR<iksP, brP, mu>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST


#define INST(iksP, brP, mu)                                                \
    extern template void HomXOR<iksP, brP, mu>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST


#define INST(iksP, brP, mu)                                                \
    extern template void HomANDNY<iksP, brP, mu>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST


#define INST(iksP, brP, mu)                                                \
    extern template void HomANDYN<iksP, brP, mu>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST


#define INST(iksP, brP, mu)                                                \
    extern template void HomORNY<iksP, brP, mu>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST


#define INST(iksP, brP, mu)                                                \
    extern template void HomORYN<iksP, brP, mu>(TLWE<typename brP::targetP> &res, \
                                        const TLWE<typename iksP::domainP> &ca, \
                                        const TLWE<typename iksP::domainP> &cb, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(INST)
#undef INST


#define INST(P)                                                   \
    extern template void HomMUX<P>(TLWE<P> & res, const TLWE<P> &cs,     \
                            const TLWE<P> &c1, const TLWE<P> &c0, \
                            const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(P)                                                    \
    extern template void HomNMUX<P>(TLWE<P> & res, const TLWE<P> &cs,     \
                             const TLWE<P> &c1, const TLWE<P> &c0, \
                             const EvalKey &ek)
INST(lvl1param);
INST(lvl0param);
#undef INST

#define INST(bkP)                                                              \
    extern template void HomMUXwoIKSandSE<bkP>(TRLWE<typename bkP::targetP> & res,    \
                                        const TLWE<typename bkP::domainP> &cs, \
                                        const TLWE<typename bkP::domainP> &c1, \
                                        const TLWE<typename bkP::domainP> &c0, \
                                        const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(INST)
#undef INST

#define INST(iksP, bkP)                         \
    extern template void HomMUXwoSE<iksP, bkP>(        \
        TRLWE<typename bkP::targetP> & res,     \
        const TLWE<typename iksP::domainP> &cs, \
        const TLWE<typename iksP::domainP> &c1, \
        const TLWE<typename iksP::domainP> &c0, const EvalKey &ek)
TFHEPP_EXPLICIT_INSTANTIATION_GATE(INST)
#undef INST

}