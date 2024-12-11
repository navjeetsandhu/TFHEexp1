#include "decomposition.hpp"
#include "trgsw.hpp"
#include "externalproduct.hpp"

namespace TFHEpp {
#define INST(P)                                                       \
    template void Decomposition<P>(DecomposedPolynomial<P> & decpoly, \
                                   const Polynomial<P> &poly,         \
                                   typename P::T randbits)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

#define INST(P)                               \
    template void trgswfftExternalProduct<P>( \
        TRLWE<P> & res, const TRLWE<P> &trlwe, const TRGSWFFT<P> &trgswfft)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST


#define INST(P)                                                  \
    template TRGSW<P> trgswSymEncrypt<P>(const Polynomial<P> &p, \
                                         const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

#define INST(P)                                                        \
    template TRGSWFFT<P> trgswfftSymEncrypt<P>(const Polynomial<P> &p, \
                                               const Key<P> &key)
TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(INST)
#undef INST

}  // namespace TFHEpp
