#pragma once

#include <array>
#include <cstdint>

#include "mulfft.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "decomposition.hpp"
#include "trgsw.hpp"
namespace TFHEpp {

template <class P>
TRGSWNTT<P> ApplyNTT2trgsw(const TRGSW<P> &trgsw)
{
    TRGSWNTT<P> trgswntt;
    for (int i = 0; i < (P::k + 1) * P::l; i++)
        for (int j = 0; j < P::k + 1; j++)
            TwistINTT<P>(trgswntt[i][j], trgsw[i][j]);
    return trgswntt;
}

template <class P>
TRGSWRAINTT<P> ApplyRAINTT2trgsw(const TRGSW<P> &trgsw)
{
    constexpr uint8_t remainder = ((P::nbit - 1) % 3) + 1;
    TRGSWRAINTT<P> trgswntt;
    for (int i = 0; i < (P::k + 1) * P::l; i++)
        for (int j = 0; j < P::k + 1; j++) {
            raintt::TwistINTT<typename P::T, P::nbit, true>(
                trgswntt[i][j], trgsw[i][j], (*raintttable)[1],
                (*raintttwist)[1]);
            for (int k = 0; k < P::n; k++)
                if ((k & ((1 << remainder) - 1)) > 1)
                    trgswntt[i][j][k] =
                        raintt::MulSREDC(trgswntt[i][j][k], raintt::R4);
                else
                    trgswntt[i][j][k] =
                        raintt::MulSREDC(trgswntt[i][j][k], raintt::R2);
        }
    return trgswntt;
}

template <class P>
TRGSWNTT<P> TRGSW2NTT(const TRGSW<P> &trgsw)
{
    TRGSWNTT<P> trgswntt;
    for (int i = 0; i < (P::k + 1) * P::l; i++)
        for (int j = 0; j < P::k + 1; j++) {
            PolynomialNTT<P> temp;
            TwistINTT<P>(temp, trgsw[i][j]);
            for (uint32_t k = 0; k < P::n; k++)
                trgswntt[i][j][k] = temp[cuHEpp::BitReverse<P::nbit>(k)];
        }
    return trgswntt;
}

template <class P>
TRGSWNTT<P> trgswnttSymEncrypt(const Polynomial<P> &p, const double alpha,
                               const Key<P> &key)
{
    TRGSW<P> trgsw = trgswSymEncrypt<P>(p, alpha, key);
    return ApplyNTT2trgsw<P>(trgsw);
}

template <class P>
TRGSWNTT<P> trgswnttSymEncrypt(const Polynomial<P> &p, const uint eta,
                               const Key<P> &key)
{
    TRGSW<P> trgsw = trgswSymEncrypt<P>(p, eta, key);
    return ApplyNTT2trgsw<P>(trgsw);
}

template <class P>
TRGSWNTT<P> trgswnttSymEncrypt(const Polynomial<P> &p, const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return trgswnttSymEncrypt<P>(p, P::alpha, key);
    else
        return trgswnttSymEncrypt<P>(p, P::eta, key);
}

template <class P>
TRGSWRAINTT<P> trgswrainttSymEncrypt(const Polynomial<P> &p, const double alpha,
                                     const Key<P> &key)
{
    TRGSW<P> trgsw = trgswSymEncrypt<P>(p, alpha, key);
    return ApplyRAINTT2trgsw<P>(trgsw);
}

template <class P>
TRGSWRAINTT<P> trgswrainttSymEncrypt(const Polynomial<P> &p, const uint eta,
                                     const Key<P> &key)
{
    if constexpr (hasq<P>::value && P::q == raintt::P) {
        constexpr uint8_t remainder = ((P::nbit - 1) % 3) + 1;
        constexpr std::array<typename P::T, P::l> h = hgen<P>();
        TRGSWRAINTT<P> trgswraintt;
        for (TRLWERAINTT<P> &trlweraintt : trgswraintt) {
            trlweraintt = trlwerainttSymEncryptZero<P>(eta, key);
            for (int k = 0; k <= P::k; k++)
                for (int j = 0; j < P::n; j++)
                    trlweraintt[k][j] = raintt::MulSREDC(
                        trlweraintt[k][j], ((j & ((1 << remainder) - 1)) > 1)
                                               ? raintt::R3
                                               : raintt::R2);
        }
        for (int i = 0; i < P::l; i++) {
            Polynomial<P> pscaled;
            for (int j = 0; j < P::n; j++) {
                pscaled[j] =
                    static_cast<uint64_t>(static_cast<int32_t>(p[j]) >= 0
                                              ? p[j]
                                              : (P::q + p[j])) *
                    h[i] % P::q;
            }
            for (int k = 0; k < P::k + 1; k++) {
                PolynomialRAINTT<P> praintt;
                raintt::TwistINTT<typename P::T, P::nbit, false>(
                    praintt, pscaled, (*raintttable)[1], (*raintttwist)[1]);
                for (int j = 0; j < P::n; j++)
                    if ((j & ((1 << remainder) - 1)) > 1)
                        praintt[j] = raintt::MulSREDC(praintt[j], raintt::R4);
                    else
                        praintt[j] = raintt::MulSREDC(praintt[j], raintt::R2);
                for (int j = 0; j < P::n; j++)
                    trgswraintt[i + k * P::l][k][j] = raintt::AddMod(
                        trgswraintt[i + k * P::l][k][j], praintt[j]);
            }
        }
        return trgswraintt;
    }
    else {
        TRGSW<P> trgsw = trgswSymEncrypt<P>(p, eta, key);
        return ApplyRAINTT2trgsw<P>(trgsw);
    }
}

template <class P>
TRGSWRAINTT<P> trgswrainttSymEncrypt(const Polynomial<P> &p, const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        return trgswrainttSymEncrypt<P>(p, P::alpha, key);
    else
        return trgswrainttSymEncrypt<P>(p, P::eta, key);
}

}  // namespace TFHEpp