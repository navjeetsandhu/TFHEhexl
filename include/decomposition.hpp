#pragma once

#include <array>
#include <cstdint>

#include "mulfft.hpp"
#include "params.hpp"
#include "trlwe.hpp"


/*
 * https://www.zama.ai/post/tfhe-deep-dive-part-3
 * Decompose large coefficient into small base (Bg)
 */

namespace TFHEpp {


template <class P>
constexpr std::array<typename P::T, P::l> hgen()
{
    std::array<typename P::T, P::l> h{};
    if constexpr (hasq<P>::value)
        for (int i = 0; i < P::l; i++)
            h[i] = (P::q + (1ULL << ((i + 1) * P::Bgbit - 1))) >>
                   ((i + 1) * P::Bgbit);
    else
        for (int i = 0; i < P::l; i++)
            h[i] = 1ULL << (std::numeric_limits<typename P::T>::digits -
                            (i + 1) * P::Bgbit);
    return h;
}


// lvl1param  offset = 0x82080000,         2181562368
// lvl2param  offset = 0x8040201000000000, 9241421688455823360
template <class P>
constexpr typename P::T offsetgen()
{
    constexpr uint32_t  max_digits = std::numeric_limits<typename P::T>::digits;
    typename P::T offset = 0;
    // lvl1param   P::l = 3  P::Bgbit 6
    // lvl2param   P::l = 4  P::Bgbit 9
    for (int i = 1; i <= P::l; i++)
        offset += P::Bg / 2 * (1ULL << (max_digits - i * P::Bgbit));
    return offset;
}

template <class P>
constexpr typename P::T totalOffset()
{
    // lvl1 offset 0x82080000,         roundoffsetBit: 13, roundoffset: 0x2000,    totaloffset 0x82082000
    // lvl2 offset 0x8040201000000000, roundoffsetBit: 27, roundoffset: 0x8000000, totaloffset 0x8040201008000000
    constexpr typename P::T offset = offsetgen<P>();
    constexpr uint32_t roundoffsetBit = std::numeric_limits<typename P::T>::digits - P::l * P::Bgbit - 1;
    constexpr typename P::T roundoffset = 1ULL << roundoffsetBit;
    constexpr typename P::T totaloffset = offset + roundoffset;
    return  totaloffset;
}


template <class P>
inline void Decomposition(DecomposedPolynomial<P> &decpoly,
                          const Polynomial<P> &poly, typename P::T randbits = 0)
{

    constexpr typename P::T total_offset = totalOffset<P>();

    constexpr auto mask = static_cast<typename P::T>((1ULL << P::Bgbit) - 1);
    constexpr typename P::T halfBg = (1ULL << (P::Bgbit - 1));
    constexpr uint32_t maxDigits = std::numeric_limits<typename P::T>::digits;

    // lvl1   maxDigits 32, l 3  Bgbit 6, Bg 64,   halfBg: 32,    mask: 0x3f (6 digits)
    // lvl2   maxDigits 64, l 4  Bgbit 9, Bg 512,  halfBg: 256,   mask: 0x1ff (9 digits)

    for (int i = 0; i < P::n; i++) {
        auto valuePlusOffset = poly[i] + total_offset;
        for (int ii = 0; ii < P::l; ii++) {
            auto digitsToShift = maxDigits - (ii + 1) * P::Bgbit;
            auto shiftedValue = valuePlusOffset >> digitsToShift;
            auto maskedValue = shiftedValue & mask;
            decpoly[ii][i] = maskedValue - halfBg;
        }
    }
}


template <class P>
void DecompositionRAINTT(DecomposedPolynomialRAINTT<P> &decpolyntt,
                         const Polynomial<P> &poly, typename P::T randbits = 0)
{
    DecomposedPolynomial<P> decpoly;
    Decomposition<P>(decpoly, poly);
    for (int i = 0; i < P::l; i++)
        raintt::TwistINTT<typename P::T, P::nbit, false>(
            decpolyntt[i], decpoly[i], (*raintttable)[1], (*raintttwist)[1]);
}

}  // namespace TFHEpp