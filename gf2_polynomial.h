#pragma once

#include "polynomial.h"
#include <utility>

namespace mr {

    template<unsigned max_order>
    struct gf2_polynomial : polynomial<bit_t, max_order>
    {
        using base = polynomial<bit_t, max_order>;

        template<typename...Powers>
        constexpr gf2_polynomial(const unsigned &power, Powers &&...powers)
            : gf2_polynomial::gf2_polynomial(std::forward<Powers>(powers)...)
        {
            base::coeffs[power] = 1;
        }

    private:
        using base::base; // hide the inherited ctor
    };

    template<unsigned...Pwrs>
    struct static_gf2_polynomial : gf2_polynomial<select_last_unsigned<Pwrs...>::value>
    {
        constexpr static const unsigned max_order = select_last_unsigned<Pwrs...>::value;
        using base = gf2_polynomial<max_order>;

        constexpr static_gf2_polynomial()
            : base(Pwrs...) {} // represented value fixed by template args

    private:
        using base::base; // hide the inherited ctor
    };

// detect 128 bit builtin types
#ifdef __SIZEOF_INT128__
    using int128_t  = __int128_t;
    using uint128_t = __uint128_t;

    // specialize polynomials with coeffs in GF(2) with casting as uint128_t

    template<unsigned max_order>
    constexpr void gf2_polynomial_to_decimal(const polynomial<bit_t, max_order> &p, uint128_t &d)
    {
        const auto deg = p.degree();
        d = 0;
        for(auto i=0; i<=deg; i++)
            d |= (p[i] << i);
    }

    template<unsigned max_order>
    constexpr void decimal_to_gf2_polynomial(const uint128_t &d, polynomial<bit_t, max_order> &p)
    {
        // degree not known without clz
        for(auto i=0; i<max_order; i++) {
            const auto dbg = d & (1U << i);
            p[i] = dbg;
        }
    }
#endif

}
