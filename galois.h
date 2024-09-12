#pragma once

#include "gf2_polynomial.h"
#include <optional>
#include <array>
#include <assert.h>

namespace mr {

    template<unsigned M,
             unsigned...PrimPwrs>
    struct binary_galois_field {
        constexpr static const auto m = M;
        constexpr static const auto n = (1 << m) - 1;

        using init_field_polynomial_type = polynomial<bit_t, n>; // full range up to the nth order, used to represent dividend of modulo operation to get the residual for GF element polynomial representation
        using primitive_poly_type = polynomial<bit_t, m>; // only 2^m+1 divides it with 0 rest
        using minimal_polynomial_type = primitive_poly_type; // same as primitive poly

        using element_polynomial_type = polynomial<bit_t, m-1>; // no element has polynomial of m'th order so it's 1 less than primitive one
        using element_power_type = unsigned;

        constexpr static const auto primitive_polynomial = static_gf2_polynomial<PrimPwrs...>();

        // each element has two useful representations

        // 1. powers of alpha are just like indexes going from 0,1,2,...n. Useful in element multiplication

        // 2. polynomial representation is a unique polynomial in GF(2) for every alpha^i. Useful in element addition.

        constexpr static element_polynomial_type poly_by_power(const element_power_type &powr)
        {
            // TODO: check if powr is available
            assert(powr < n);
            return at(powr).poly();
        }

        constexpr static std::optional<element_power_type> power_by_poly(const element_polynomial_type &poly)
        {
            // may return no value for zero polynomial (result of arithmetic)

            // TODO: use lookup tables?

            for(const auto &element : elements)
                if(element.poly() == poly)
                    return element.exponent();
            return {};
        }

        struct element {
            std::optional<element_power_type> powr_rep; // aka simply the element alpha exponent, always going from 0...(2^m)-1
            element_polynomial_type  poly_rep = {}; // order of polynomials depends on generating the primitive_polynomial

            constexpr std::optional<element_power_type> exponent() const { return powr_rep; }
            constexpr element_polynomial_type poly()      const { return poly_rep; }

            constexpr element() {}

            constexpr element(const element_polynomial_type &poly_rep_, const std::optional<element_power_type> &powr_rep_)
                : powr_rep(powr_rep_)
                , poly_rep(poly_rep_) {}

            constexpr ~element() {}

            constexpr operator bool() const {
                return powr_rep.has_value()
                       && poly_rep != element_polynomial_type(0, 0);
            }

            constexpr bool operator == (const element &other) const
            {
                return powr_rep == other.powr_rep
                       && poly_rep == other.poly_rep;
            }

            constexpr element operator + (const element &other) const
            {
                const auto new_poly_rep = poly_rep + other.poly_rep; // simply add the polynomial coefficients, it'll automatically reduce itself with modulo2 arithmetic
                const auto new_powr_rep = binary_galois_field::power_by_poly(new_poly_rep);
                return element(new_poly_rep, new_powr_rep);
            }

            constexpr element& operator += (const element &other)
            {
                return (*this = *this + other);
            }

            constexpr element operator - (const element &other) const
            {
                const auto new_poly_rep = poly_rep - other.poly_rep; // simply add (it's not different from subtracting) the polynomial coefficients, it'll automatically reduce itself with modulo2 arithmetic
                const auto new_powr_rep = binary_galois_field::power_by_poly(new_poly_rep);
                return element(new_poly_rep, new_powr_rep);
            }

            constexpr element& operator -= (const element &other)
            {
                return (*this = *this - other);
            }

            constexpr element operator * (const element &other) const
            {
                if(!powr_rep.has_value()
                    || !other.powr_rep.has_value())
                    return binary_galois_field::elements[0];

                const auto new_powr_rep = (powr_rep.value() + other.powr_rep.value()) % n; // add the powers (shift the element) and modulo n to wrap around the cyclic stuff
                const auto new_poly_rep = binary_galois_field::poly_by_power(new_powr_rep);
                return element(new_poly_rep, new_powr_rep);
            }

            constexpr element& operator *= (const element &other)
            {
                return (*this = *this * other);
            }

            // taking the power

            constexpr element operator ^ (const int &_pow) const
            {
                if(!powr_rep.has_value())
                    return binary_galois_field::elements[0];

                if(_pow == 0)
                    return binary_galois_field::elements[1];

                const auto negative_pow = _pow < 0;

                const auto tmp_pow = negative_pow ? -_pow : _pow;

                const auto new_powr_rep = (powr_rep.value() * tmp_pow) % n; // multiply the powers (shift the element) and modulo n to wrap around the cyclic stuff
                const auto new_poly_rep = binary_galois_field::poly_by_power(new_powr_rep);
                const auto out = element(new_poly_rep, new_powr_rep);

                return negative_pow ? out.inverse() : out;
            }

            constexpr element& operator ^= (const int &pow)
            {
                return (*this = *this ^ pow);
            }

            constexpr element inverse() const {
                assert(exponent().has_value());

                const auto e = *exponent();
                const auto new_powr_rep = (n - e) % n;
                const auto new_poly_rep = binary_galois_field::poly_by_power(new_powr_rep);
                return element(new_poly_rep, new_powr_rep);
            }
        };

        using elements_type = std::array<element, n+2>; // +2 elements, one to represent 0, second to represent the cyclic wrap-around representation of 1 (with same value as element next to 0)

        constexpr static elements_type make_elements() {
            // extract the top coeff from primitive polynomial
            // and construct the replacement element for x^top = x^next + ...

            elements_type _elements;

            // trimming the primitive poly yields the required replacement polynomial

            element_polynomial_type gf_top_poly_replacement;
            polynomial_copy_unsafe(primitive_polynomial, gf_top_poly_replacement);

            const element_polynomial_type _0_poly;
            const element_polynomial_type _1_poly(1, 0);
            const element_polynomial_type _x_poly(1, 1);

            _elements[0] = element(_0_poly, {}); // include the root zero element
            _elements[1] = element(_1_poly, 0);

            for(auto i=2; i<n+2; i++) {
                const element_polynomial_type _prev = _elements[i-1].poly();
                const auto mul_result = _prev * _x_poly;

                const element_polynomial_type _next =
                    mul_result.trimmed()
                    + (mul_result.overflow() ? gf_top_poly_replacement : _0_poly);

                _elements[i] = element(_next, i-1); // actual element alpha power is i-1
            }

            return _elements;
        }

        constexpr static const elements_type elements = make_elements();

        constexpr static const element& at(const element_power_type &power_index) {
            return elements[power_index + 1]; // +1 for zero first element
        }

        constexpr const element& operator[] (const element_power_type & power_index) const {
            return at(power_index);
        }

        using minimal_polynomial_finder_type = polynomial<element, m>;

        constexpr static minimal_polynomial_type reduce_coefficients(const minimal_polynomial_finder_type &finder) {
            minimal_polynomial_type out;
            for(size_t i=0; i<out.num_coeffs; i++) {
                const auto has_value = finder.coeffs[i].exponent().has_value();
                const auto poly_is_zero = finder.is_zero();
                out.coeffs[i] = has_value
                                && !poly_is_zero;
            }
            return out;
        }

        static void print() {
            printf("Galois field (2^m), m=%d, n=2^m-1=%d:\n", m, n);

            for(auto i=0; i<=n; i++) {
                const auto element_poly = init_field_polynomial_type(1, i);

                printf("[a^%3d]: %s %% %s = %s%s\n",
                    i,
                    element_poly.to_string().c_str(),
                    primitive_polynomial.to_string().c_str(),
                    at(i).poly().to_string().c_str(),
                    i==n ? " (shown for clarity)" : ""
                );
            }
        }
    };

}
