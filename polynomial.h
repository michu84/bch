#pragma once

#include <stddef.h>
#include <iostream>
#include <assert.h>
#include <string.h>
#include <cstdint>
#include <algorithm>

using bit_t = bool;

namespace mr {

    namespace {
        template<unsigned...>
        struct seq{};

        template<unsigned N, unsigned...Ids>
        struct gen_seq : gen_seq<N-1, N-1, Ids...> {};

        template<unsigned...Ids>
        struct gen_seq<0,Ids...> : seq<Ids...> {};

        template <class F, class... Args>
        constexpr void constexpr_for(F&& f, Args&&... args)
        {
            ( std::forward<F>(f) (std::forward<Args>(args)), ... );
        }

        template <class F, unsigned...Ids>
        constexpr void constexpr_for(F&& f, seq<Ids...>)
        {
            ( std::forward<F>(f) (Ids), ... );
        }

        constexpr size_t _log2(size_t n)
        {
            return ( (n<2) ? 1 : 1 + _log2(n/2));
        }

        constexpr auto poly_coeff_add(const auto &a, const auto &b)
        {
            return a + b;
        }

        constexpr auto poly_coeff_sub(const auto &a, const auto &b)
        {
            return a - b;
        }

        constexpr auto poly_coeff_mul(const auto &a, const auto &b)
        {
            return a * b;
        }

        constexpr bit_t poly_coeff_add(const bit_t &a, const bit_t &b)
        {
            return a ^ b; // modulo 2
        }

        constexpr bit_t poly_coeff_sub(const bit_t &a, const bit_t &b)
        {
            return poly_coeff_add(a, b); // same as addition, modulo 2
        }

        constexpr bit_t poly_coeff_mul(const bit_t &a, const bit_t &b)
        {
            return a & b; // modulo 2
        }

        template<unsigned...Is>
        struct select_last_unsigned
        {
            template<unsigned J, unsigned...Js>
            struct _ {
                constexpr static unsigned value = _<Js...>::value;
            };

            template<unsigned J>
            struct _<J> {
                constexpr static unsigned value = J;
            };

            constexpr static unsigned value = _<Is...>::value;
        };

    }

    template<typename C, unsigned max_order>
    struct polynomial;

    template<typename T, unsigned from_order, unsigned to_order>
    constexpr void polynomial_copy_unsafe(const polynomial<T,from_order> &from, polynomial<T,to_order> &to)
    {
        if constexpr (from_order < to_order) {
            constexpr auto copied = from_order + 1; // +1 for coeff x^0
            std::copy(from.coeffs, from.coeffs + copied, to.coeffs);
            std::fill_n(to.coeffs + copied, to_order - copied + 1, T{});
        } else {
            constexpr auto copied = to_order + 1; // +1 for coeff x^0
            std::copy(from.coeffs, from.coeffs + copied, to.coeffs);
        }

    }

    template<typename T, unsigned from_order, unsigned to_order>
    constexpr void polynomial_copy(const polynomial<T,from_order> &from, polynomial<T,to_order> &to)
    {
        static_assert(from_order <= to_order, "copy destination polynomial order not enough to contain all potential source polynomial data!");

        polynomial_copy_unsafe(from ,to);
    }

    template<typename T, unsigned from_order, unsigned to_order>
    constexpr void polynomial_move_unsafe(const polynomial<T,from_order> &from, polynomial<T,to_order> &to)
    {
        if constexpr (from_order < to_order) {
            constexpr auto moved = from_order + 1; // +1 for coeff x^0
            std::move(from.coeffs, from.coeffs + moved, to.coeffs);
            std::fill_n(to.coeffs + moved, to_order - moved + 1, T{});
        } else {
            constexpr auto moved = to_order + 1; // +1 for coeff x^0
            std::move(from.coeffs, from.coeffs + moved, to.coeffs);
        }
    }

    template<typename T, unsigned from_order, unsigned to_order>
    constexpr void polynomial_move(const polynomial<T,from_order> &from, polynomial<T,to_order> &to)
    {
        static_assert(from_order <= to_order, "copy destination polynomial order not enough to contain all potential source polynomial data!");

        polynomial_move_unsafe(from ,to);
    }

    // polynomial of type P(x) = C0*x^0 + C1*x^1 + C2*x^2 + ... + Cn*x^n

    template<typename C, unsigned MaxOrder>
    struct polynomial {

        static_assert(!std::is_same<C, unsigned>::value, "watch out!" );

        using coeff_type = C;

        constexpr static auto max_order = MaxOrder;
        constexpr static auto num_coeffs = max_order + 1;

        coeff_type coeffs[num_coeffs] = {};

        constexpr polynomial() {}

        constexpr polynomial(const coeff_type &v) {
            coeffs[0] = v;
        }

        constexpr polynomial(coeff_type &&v) {
            coeffs[0] = std::move(v);
        }

        template<typename...Args>
        constexpr polynomial(const coeff_type &coeff_x, const unsigned &x_power, Args &&...args)
            : polynomial(std::forward<Args>(args)...)
        {
            coeffs[x_power] = coeff_x;
        }

        template<unsigned other_max_order>
        constexpr polynomial(const polynomial<coeff_type, other_max_order> &other) {
            static_assert(other_max_order <= max_order, "detected potential polynomial overflow in copy ctor, a potential data loss!");
            // assert(other_max_order <= this->degree());
            polynomial_copy(other, *this);
        }

        template<unsigned other_max_order>
        constexpr polynomial(polynomial<coeff_type, other_max_order> &&other) {
            static_assert(other_max_order <= max_order, "detected potential polynomial overflow in copy ctor, a potential data loss!");
            // assert(other_max_order <= this->degree());
            polynomial_move(other, *this);
        }

        constexpr ~polynomial() {}

        constexpr polynomial(const polynomial &other)
        {
            polynomial_copy(other, *this);
        }

        constexpr polynomial(polynomial &&other)
        {
            polynomial_move(other, *this);
        }

        constexpr polynomial& operator = (const polynomial &other)
        {
            polynomial_copy(other, *this);
            return *this;
        }

        constexpr polynomial& operator = (polynomial &&other)
        {
            polynomial_move(other, *this);
            return *this;
        }

        constexpr static polynomial make_from_memory(const void *ptr, size_t offset_bits = 0, size_t start_power = 0 ) {
            polynomial poly;
            const auto data_bytes = static_cast<const uint8_t *>(ptr);

            for(unsigned i=offset_bits, j=0; i<offset_bits + num_coeffs; i++, j++) {
                const auto local_bit_idx = i % 8;
                const auto global_byte_idx = i / 8;
                const bool bit_is_set = data_bytes[global_byte_idx] & (1U << local_bit_idx);

                poly[start_power + j] = bit_is_set;

#if 0
                printf("polynomial, bit[%d,%d]: %s\n", global_byte_idx, local_bit_idx, bit_is_set ? "1" : "0");
#endif
            }

            return poly;
        }

        constexpr unsigned degree() const {
#if 0
            for(unsigned i=num_coeffs-1; i>=0; i--)
                if(static_cast<bool>(coeffs[i])) // either bool or has operator bool() const
                    return i;
            return 0;
#else
            unsigned d=0;
            for(unsigned i=0; i<num_coeffs; i++) {
                if(static_cast<bool>(coeffs[i]))
                    d = i;
            }
            return d;
#endif
        }

        constexpr bool is_zero() const {
            bool has_nonzero_coeff = false;
            for(const auto &coeff : coeffs) // don't use degree() to avoid cyclic dependency loop
                has_nonzero_coeff |= static_cast<bool>(coeff);
            return !has_nonzero_coeff;
        }

        constexpr operator bool() const {
            return !is_zero();
        }

        constexpr bool operator == (const polynomial &other) const {
            for(size_t i=0; i<num_coeffs; i++)
                if(coeffs[i] != other[i])
                    return false;

            return true;
        }

        constexpr const coeff_type& operator [] (const size_t &idx) const {
            return coeffs[idx];
        }

        constexpr coeff_type& operator [] (const size_t &idx) {
            return coeffs[idx];
        }

        constexpr polynomial operator + (const polynomial &other) const
        {
            // add coeffs

            polynomial out;
            const auto limit = std::max(degree(), other.degree());

            for(size_t i=0; i<limit+1; i++) {
                out[i] = poly_coeff_add( // if coeffs are bool, it evaluates as modulo2
                    coeffs[i],
                    other.coeffs[i]
                );
            }

            return out;
        }

        constexpr polynomial& operator += (const polynomial &other)
        {
            return (*this = (*this) + other);
        }

        constexpr polynomial operator - (const polynomial &other) const
        {
            // subtract coeffs

            polynomial out;
            const auto limit = std::max(degree(), other.degree());

            for(size_t i=0; i<limit+1; i++) {
                out[i] = poly_coeff_sub( // if coeffs are bool, it evaluates as modulo2
                    coeffs[i],
                    other.coeffs[i]
                );
            }

            return out;
        }

        constexpr polynomial& operator -= (const polynomial &other)
        {
            return (*this = (*this) - other);
        }

        // prevent overflow by doubling the precision in result
        // reducing afterwards brings result to the original max_order

        struct mul_result : polynomial<coeff_type, 2*max_order>
        {
            using trimmed_type = polynomial<coeff_type, max_order>;

            constexpr bool overflow() const {
                return this->degree() > max_order;
            }

            constexpr trimmed_type trimmed() const
            {
                trimmed_type out;

                polynomial_copy_unsafe(*this, out);

                return out;
            }
        };

        constexpr mul_result operator * (const polynomial &other) const
        {
            // convolve two polys coefficents

            mul_result out;

            for(size_t i=0; i<num_coeffs; i++)
                for(size_t j=0; j<other.num_coeffs; j++) {
                    const auto out_idx = i + j; // here we multiply
                    out[out_idx] = poly_coeff_add( // if coeffs are bool, it evaluates as modulo2
                        poly_coeff_mul(
                            coeffs[i],
                            other.coeffs[j]
                        ),
                        out[out_idx]
                    );
                }

            return out;
        }

        constexpr polynomial& operator *= (const polynomial &other)
        {
            return (*this = (((*this) * other).trimmed()));
        }

        struct div_result {
            polynomial q, r; // can safely assume max order same or lower than (*this) one?
        };

        constexpr div_result operator / (const polynomial &divisor) const
        {
            // static_assert(false, "figure out constexpr version of this function");
            assert(!divisor.is_zero());

            div_result out;

            // q = 0
            out.r = *this;

            while(!out.r.is_zero()
                  && (out.r.degree() >= divisor.degree()))
            {
                const auto powr_diff = out.r.degree() - divisor.degree();
                const polynomial dq(1, powr_diff);

                out.q += dq;
                out.r -= (dq * divisor).trimmed(); // TODO: is trimming needed? (optimization)
            }

            return out;
        }

        constexpr polynomial& operator /= (const polynomial &other)
        {
            return (*this = (*this) / other);
        }

        constexpr polynomial operator % (const polynomial &divisor) const
        {
            return ((*this) / divisor).r;
        }

        constexpr polynomial& operator %= (const polynomial &other)
        {
            return (*this = (*this) % other);
        }

        /*constexpr*/ std::string to_string(bool trim_leading_zeros = false) const
        {
            char bits[num_coeffs+1] = {};
            const auto s = trim_leading_zeros ? degree()+1 : num_coeffs;
            for(size_t i=0; i<s; i++)
                bits[i] = 48 + coeffs[s-1 - i];
            return std::string(bits);
        }

        /*constexpr*/ operator std::string () const
        {
            return to_string();
        }
    };

    template<typename T, unsigned max_order>
    std::ostream& operator << (std::ostream &s, const polynomial<T,max_order> &p)
    {
        for(size_t i=0; i<p.degree()+1; i++)
            s << p[p.num_coeffs-1 - i];
        return s;
    }

}
