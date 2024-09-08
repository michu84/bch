#pragma once

#include "polynomial.h"
#include "galois.h"
#include <cstring>

namespace mr {

    template<unsigned M, // power of 2 factor in the expressions for code length and number of elements in the Galois field
             unsigned ECC, // error correction capability
             unsigned...PrimPolyNonzeroPowers> // unit powers of primitive polynomial to generate Galois field
    struct bch {
        static_assert(M > 1, "sanity check!");

        constexpr static unsigned m = M;
        constexpr static unsigned t = ECC;
        constexpr static unsigned n = (1 << m) - 1; // code length in bits

        constexpr static unsigned n_bits_last_byte = n % 8;
        constexpr static bool has_trailing_bits = n_bits_last_byte != 0;
        constexpr static unsigned n_bytes = n / 8 + has_trailing_bits;

        using galois_field_type = binary_galois_field<m, PrimPolyNonzeroPowers...>;

        using galois_field_element_type = typename galois_field_type::element;
        using galois_field_element_polynomial_type = typename galois_field_type::element_polynomial_type;

        using minimal_polynomial_finder_type = polynomial<galois_field_element_type, m>; // each coefficient is an gf element polynomial
        using minimal_polynomial_type = typename galois_field_type::minimal_polynomial_type;

        using generator_polynomial_type = polynomial<bit_t, n>;

        using error_locator_polynomial_type = polynomial<galois_field_element_type, n-1>; // runtime use only, allow max theoretical error count equal to n just for sanity

        using error_pattern_polynomial_type = polynomial<bit_t, n-1>;

        // can't return constexpr std::vector, need to improvise with a struct, maybe std::array?
        using coset_power_type = unsigned int;

        struct cyclotomic_cosets_type {
            struct coset {
                coset_power_type base;
                coset_power_type power;
            };
            std::optional<coset> cosets[t][m]; // std::optional allows to fill in only valid fields with values, this may differ between bch<m, t, ...>, some odd powers may be missing because of cyclic redundancy
        };

        // can't return constexpr std::vector, need to improvise with a struct, maybe std::array?
        struct minimal_polynomials_type {
            struct coset_poly {
                coset_power_type coset_base;
                minimal_polynomial_type poly;
            };
            std::optional<coset_poly> polynomials[t]; // we must have t polynomials to correct t errors, but for syndromes we need 2*t with some repeating among available slots (cosets solve which apply where)
        };

        constexpr static cyclotomic_cosets_type make_cyclotomic_cosets() {
            cyclotomic_cosets_type cyclotomic;
            bool skip_redundant_mask[n] = {{0}}; // skip redundant alpha powers

            for(unsigned i=1, ci=0; i<2*t; i+=2) { // skip 0 as it's not effective anyway, step=2 since even powers still are ineffective
                auto &alpha_power_cosets = cyclotomic.cosets[ci];
                bool cosets_found = false;

                for(size_t j=0, cj=0; j<m; j++) {
                    const auto alpha_power = (i << j) % n;
                    if(!alpha_power // skip zero powers
                        || skip_redundant_mask[alpha_power]) { // don't repeat alpha_power values
                        continue;
                    }

                    alpha_power_cosets[cj++] = {i, alpha_power};
                    skip_redundant_mask[alpha_power] = true;
                    cosets_found |= true;
                }

                ci += cosets_found;
            }

            return cyclotomic;
        }

        constexpr static minimal_polynomials_type make_minimal_polynomials() {
            minimal_polynomials_type minimal_polys;

            const galois_field_element_polynomial_type coeff_1_poly(1, 0); // 1 as gf element polynomial
            const galois_field_element_type coeff_1_element(coeff_1_poly, 0); // 1 as gf element with gf element polynomial as coefficient
            const minimal_polynomial_finder_type finder_0;
            const minimal_polynomial_finder_type finder_1(coeff_1_element, 0); // 1*x^0 = 1*1 = 1, 1 as poly finder type
            const minimal_polynomial_finder_type finder_x(coeff_1_element, 1); // 1*x^1 = 1*x = x, x as poly finder type

            for(unsigned i=0, ci=0; i<t; i++) {
                const auto &coset_items = cyclotomic_cosets.cosets[i];

                const auto first_coset_item = coset_items[0];
                if(!first_coset_item.has_value())
                    break;

                const auto coset_base = first_coset_item->base;

                // get minimal poly

                auto this_minimal_poly_finder = finder_1;

                for(size_t j=0; j < m && coset_items[j].has_value(); j++) {
                    const auto &coset_item = coset_items[j];
                    const auto &alpha_element = galois_field_type::at(coset_item->power);

                    const minimal_polynomial_finder_type finder_alpha(alpha_element, 0);
                    assert(finder_alpha != finder_0);
                    assert(finder_x != finder_alpha);
                    this_minimal_poly_finder *= finder_x - finder_alpha;
                }

                const auto min_poly = galois_field_type::reduce_coefficients(this_minimal_poly_finder);

                minimal_polys.polynomials[ci++] = {coset_base, min_poly};
            }

            return minimal_polys;
        }

        constexpr static generator_polynomial_type make_generator_polynomial() {
            generator_polynomial_type gen_poly(1, 0); // represent scalar 1 to init the multiplications

            for(const auto &min_poly : minimal_polynomials.polynomials) {
                if(!min_poly.has_value())
                    break; // no more polys

                generator_polynomial_type min_as_gen_type;
                polynomial_copy(min_poly->poly, min_as_gen_type);

                gen_poly *= min_as_gen_type;
            }

            return gen_poly;
        }

        constexpr static size_t data_size_bits() {
            return n - generator_polynomial.degree();
        }

        constexpr static size_t parity_size_bits() {
            return generator_polynomial.degree();
        }

        constexpr static size_t word_size_bits() {
            return n;
        }

        constexpr static const cyclotomic_cosets_type cyclotomic_cosets = make_cyclotomic_cosets();
        constexpr static const minimal_polynomials_type minimal_polynomials = make_minimal_polynomials();
        constexpr static const generator_polynomial_type generator_polynomial = make_generator_polynomial();

        constexpr static const auto data_bits = data_size_bits();
        constexpr static const auto parity_bits = parity_size_bits();
        constexpr static const auto max_ecc_bits = m * t;

        static_assert(parity_bits <= max_ecc_bits, "effective number of configured parity bits exceed the maximum limit!");

        constexpr static unsigned n_data_bits_last_byte = data_bits % 8;
        constexpr static bool data_has_trailing_bits = n_data_bits_last_byte != 0;
        constexpr static unsigned n_data_bytes = data_bits / 8 + data_has_trailing_bits;

        struct syndrome_polynomials_gen_type {
            size_t global[2*t];
        };

        constexpr static std::optional<size_t> minimal_polynomial_global_index(const size_t &syndrome_index) {
            for(size_t c=0; c<t; c++) {
                const auto &cosets = cyclotomic_cosets.cosets[c];
                for(size_t s=0; (s < m) && cosets[s].has_value(); s++) {
                    const auto coset_power_match = cosets[s]->power == (syndrome_index+1);
                    const auto &poly = minimal_polynomials.polynomials[c];
                    const auto poly_coset_base_match = poly->coset_base == cosets[s]->base;
                    if(coset_power_match
                        && poly_coset_base_match)
                        return c;
                }
            }
            return {};
        }

        constexpr static syndrome_polynomials_gen_type make_syndrome_polynomials_gen() {
            syndrome_polynomials_gen_type syndrome_polys;
            for(unsigned i=0; i<2*t; i++) {
                const auto result = minimal_polynomial_global_index(i);
                if(result.has_value())
                    syndrome_polys.global[i] = result.value();
            }
            return syndrome_polys;
        }

        constexpr static syndrome_polynomials_gen_type syndrome_generator_index = make_syndrome_polynomials_gen();

        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        struct encoded_frame {
            constexpr static auto n_bytes = bch::n_bytes;

            uint8_t data_bytes[n_bytes]; // data + parity, scrambling mode doesn't affect size
        };

        static std::optional<encoded_frame> encode_codeword(const void *data_ptr, size_t offset_bits = 0) {

            // it processes bch_type::data_bits bits so make sure this amount is available after data

            const auto data_bytes = static_cast<const uint8_t *>(data_ptr);

            encoded_frame enc_frame;
            std::memset(enc_frame.data_bytes, 0, n_bytes);

            static_assert(generator_polynomial.degree() > 0, "");

            const auto _x = polynomial<bit_t, n>::make_from_memory(data_ptr, offset_bits);
            constexpr const polynomial<bit_t, n> _x_shift_degree(1, generator_polynomial.degree());
            const auto _input_shifted = (_x_shift_degree * _x).trimmed();
            const auto dd = _input_shifted / generator_polynomial;
            const auto encoded_parity_polynomial = dd.r; // ignoring division quotient, care only about residual

            // TODO: add proper schemes for configuring the encoded bits layout (scramble or not etc.)

    #ifdef PACK_UNSCRAMBLED_DATA_FIRST
            // MSB...LSB: [parityN, ..., parity0, dataN, ..., data0]
            // MSB...LSB: [byteM, ..., byte0]

            // fill data bytes (easier decoding if data first)

            const auto data_memcpy_bytes = data_bits / 8;
            std::memcpy(enc_frame.data_bytes, data_ptr, data_memcpy_bytes);

            for(unsigned i=0; i<data_bits; i++) {
                const auto global_byte_idx = i / 8;
                const auto local_bit_idx = i % 8;

                const bool data_bit_set = data_bytes[global_byte_idx] & (1U << local_bit_idx);
                enc_frame.data_bytes[global_byte_idx] |= (data_bit_set << local_bit_idx);
            }

            // fill parity bytes

            for(unsigned i=data_bits, j=0; i<n; i++, j++) {
                const auto global_enc_byte_idx = i / 8;
                const auto local_enc_bit_idx = i % 8;

                enc_frame.data_bytes[global_enc_byte_idx] |= (encoded_parity_polynomial[j] << local_enc_bit_idx);
            }
    #elif defined(PACK_UNSCRAMBLED_PARITY_FIRST)
            // MSB...LSB: [dataN, ..., data0, parityN, ..., parity0]
            // MSB...LSB: [byteM, ..., byte0]

            // fill parity bytes

            for(unsigned i=0; i<parity_bits; i++) {
                const auto global_enc_byte_idx = i / 8;
                const auto local_enc_bit_idx = i % 8;

                enc_frame.data_bytes[global_enc_byte_idx] |= (encoded_parity_polynomial[i] << local_enc_bit_idx);
            }

            // fill data bytes

            for(unsigned i=parity_bits, j=0; i<n; i++, j++) {
                const auto global_enc_byte_idx = i / 8;
                const auto local_enc_bit_idx = i % 8;

                const auto global_data_byte_idx = j / 8;
                const auto local_data_bit_idx = j % 8;

                const bool data_bit_set = data_bytes[global_data_byte_idx] & (1U << local_data_bit_idx);
                enc_frame.data_bytes[global_enc_byte_idx] |= (data_bit_set << local_enc_bit_idx);
            }

    #else // PACK_SCRAMBLED
            static_assert(false, "TODO");
    #endif
            return enc_frame;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        constexpr static uint8_t make_mask_last_byte() {
            uint8_t mask = 0;
            for(auto i=0; i<n_bits_last_byte; i++)
                mask |= (1 << i);
            return mask;
        }

        constexpr static const auto mask_last_byte = make_mask_last_byte();

        static void decode_output(const encoded_frame &encoded, void *out) {
            // copy all the full bytes right away

            // TODO: maybe act depending on packing strategy
            // packing_strategy::unpack(encoded, out);

            auto *data_bytes = static_cast<uint8_t*>(out);

            if(data_has_trailing_bits)
                std::memset(&data_bytes[n_data_bytes-1], 0, 1); // clear last byte so we don't have to clear trailing empty bits individually

    #ifdef PACK_UNSCRAMBLED_DATA_FIRST
            const auto memcpy_bytes = n_bytes - has_trailing_bits;
            std::memcpy(out, encoded.data_bytes, memcpy_bytes); // if there are no trailing bits, copy all in single pass, otherwise

            // finish with trailing bits only

            if(has_trailing_bits) {
                const auto *trailing_byte_in_ptr = encoded.data_bytes + memcpy_bytes;
                auto *trailing_byte_out_ptr = static_cast<uint8_t*>(out) + memcpy_bytes;
                for(size_t b=0; b<n_bits_last_byte; b++)
                    *trailing_byte_out_ptr |= *trailing_byte_in_ptr & (1 << b);
            }
    #elif defined(PACK_UNSCRAMBLED_PARITY_FIRST)
            for(unsigned i=parity_bits, j=0; i<n; i++, j++) {
                const auto global_enc_byte_idx = i / 8;
                const auto local_enc_bit_idx = i % 8;

                const auto global_data_byte_idx = j / 8;
                const auto local_data_bit_idx = j % 8;

                const bool data_bit_set = encoded.data_bytes[global_enc_byte_idx] & (1U << local_enc_bit_idx);
                data_bytes[global_data_byte_idx] |= (data_bit_set << local_data_bit_idx);
            }

    #else // PACK_SCRAMBLED
            static_assert(false):
    #endif
        }

        using syndrome_polynomial_type = polynomial<bit_t, m-1>;
        using syndrome_polynomials_type = std::array<syndrome_polynomial_type, 2*t>;
        using syndrome_elements_type = std::array<galois_field_element_type, 2*t>;

        static int decode_codeword(encoded_frame &encoded, void *out) {
            using codeword_poly_type = polynomial<bit_t, n-1>;
            const auto v = codeword_poly_type::make_from_memory(encoded.data_bytes);

    #ifdef DEBUG_VERBOSE
            printf("decoding poly:\t%s\n", v.to_string().c_str());
    #endif

            syndrome_polynomials_type syndrome_polys;
            syndrome_elements_type syndrome_elements;

            unsigned debug_errors_detected = 0;
            unsigned errors_corrected = 0;

            for (size_t c=0; c<2*t; c++) {
                const auto &idx = syndrome_generator_index.global[c];
                const auto &syndrome_gen_poly = minimal_polynomials.polynomials[idx];

                const auto tmp = v / syndrome_gen_poly->poly;

                assert(syndrome_polys[c].max_order >= tmp.r.degree());

                // unsafe copy trims tmp.r while copying to syndrome_polys[c]
                polynomial_copy_unsafe(tmp.r, syndrome_polys[c]); // ignoring the quotient, care about remainder only

                debug_errors_detected += !syndrome_polys[c].is_zero(); // nonzero syndrome indicate an error
            }

            // no errors only if all syndromes are zero

            if(debug_errors_detected) {

                // correct errors or return error if not recoverable

    #ifdef DEBUG_VERBOSE
                printf("decoding bit errors...\n");
    #endif

                // evaluate syndrome polynomials with respective gf elements

                for (size_t i=0; i<2*t; i++) {

                    // S[i](x) --> S[i](a^i) == x --> a^i, i=1...2*t

                    const auto &syndrome_poly = syndrome_polys[i];
                    auto &syndrome_element = syndrome_elements[i];

                    const auto &x = galois_field_type::at(i+1);

                    for(int j=0; j<signed(m); j++) {
                        if(syndrome_poly[j])
                            syndrome_element += x^j; // gf element power operator
                    }

    #ifdef DEBUG_VERBOSE
                    if(syndrome_poly.is_zero())
                        printf("syndrome[%li]: %s --> 0\n", i, syndrome_polys[i].to_string().c_str());
                    else
                        printf("syndrome[%li]: %s --> a^%3d\n", i, syndrome_polys[i].to_string().c_str(), *syndrome_element.exponent());
    #endif
                }

                // construct error locator polynomial

                // Berlekamp-Massey algorithm

                using elp_type = polynomial<galois_field_element_type, n-1>;

                size_t L = 0,
                       m = 1;
                galois_field_element_type d,
                                          b(galois_field_type::poly_by_power(0), 0); // 0
                elp_type C(galois_field_element_type(galois_field_type::poly_by_power(0), 0), 0),
                         B = C; // 1

                auto &S = syndrome_elements;
                const auto num_syndromes = 2*t;

                for(size_t k=0; k<num_syndromes; k++) {

                    // skip odd iterations for binary GF2 field

                    if ((k & 1) != 0) {
                        ++m;
                        continue;
                    }

                    // calculate the discrepancy

                    d = {};
                    for(size_t l=0; l<=L; l++)
                        d += C[l] * S[k-l];

                    if(d.poly().is_zero()) {

                        // we're good, just move on to tne next iteration

                        ++m;
                        continue;
                    }

                    // update C

                    const auto d_b = d * b.inverse();
                    const auto exponent = elp_type(d_b, m);
                    const auto delta = exponent * B;

                    if (2*L <= k){
                        B = C; // save for later

                        C -= delta.trimmed();
                        m = 1;

                        b = d; // save for later
                        L = k + 1 - L;

                    } else {
                        C -= delta.trimmed();
                        ++m;
                    }
                }
                const auto &error_locator_polynomial = C;

                const auto num_errors_found = error_locator_polynomial.degree();

                assert(num_errors_found);

                // roots of elp define the nonzero error location pattern coeffs

                codeword_poly_type errors_mask;

    #ifndef USE_BRUTE_FORCE_ELP_ROOTS_SEARCH

                // Chien's Search

                auto tmp_elp = error_locator_polynomial;

                for(size_t i=0; i<n; i++) {
                    galois_field_element_type sum;

                    for(size_t j=0; j<error_locator_polynomial.num_coeffs; j++)
                        sum += tmp_elp[j];

                    if(sum.poly().is_zero()) {
                        const auto root = galois_field_type::at(i).inverse();
                        const auto error_location = *root.exponent();

    #ifdef DEBUG_VERBOSE
                        printf("found error locator polynomial root @ a^%3ld --> error location: a^%3d\n", i, error_location);
    #endif

                        errors_mask[error_location] = true;
                        ++errors_corrected;
                    }

                    for(size_t j=0; i<n-1 && j<error_locator_polynomial.num_coeffs; j++)
                        tmp_elp[j] *= galois_field_type::at(j);
                }
    #else

                // brute force

                for(size_t i=0; i<n; i++) {
                    // see which elements are roots by substitution

                    galois_field_element_type root_candidate(galois_field_type::at(i));

                    // substitute into elp

                    galois_field_element_type result;

                    // int j because of galois_field_element_type::operator ^
                    for(int j=0; j<static_cast<const signed &>(error_locator_polynomial.num_coeffs); j++)
                        result += error_locator_polynomial[j] * (root_candidate ^ j);

                    if(result.poly().is_zero()) {
                        const auto root = root_candidate.inverse();
                        const auto error_location = *root.exponent();

    #ifdef DEBUG_VERBOSE
                        printf("found error locator polynomial root @ a^%3ld --> error location: a^%3d\n", i, error_location);
    #endif

                        errors_mask[error_location] = true;
                        ++errors_corrected;
                    }
                }
    #endif

                // correct error bits

    #ifdef DEBUG_VERBOSE
                printf("corrupted (%d errs):\t%s\nerror mask:\t\t%s\ncorrected:\t\t%s\n",
                    num_errors_found,
                    v.to_string().c_str(),
                    errors_mask.to_string().c_str(),
                    (v + errors_mask).to_string().c_str()
                );
    #endif

                for(size_t i=0; i<n; i++) {
                    if(errors_mask[i]) {
                        const auto global_byte_idx = i / 8;
                        const auto local_bit_idx = i % 8;
                        encoded.data_bytes[global_byte_idx] ^= (1U << local_bit_idx); // TODO: pass errors_mask to decode_output?
                    }
                }

                if(errors_corrected != num_errors_found)
                    return -num_errors_found;
            }

            decode_output(encoded, /*error_locator_polynomial, */out); // TODO: include error locator

            return errors_corrected; // no error or number of corrected errors
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////

        static void print_bch_info() {
            printf("BCH (m=%d, n=%d, t=ECC=%d, data_bits=%ld, parity_bits=%ld, primitive_polynomial=%s):\n",
                   m, n, t, data_size_bits(), parity_size_bits(), galois_field_type::primitive_polynomial.to_string().c_str());
        }

        static void print_galois_field() {
            galois_field_type::print();
        }

        static void print_primitive_polynomial() {
            printf("primitive polynomial {msb...lsb}: %s\n", galois_field_type::primitive_polynomial.to_string().c_str());
        }

        static void print_cyclotomic_cosets() {
            for(const auto &cosets : cyclotomic_cosets.cosets) {
                const auto &first_coset = cosets[0];
                if(!first_coset.has_value())
                    continue;

                printf("cyclotomic cosets [a^%3d]: ", first_coset->base);
                for(const auto &coset : cosets) {
                    if(!coset.has_value())
                        break;

                    printf("%u, ", coset->power);
                }
                printf("\n");
            }
        }

        static void print_minimal_polynomials() {
            for(const auto &poly : minimal_polynomials.polynomials) {
                if(!poly.has_value())
                    break;
                printf("minimal polynomial a^[%3d] {msb...lsb}: %s\n", poly->coset_base, poly->poly.to_string().c_str());
            }
        }

        static void print_generator_polynomial() {
            printf("generator polynomial (%dth order) {msb...lsb}: %s\n", generator_polynomial.degree(), generator_polynomial.to_string(true).c_str());
        }

        static void print_info() {
            print_bch_info();
            print_galois_field();
            print_primitive_polynomial();
            print_cyclotomic_cosets();
            print_minimal_polynomials();
            print_generator_polynomial();
        }
    };

}
