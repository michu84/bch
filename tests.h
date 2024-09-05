#pragma once

#include "bch.h"
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#ifdef DEBUG_VERBOSE_ENC_DEC
#include "measure_time.h"
#endif

void test_m5();
void test_m7();

void corrupt_bit_in_bytes(uint8_t *bytes, const auto corrupted_bit_idx) {
    const auto byte_idx = corrupted_bit_idx / 8;

#ifdef DEBUG_VERBOSE
    printf("testing by corrupting message bit %i @ byte %i\n", corrupted_bit_idx, byte_idx);
#endif

    bytes[byte_idx] ^= 1U << (corrupted_bit_idx % 8);
}

void corrupt_encoded_frame(auto &encoded, const auto add_errors, const auto msg_size_bits) {
    bool mask[msg_size_bits-1];
    unsigned target_bit_idx;
    for(unsigned i=0; i<add_errors; i++) {
        
        // do not repeat same bit multiple times
        
        do {
            target_bit_idx = rand() % msg_size_bits;
        } while(mask[target_bit_idx]); // try until we hit free bit

        mask[target_bit_idx] = true;

        corrupt_bit_in_bytes(encoded.data_bytes, target_bit_idx);        
    }

// #ifdef DEBUG_VERBOSE
//     using encoded_type = typename std::decay<decltype(encoded)>::type;
//     printf("test poly:\t\t%s\n", polynomial<bit_t, encoded_type::n_bytes>::make_from_memory(encoded.data_bytes).to_string().c_str());
// #endif
}

void corrupt_encoded_frame_example(auto &encoded) {
    corrupt_bit_in_bytes(encoded.data_bytes, 9);
    corrupt_bit_in_bytes(encoded.data_bytes, 22);
    corrupt_bit_in_bytes(encoded.data_bytes, 27);
}

void corrupt_encoded_frame_example_zero_syndrome_anomaly(auto &encoded) {
    corrupt_bit_in_bytes(encoded.data_bytes, 20);
    corrupt_bit_in_bytes(encoded.data_bytes, 21);
    corrupt_bit_in_bytes(encoded.data_bytes, 7);
}

template<typename T>
void print_as_hex(T &&x, size_t n, bool end_line = false) {
    printf("(hex:");

    for(size_t i=0; i<n; i++)
        printf(" %.2x", static_cast<const uint8_t &>(std::forward<T>(x)[i]));

    printf(")");

    if(end_line)
        printf("\n");
}

template<unsigned...>
struct print_poly_coeffs;

template<unsigned C, unsigned...Cs>
struct print_poly_coeffs<C, Cs...> {
    constexpr void operator() (std::stringstream &ss) const {
        print_poly_coeffs<C>{}
        (ss);

        if constexpr (sizeof...(Cs))
            ss << ", ";

        print_poly_coeffs<Cs...>{}
        (ss);
    }
};

template<unsigned C>
struct print_poly_coeffs<C> {
    constexpr void operator() (std::stringstream &ss) const {
        ss << C;
    }
};

template<unsigned m, unsigned t, unsigned...poly>
struct print_bch_type {
    constexpr void operator() (std::stringstream &ss) const {
        ss << "mr::bch<" << m << ", " << t << ", ";

        print_poly_coeffs<poly...>{}
        (ss);

        ss << ">";
    }
};

template<unsigned m, 
         unsigned t,
         unsigned...PrimitivePolynomialCoeffs>
int test_m_t_coeffs(unsigned add_errors) {
    using bch_type = mr::bch<m, t, PrimitivePolynomialCoeffs...>;

#ifdef DEBUG_VERBOSE
    bch_type::print_info();
#endif

    const std::string msg ="Hello\0";

    const auto msg_bytes = msg.c_str();
    const auto msg_size_bits = bch_type::data_bits + bch_type::parity_bits;
    const auto msg_size_bytes = msg_size_bits / 8 + (msg_size_bits % 8 != 0);

    std::stringstream bch_type_ss;
    print_bch_type<m, t, PrimitivePolynomialCoeffs...>{} (bch_type_ss);

    const auto bch_type_str = bch_type_ss.str();

    std::stringstream encode_type_ss;
    encode_type_ss << bch_type_str << "::encode_codeword...";

    auto encoded =
#ifdef DEBUG_VERBOSE_ENC_DEC
        mr::measure_time{} (encode_type_ss.str(), [&] {
            return bch_type::encode_codeword(msg_bytes);
        });
#else
        bch_type::encode_codeword(msg_bytes);
#endif

    if(add_errors > 0
        && encoded.has_value()) {

#ifdef DEBUG_VERBOSE
        printf("test poly:\t\t%s\n", mr::polynomial<bit_t, bch_type::n-1>::make_from_memory(encoded->data_bytes).to_string().c_str());
#endif

#if 1
        corrupt_encoded_frame(*encoded, add_errors, msg_size_bits);
#else
        corrupt_encoded_frame_example_zero_syndrome_anomaly(*encoded); // rare anomaly that broke the previous decoder implementation for m=5, t=3, p=0,2,5
        // corrupt_encoded_frame_example(*encoded); // just a working example
#endif

#if 0
#ifdef DEBUG_VERBOSE
        printf("test poly corrupted:   %s\n", polynomial<bit_t, bch_type::n-1>::make_from_memory(encoded->data_bytes).to_string().c_str());
#endif
#endif
    }

    // uint8_t corrupted_copy[bch_type::n_bytes];
    // std::copy(encoded->data_bytes, encoded->data_bytes + bch_type::n_bytes, corrupted_copy);
    std::vector<uint8_t> corrupted_copy(encoded->data_bytes, encoded->data_bytes + bch_type::n_bytes);
    
    char decoded[msg_size_bytes];
    std::memset(decoded, 0, msg_size_bytes);

    if(encoded.has_value()) {
        std::stringstream decode_type_ss;
        decode_type_ss << bch_type_str << "::decode_codeword...";

        const auto error_code =
#ifdef DEBUG_VERBOSE_ENC_DEC
            mr::measure_time{} (decode_type_ss.str(), [&] {
                return bch_type::decode_codeword(*encoded, &decoded);
            });
#else
            bch_type::decode_codeword(*encoded, &decoded);
#endif

        if(error_code < 0) {
#ifdef DEBUG_VERBOSE
            printf("detected %d errors\n", -error_code);
#endif
            assert(false);
        }

#ifdef DEBUG_VERBOSE
        printf("%s test result:\n", bch_type_str.c_str());

        printf("input:\t\t");
        print_as_hex(msg.c_str(), bch_type::n_data_bytes, false);
        printf(" \"%s\"\n", msg.c_str());

        printf("encoded:\t");
        print_as_hex(encoded->data_bytes, encoded->n_bytes, true);

        printf("corrupted:\t");
        print_as_hex(corrupted_copy, encoded->n_bytes, true);

        printf("decoded:\t");
        print_as_hex(decoded, bch_type::n_data_bytes, false);
        printf(" \"%s\"\n", decoded);
#endif

#if 1
        for(size_t bit=0; bit<bch_type::data_bits; bit++) {
            const auto byte_idx = bit / 8;
            const auto bit_idx = bit % 8;
            const auto lhs = decoded[byte_idx] & (1 << bit_idx);
            const auto rhs = msg_bytes[byte_idx] & (1 << bit_idx);

            if(lhs != rhs) {
                printf("%s != %s (FAIL)\n", decoded, msg_bytes);
                assert(false);
            }
        }
#endif
    }

    return 0;
}

template<unsigned m, 
         unsigned t,
         unsigned...PrimitivePolynomialCoeffs>
int test_m_t_coeffs_iterator(unsigned n_random_errors, unsigned n_times) {
    for(size_t j=n_random_errors; j>0; j--)
    for(size_t i=0; i<n_times; i++) {
        const auto result = test_m_t_coeffs<m, t, PrimitivePolynomialCoeffs...>(j);
        assert(result >= 0);
    }
    
    return 0;
}

template<unsigned m,
         unsigned t,
         unsigned...PrimitivePolynomialCoeffs>
int test_m_t_coeffs(unsigned n_random_errors, unsigned n_times) {
    for(size_t i=0; i<n_times; i++) {
        const auto result = test_m_t_coeffs<m, t, PrimitivePolynomialCoeffs...>(n_random_errors);
        assert(result >= 0);
    }

    return 0;
}

template<unsigned m, 
         unsigned t,
         unsigned...PrimitivePolynomialCoeffs,
         unsigned...S>
int test_m_t_coeffs_seq(mr::seq<S...>, unsigned n_times = 1) {
    (test_m_t_coeffs<m, t, PrimitivePolynomialCoeffs...>(S+1, n_times), ...);
    
    return 0;
}

void test_all();
