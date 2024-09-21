#pragma once

#include "bch.h"
#include "measure_time.h"

#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

void test_m5();
void test_m7();

void corrupt_bit_in_bytes(uint8_t *bytes, const auto corrupted_bit_idx) {
    const auto byte_idx = corrupted_bit_idx / 8;

#ifdef DEBUG_VERBOSE
    std::cout << "testing by corrupting message bit " << corrupted_bit_idx
              << " @ byte " << byte_idx << std::endl;
#endif

    bytes[byte_idx] ^= 1U << (corrupted_bit_idx % 8);
}

void corrupt_encoded_frame(auto &encoded, const auto add_errors, const auto msg_size_bits) {
    bool mask[msg_size_bits] = {};

    unsigned target_bit_idx;
    for(unsigned i=0; i<add_errors; i++) {

        // do not repeat same bit multiple times

        do {
            target_bit_idx = rand() % msg_size_bits;
        } while(mask[target_bit_idx]); // try until we hit free bit

        mask[target_bit_idx] = true;

        corrupt_bit_in_bytes(encoded.data_bytes, target_bit_idx);
    }
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
    std::cout << "(hex:" << std::hex << std::setfill('0');

    for(size_t i=0; i<n; i++)
        std::cout << " " << std::setw(2) << static_cast<int>(std::forward<T>(x)[i]);

    std::cout << std::setw(1) << std::dec << std::setfill(' ') << ")";

    if(end_line)
        std::cout << std::endl;
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

template<typename>
struct print_bch;

template<unsigned m, unsigned t, unsigned...poly>
struct print_bch<mr::bch<m, t, poly...>> {
    constexpr void operator() (std::stringstream &ss) const {
        ss << "mr::bch<" << m << ", " << t << ", ";

        print_poly_coeffs<poly...>{}
        (ss);

        ss << ">";
    }
};

template<typename>
struct test_bch;

template<unsigned m,
         unsigned t,
         unsigned...PrimitivePolynomialCoeffs>
struct test_bch<mr::bch<m, t, PrimitivePolynomialCoeffs...>> {
    using bch_type = mr::bch<m, t, PrimitivePolynomialCoeffs...>;

    void operator() (unsigned add_errors) const {

#ifdef DEBUG_VERBOSE
        bch_type::print_info();
#endif

        // initialize a full length buffer so codec won't access unrelated memory if data size less than full capacity

        char msg_buffer[bch_type::n_bytes] = {};
        std::snprintf(msg_buffer, 5+1,  "Hello");
        const std::string msg(msg_buffer, bch_type::n_data_bytes);

        const auto msg_bytes = msg.c_str();
        constexpr auto msg_size_bits = bch_type::data_bits + bch_type::parity_bits;
        constexpr auto msg_size_bytes = msg_size_bits / 8 + (msg_size_bits % 8 != 0);

        std::stringstream bch_type_ss;
        print_bch<bch_type>{} (bch_type_ss);

        const auto bch_type_str = bch_type_ss.str();

        std::stringstream encode_type_ss;
        encode_type_ss << bch_type_str << "::encode_codeword...";

#ifdef DEBUG_VERBOSE_ENC_DEC
        const auto measure_result = mr::measure_time{} (encode_type_ss.str(), [&] {
            return bch_type::encode_codeword(msg_bytes);
        });
        auto encoded = measure_result.result;
#else
        auto encoded = bch_type::encode_codeword(msg_bytes);
#endif

        if(add_errors > 0
            && encoded.has_value()) {

#ifdef DEBUG_VERBOSE
            std::cout << "test poly:\t\t" << mr::polynomial<bit_t, bch_type::n-1>::make_from_memory(encoded->data_bytes).to_string() << std::endl;
#endif

#if 1
            corrupt_encoded_frame(*encoded, add_errors, msg_size_bits);
#else
            corrupt_encoded_frame_example_zero_syndrome_anomaly(*encoded); // rare anomaly that broke the previous decoder implementation for m=5, t=3, p=0,2,5
            // corrupt_encoded_frame_example(*encoded); // just a working example
#endif
        }

        const std::vector<uint8_t> corrupted_copy(encoded->data_bytes, encoded->data_bytes + bch_type::n_bytes);

        char decoded[msg_size_bytes] = {};

        if(encoded.has_value()) {
            std::stringstream decode_type_ss;
            decode_type_ss << bch_type_str << "::decode_codeword...";


#ifdef DEBUG_VERBOSE_ENC_DEC
            const auto measure_result = mr::measure_time{} (decode_type_ss.str(), [&] {
                return bch_type::decode_codeword(*encoded, &decoded);
            });
            const auto error_code = measure_result.result;
#else
            const auto error_code = bch_type::decode_codeword(*encoded, &decoded);
#endif

            if(error_code < 0) {
#ifdef DEBUG_VERBOSE
                std::cout << "detected " << -error_code << " errors" << std::endl;
#endif
                assert(false);
            }

#ifdef DEBUG_VERBOSE
            std::cout << bch_type_str << " test results:" << std::endl;

            std::cout << "input:\t\t";
            print_as_hex(msg.c_str(), bch_type::n_data_bytes, false);
            std::cout << " \"" << msg << "\"" << std::endl;

            std::cout << "encoded:\t";
            print_as_hex(encoded->data_bytes, encoded->n_bytes, true);

            std::cout << "corrupted:\t";
            print_as_hex(corrupted_copy, encoded->n_bytes, true);

            std::cout << "decoded:\t";
            print_as_hex(decoded, bch_type::n_data_bytes, false);
            std::cout << " \"" << decoded << "\"" << std::endl;
#endif

#if 1
            for(size_t bit=0; bit<bch_type::data_bits; bit++) {
                const auto byte_idx = bit / 8;
                const auto bit_idx = bit % 8;
                const auto lhs = decoded[byte_idx] & (1 << bit_idx);
                const auto rhs = msg_bytes[byte_idx] & (1 << bit_idx);

                if(lhs != rhs) {
                    std::cout << decoded << " != " << msg_bytes << " (FAIL)" << std::endl;
                    assert(false);
                }
            }
#endif
        }
    }

    void operator() (unsigned n_random_errors, unsigned n_times) const {
        for(size_t i=0; i<n_times; i++)
            (*this)(n_random_errors);
    }
};

template<typename>
struct test_bch_iterator;

template<unsigned m,
         unsigned t,
         unsigned...PrimitivePolynomialCoeffs>
struct test_bch_iterator<mr::bch<m, t, PrimitivePolynomialCoeffs...>> {
    using bch_type = mr::bch<m, t, PrimitivePolynomialCoeffs...>;
    using test_type = test_bch<bch_type>;

    void operator() (unsigned n_random_errors, unsigned n_times) const {
        for(size_t j=n_random_errors; j>0; j--)
            test_type {}
            (j, n_times);
    }
};

template<typename>
struct bch_test_benchmark;

template<unsigned m, unsigned t, unsigned...poly>
struct bch_test_benchmark<mr::bch<m, t, poly...>> {
    using bch_type = mr::bch<m, t, poly...>;
    using test_bch_type = test_bch<bch_type>;
    using print_tested_bch_type = print_bch<bch_type>;

    void operator() (auto repeats) const {
        using measurement_type = mr::measure_time<>;

        const auto started = measurement_type::start();

        test_bch_type{}
        (t, repeats);

        const auto elapsed = measurement_type::end(started);

        std::stringstream bch_type_ss;

        print_tested_bch_type{}
        (bch_type_ss);

        std::cout << bch_type_ss.str() << " " << repeats << " test iterations (encode -> add errors -> decode) took "
                  << measurement_type::seconds(elapsed) << " [s] (avg: "
                  << measurement_type::milliseconds(elapsed) / repeats << " [ms])" << std::endl << std::endl;
    }
};

void test_all();

template<typename>
struct measure_time_test;

template<typename Resolution, typename ClockType>
struct measure_time_test<mr::measure_time<Resolution, ClockType>> {
    using tested_type = mr::measure_time<Resolution, ClockType>;

    constexpr static auto seconds_to_wait = 1;
    constexpr static auto tolerance = 0.1;

    void operator() (unsigned repeats = 1) const {
        do {
            std::stringstream ss;
            ss << "testing time measurement with a " << seconds_to_wait << "s (" << tolerance << "s tolerance) pause...";

            const auto result = tested_type{} (ss.str(), [] {
                std::this_thread::sleep_for(std::chrono::seconds(seconds_to_wait)); // just one second
                return 0;
            });

            [[ maybe_unused ]] const auto time_diff = tested_type::seconds(result.duration) - seconds_to_wait;

            assert(time_diff > 0); // can't be faster than expected, too rigorous?
            assert(std::fabs(time_diff) < tolerance); // allow 100ms drift...
        } while(--repeats > 0);

        std::cout << std::endl;
    }
};
