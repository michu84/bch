#include <cassert>
#include <cstdlib>
#include <cstring>

#include "tests.h"

// ATSC A/336 audio watermarking BCH(127,50,13) specification: m=7, n=127, k=50, t=13, primitive = 1 + x^6 + x^7, generator is a very long 77th order polynomial:
// 1 + x^2 + x^5 + x^9 + x^12 + x^13 +x^17 + x^18 + x^19 + x^20 + x^21 + x^26 + x^29 + x^30 + x^32
//   + x^34 + x^35 + x^39 + x^40 + x^41 + x^42 + x^44 + x^49 + x^50 + x^51 + x^59 + x^60 + x^62
//   + x^63 + x^64 + x^66 + x^67 + x^68 + x^71 + x^72 + x^74 + x^75 + x^76 + x^77

using ATSC_A336_codec_type = mr::bch<7, 13, 0, 6, 7>;

// Phobos Lander telemetry BCH(127,113,2) codes: m=7, n=127, k=113, t=2, primitive = 1 + x^3 + x^7, generator = 1 + x^1 + x^2 + x^4 + x^5 + x^6 + x^8 + x^9 + x^14

using phobos_lander_codec_type = mr::bch<7, 2, 0, 3, 7>;

int main(int argc, char* argv[])
{
    const auto rand_time = time(nullptr);

#ifdef DEBUG_VERBOSE
    std::cout << "rand_time: " << rand_time << std::endl;
#endif

    srand(rand_time);

    const auto repeats = argc == 2 ? std::atoi(argv[1]) : 100;

    // measure_time_test<mr::measure_time<>>{}
    // (3);

    bch_tester<ATSC_A336_codec_type>{}
    (repeats);

    // bch_tester<phobos_lander_codec_type>{}
    // (repeats);

    return 0;
}
