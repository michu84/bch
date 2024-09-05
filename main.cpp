#include <cassert>
#include <cstring>

#include "tests.h"
#include "measure_time.h"

int main(int argc, char* argv[])
{
    const auto rand_time = time(nullptr);

#ifdef DEBUG_VERBOSE
    printf("rand_time: %li\n", rand_time);
#endif

    srand(rand_time);

    constexpr const auto m = 7, t = 13;
    const auto repeats = argc == 2 ? std::atoi(argv[1]) : 100;

    const auto started = mr::measure_time::start();

    test_m_t_coeffs<m, t, 0, 6, 7>(t, repeats);

    const auto elapsed = mr::measure_time::end(started);

    std::stringstream bch_type_ss;
    print_bch_type<m, t, 0, 6, 7>{} (bch_type_ss);

    printf("%s %d test iterations (encode -> add errors -> decode) took %f [s] (avg: %f [ms])\n",
           bch_type_ss.str().c_str(),
           repeats,
           elapsed.count() / 1000000.0,
           elapsed.count() / double(repeats) / 1000);

    return 0;
}
