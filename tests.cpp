#include "tests.h"
#include <assert.h>
#include "polynomial.h"

#define __use_constexpr_poly_division

#ifdef __use_constexpr_poly_division
#define __assert_type static_assert
#define __call_ctx constexpr
#else
#define __assert_type assert
#define __call_ctx
#endif

void test_m5()
{
    constexpr auto m = 5;
    using poly_type = mr::polynomial<bit_t, m>;
    constexpr poly_type _x5_x2_1(
                            1, 0,
                            1, 2,
                            1, 5
                        ),
                        _0,
                        _1(1, 0),
                        _x(1, 1),
                        _x5(1, 5);

    static_assert(_x5_x2_1.degree() == m);
    static_assert(_x.degree() == 1);
    static_assert(_1.degree() == 0);
    static_assert(_x5.degree() == 5);

    __call_ctx auto div_result3 = _x / _x5_x2_1;

    __assert_type(div_result3.q == _0);
    __assert_type(div_result3.r == _x);

    __call_ctx auto div_result2 = _x5_x2_1 / _1;

    __assert_type(div_result2.q == _x5_x2_1);
    __assert_type(div_result2.r == _0);

    __call_ctx auto div_result_ = _x5 / _x5_x2_1;
    constexpr auto expected_div_remainder = poly_type(1, 0, 1, 2);

    __assert_type(div_result_.q == _1);
    __assert_type(div_result_.r == expected_div_remainder);

    std::cout << "hello" << std::endl;
    std::cout << _x5_x2_1 << std::endl;
    std::cout << std::endl;

    using gf_poly_type = mr::polynomial<bit_t, 31>;
    constexpr gf_poly_type gf_x5_x2_1(1,0, 1,2, 1,5);
    for(unsigned i=0; i<(1 << m); i++) {
        const auto element_poly = gf_poly_type(1, i);
        const auto test = element_poly / gf_x5_x2_1;
#if 0
        std::cout << std::format("[{}]: ", i) << element_poly << " / " << gf_x5_x2_1 << ": q=" << test.q << ", r=" << test.r << std::endl;
#else
        mr::polynomial<bit_t,m-1> final_r;
        polynomial_copy_unsafe(test.r, final_r);

        printf("[a^%2d]: %s / %s: q=%s, r=%s\n",
            i,
            element_poly.to_string().c_str(),
            gf_x5_x2_1.to_string().c_str(),
            test.q.to_string().c_str(),
            final_r.to_string().c_str()
        );
#endif
    }
    std::cout << std::endl;
}

void test_m7()
{
    constexpr auto m = 7;
    using gf_poly_type = mr::polynomial<bit_t, 127>;
    constexpr gf_poly_type gf_x7_x6_1(1,0, 1,6, 1,7);
    for(unsigned i=0; i<(1 << m); i++) {
        const auto element_poly = gf_poly_type(1, i);
        const auto test = element_poly / gf_x7_x6_1;
#if 0
            std::cout << std::format("[{}]: ", i) << element_poly << " / " << gf_x7_x6_1 << ": q=" << test.q << ", r=" << test.r << std::endl;
#else
        mr::polynomial<bit_t,m-1> final_r;
        polynomial_copy_unsafe(test.r, final_r);

        printf("[a^%3d]: %s / %s: q=%s, r=%s\n",
            i,
            element_poly.to_string().c_str(),
            gf_x7_x6_1.to_string().c_str(),
            test.q.to_string().c_str(),
            final_r.to_string().c_str()
        );
#endif
    }
    std::cout << std::endl;
}

template<unsigned m, unsigned t>
void test_all() {
    constexpr mr::gen_seq<m> m_seq;
    constexpr mr::gen_seq<t> t_seq;

    // TODO: need an automatic iteration over all m and t with finding primitive poly for gf in each separately
    assert(false);
}
