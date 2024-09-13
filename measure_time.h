#ifndef MEASURE_TIME_H
#define MEASURE_TIME_H

#include <chrono>
#include <iostream>

namespace mr {

    using namespace std::chrono_literals;

    template<typename ResolutionType = std::chrono::nanoseconds,
             typename ClockType = std::chrono::high_resolution_clock>
    struct measure_time {
        using resolution_type = ResolutionType;
        using clock_type = ClockType;
        using result_value_type = double;

        static auto start() {
            return clock_type::now();
        }

        static auto end(const auto &start) {
            return clock_type::now() - start;
        }

        static auto seconds(const auto &elapsed) {
            return std::chrono::duration<result_value_type>(elapsed).count();
        }

        static auto milliseconds(const auto &elapsed) {
            return std::chrono::duration<result_value_type, std::milli>(elapsed).count();
        }

        static auto microseconds(const auto &elapsed) {
            return std::chrono::duration<result_value_type, std::micro>(elapsed).count();
        }

        static auto nanoseconds(const auto &elapsed) {
            return std::chrono::duration<result_value_type, std::nano>(elapsed).count();
        }

        template<typename F>
        static auto execute(F &&f) {
            const auto started = start();
            const auto result = std::forward<F>(f) ();
            const auto elapsed = end(started);

            return std::make_pair(result, elapsed);
        }

        template<typename F>
        static void print(const auto &elapsed) {
            std::cout << "execution took "
                      << measure_time::microseconds(elapsed) << " us"
                      << " ("
                      << measure_time::seconds(elapsed) << " s)"
                      << std::endl;
        }

        struct print_name {
            void operator() (const std::string &name) const { std::cout << name << " "; }
        };

        struct print_none {
            constexpr void operator() (const std::string &) const {}
        };

        template<typename F, typename Printer = print_none>
        constexpr auto operator() (F &&f, const std::string &name = {}) const {
            const auto result = execute(std::forward<F>(f));

            Printer{}
            (name);

            measure_time::print<F>(result.second);

            return result.first;
        }

        template<typename F>
        auto operator() (const std::string &name, F &&f) const {
            return measure_time::operator() <F, print_name>(std::forward<F>(f), name);
        }
    };
}

#endif // MEASURE_TIME_H
