#ifndef MEASURE_TIME_H
#define MEASURE_TIME_H

#include <chrono>
#include <iostream>

namespace mr {

    using namespace std::chrono_literals;

    struct measure_time {
        static auto start() {
            return std::chrono::high_resolution_clock::now();
        }

        template<class result_t = std::chrono::milliseconds,
                 class clock_t = std::chrono::steady_clock,
                 class duration_t = std::chrono::milliseconds>
        static auto since(std::chrono::time_point<clock_t, duration_t> const &start) {
            return std::chrono::duration_cast<result_t>(clock_t::now() - start);
        }

        static auto end(const auto &start) {
            return since<std::chrono::microseconds,
                         std::chrono::high_resolution_clock,
                         std::chrono::duration<double>>(start);
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
                      << elapsed.count() << " us"
                      << " ("
                      << (elapsed.count() / 1000000.0) << " s)"
                      << std::endl;
        }

        struct print_name {
            void operator() (const std::string &name) const { std::cout << name << " "; };
        };

        struct print_none {
            void operator() (const std::string &) const {};
        };

        template<typename F, typename Printer = print_none>
        auto operator() (F &&f, const std::string &name = {}) const {
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
