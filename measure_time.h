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
        auto operator() (const std::string &name, F &&f) const {
            const auto started = start();
            const auto result = std::forward<F>(f) ();
            const auto elapsed = end(started);

            std::cout << name
                      << " execution took "
                      << elapsed.count() << " us"
                      << " ("
                      << (elapsed.count() / 1000000.0) << " s)"
                      << std::endl;

            return result;
        }

        template<typename F>
        auto operator() (F &&f) const {
            const auto started = start();
                const auto result = std::forward<F>(f) ();
            const auto elapsed = end(started);

            std::cout << "execution took "
                      << elapsed.count() << " us"
                      << " ("
                      << (elapsed.count() / 1000000.0) << " s)"
                      << std::endl;

            return result;
        }
    };

}

#endif // MEASURE_TIME_H
