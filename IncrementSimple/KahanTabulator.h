#pragma once
#ifndef KahanTabulator_H
#define KahanTabulator_H

#include <stdint.h>

namespace iki {
    /**
     * Function takes advantage of the Kahan accumulation procedure in tabulating the function with any stepping forward method (Euler, Runge-Kutta etc.)
     */
    template <typename T, typename Step, typename InputIterator>
    T kahan_tabulator(Step step, InputIterator result_begin, InputIterator result_end, T arg0, T f0, T darg, uint64_t loop_size = 1u) {
        T f = f0, arg = arg0;
        volatile T y = T(0.), t = T(0.), c = T(0.); //volatile against aggressive optimization
        uint64_t full_count = 0u;
        for (; result_begin != result_end; ++result_begin) {
            for (uint64_t loop_counter = 0; loop_counter != loop_size; ++loop_counter, ++full_count) {
                y = step(darg, arg, f) - c;
                t = f + y;
                c = (t - f) - y;
                f = t;
                
                arg = darg*full_count;
            }
            *result_begin = f;
        }
        return f;
    }
} /* iki */

#endif /* KahanTabulator_H */
