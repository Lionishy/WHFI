#pragma once
#ifndef ZFunc_H
#define ZFunc_H

#include "KahanTabulator.h"
#include "Runge.h"
#include "Table.h"

#include <vector>
#include <stdint.h>

namespace iki {
    /**
     * Function creates a vector of values of the plasma dispersion function dZ/dx = -x*Z - 1 for the arguments from 0 to the max_argument with the darg step
     * darg and max_argument are expected to be positive
     */
    template <typename T>
    UniformGridTable<T> ZFunc_tabulator(T darg, T max_argument, uint64_t loop_count = 1u) {
        using namespace std;
        uint64_t step_count = uint64_t(max_argument / (darg*loop_count)) + 1u;
        UniformGridTable<T> table({{step_count,T(0.),darg*loop_count}});
        *table.begin() = T(0);
        kahan_tabulator<T>(
           make_runge4thorder<T>([](T arg, T val)->T { return -arg * val - T(1.); })
           , ++table.begin(), table.end()
           , T(0), T(0), darg
           , loop_count
        );
        return table;
    }
}

#endif
