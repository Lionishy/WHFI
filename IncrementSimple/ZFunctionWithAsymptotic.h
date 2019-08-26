#pragma once
#ifndef ZFunctionWithAsymptotic_H
#define ZFunctionWithAsymptotic_H

#include "Table.h"
#include "MathHelper.h"

#include <vector>
#include <memory>
#include <cmath>

namespace  iki {
    template <typename T>
    class ZFuncWithAsymptotic final {
    public:
        ZFuncWithAsymptotic(std::shared_ptr<UniformGridTable<T,1u,1u> const> known_values): known_values(known_values), step(known_values->getGrid()[0].step) { }
        
        T operator()(T arg) const {
            T farg = std::abs(arg);
            auto idx = uint64_t(farg / step);
            if ((idx + 1u) < known_values->values_count()) {
                return iki::sign(arg) * ( ((*known_values)[idx + 1u] - (*known_values)[idx]) / step * (farg - step * idx) + (*known_values)[idx] );
            }
            else {
                T over = T(1.) / arg, square = over * over;
                return -over * (T(1.) + square + T(3.) * square * square);
            }
        }
        
    private:
        T const step;
        std::shared_ptr<UniformGridTable<T,1u,1u> const> known_values;
    };
} /* iki */

#endif /*ZFunctionWithAsymptotic_H*/
