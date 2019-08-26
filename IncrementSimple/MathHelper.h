#ifndef MathHelper_h
#define MathHelper_h

#include <stdint.h>
#include <optional>
#include <iostream>

namespace  iki {
    template <uint32_t P, typename T>
    T pow(T arg) {
        T res = T(1.);
        for (uint32_t counter = 0u; counter != P; ++counter)
            res *= arg;
        return res;
    }
    
    template <typename T>
    T sign(T arg) { return T(arg < T(0.) ? -1. : 1.); }
    
    template<typename T, typename FunT>
    std::optional<T> step_solve(FunT f, T begin, T end, T step) {
        uint64_t count = 0u; T curr = begin, next = begin + step, fCurr = f(curr), fNext = f(next);
        while (curr < end) {
			if (fCurr * fNext < T(0.)) {
				return T(0.5) * (curr + next);
			}
            curr = next; fCurr = fNext;
            next = begin + (++count)*step; fNext = f(next);
        }
        return std::nullopt;
    }
}

#endif /* MathHelper_h */
