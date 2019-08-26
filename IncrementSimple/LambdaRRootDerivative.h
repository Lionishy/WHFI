#pragma once
#ifndef LambdaRRootDerivative_H
#define LambdaRRootDerivative_H

#include "PhysicalParameters.h"

namespace  iki { namespace whfi {
    template <typename T, typename FunT>
    class LambdaRRootDerivative final {
    public:
        LambdaRRootDerivative(FunT Z, PhysicalParameters<T> p): Z(Z), p(p) { }
        T operator()(T omega, T k) const {
            T arg_c = (omega - T(1.)) / (k * p.betta_root_c) - p.bulk_to_term_c;
            T arg_h = (omega - T(1.)) / (k * p.betta_root_h) - p.bulk_to_term_h;
            T Zc = Z(arg_c), Zh = Z(arg_h);
            return p.nc / (k * p.betta_root_c) * (-Zc + (omega / (k * p.betta_root_c) - p.bulk_to_term_c) * (Zc * arg_c + T(1.)))
            + p.nh / (k * p.betta_root_h) * (-Zh + (omega / (k * p.betta_root_h) - p.bulk_to_term_h) * (Zh * arg_h + T(1.)));
        }
        
    private:
        FunT Z;
        PhysicalParameters<T> p;
    };
    
    template <typename T, typename FunT>
    auto make_lambdar_rootderivative(FunT Z, PhysicalParameters<T> p) {
        return LambdaRRootDerivative(Z, p);
    }
} /* whfi */ } /* iki */


#endif /*LambdaRRootDerivative*/
