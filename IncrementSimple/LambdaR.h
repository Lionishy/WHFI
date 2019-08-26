#pragma once
#ifndef LambdaR_H
#define LambdaR_H

#include "PhysicalParameters.h"

#include <iostream>

namespace iki { namespace whfi {
    template <typename T, typename FunT>
    class LambdaR final {
    public:
        LambdaR(FunT Z, PhysicalParameters<T> p): Z(Z), p(p) { }
        
        T operator()(T omega, T k) const { // 1./(omega^2) is omitted deliberately
            auto res =  T(1. / 1836.) + k * k
            - p.nc * (omega / (k * p.betta_root_c) - p.bulk_to_term_c) * Z( (omega - T(1.)) / (k * p.betta_root_c) - p.bulk_to_term_c)
            - p.nh * (omega / (k * p.betta_root_h) - p.bulk_to_term_h) * Z( (omega - T(1.)) / (k * p.betta_root_h) - p.bulk_to_term_h);
			return res;
        }
        
    private:
        FunT Z;
        PhysicalParameters<T> p;
    };
    
    /**
     * Function to construct LambdaR with the FunT parameter deduced automaticly
     */
    template <typename T, typename FunT>
    auto make_lambdar(FunT Z, PhysicalParameters<T> p) {
        return LambdaR<T, FunT>(Z, p);
    }
} /* whfi */ } /* iki */

#endif /*LambdaR_H*/
