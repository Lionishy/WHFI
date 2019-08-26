#pragma once
#ifndef VDF_h
#define VDF_h

#include "PhysicalParameters.h"
#include "MathHelper.h"

#include <cmath>

namespace  iki { namespace whfi {
    template <typename T>
    struct VDF final {
    public:
        VDF(PhysicalParameters<T> params): p(params) { }
        
        T operator()(T vperp, T vparall) const {
            T coeff_c = std::exp(-pow<2>(vperp) * T(0.5)), coeff_h = std::exp(-pow<2>(vperp) * T(0.5)*p.TcTh_ratio );
            return
            p.nc * coeff_c * std::exp(-T(0.5) * pow<2>(vparall - p.bulk_to_term_c))
            + p.nh * pow<3>(std::sqrt(p.TcTh_ratio))*coeff_h *
                std::exp(-T(0.5)  * pow<2>(vparall * std::sqrt(p.TcTh_ratio) - p.bulk_to_term_h));
        }
        
    private:
        PhysicalParameters<T> p;
    };
} /* iki */ } /* whfi */

#endif /* VDF_h */
