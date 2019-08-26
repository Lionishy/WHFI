#pragma once
#ifndef KahanVDFMoments_h
#define KahanVDFMoments_h

#include "Table.h"
#include "MathHelper.h"

#include <exception>
#include <stdexcept>
#include <iterator>

namespace iki { namespace whfi {
    template <typename T, typename Iterator>
    T gUniformGrid(Iterator begin, Iterator end, Axis<T> const vperp_axis) {
        if (begin == end || std::next(begin) == end) throw std::invalid_argument("an empty value range to integrate over");
        T g = T(0.);
        volatile T y = T(0.), t = T(0.), c = T(0.);
        uint64_t count = 0u;
        for (auto curr = begin, next = std::next(begin); next != end; ++curr, ++next, ++count) {
            T arg = vperp_axis.begin + vperp_axis.step*count;
            T integral = T(0.5)*(*curr+*next)*arg*vperp_axis.step + T(1./3.)*(*next)*pow<2>(vperp_axis.step) + T(1./6.)*(*curr)*pow<2>(vperp_axis.step);
            y = integral - c;
            t = g + y;
            c = (t - g) - y;
            g = t;
        }
        return g;
    }
    
    template <typename T, typename Iterator>
    T GUniformGrid(Iterator begin, Iterator end, Axis<T> const vperp_axis) {
        if (begin == end || std::next(begin) == end) throw std::invalid_argument("an empty value range to integrate over");
        T G = T(0.);
        volatile T y = T(0.), t = T(0.), c = T(0.);
        uint64_t count = 0u;
        for (auto curr = begin, next = std::next(begin); next != end; ++curr, ++next, ++count) {
            T arg = vperp_axis.begin + vperp_axis.step*count;
            
            T integral = T(1./40.)*vperp_axis.step*(
                T(10.)*pow<3>(arg)*(*curr+*next) + T(10.)*pow<2>(arg)*vperp_axis.step*(*curr + *next*T(2.)) + T(5.)*arg*pow<2>(vperp_axis.step)*(*curr + *next*T(3.)) + pow<3>(vperp_axis.step)*(*curr + *next*T(4.))
                                                    
            );
            
            //T integral = T(0.5)*(*curr*T(0.5)*pow<3>(arg) + *next*T(0.5)*pow<3>(vperp_axis.begin + vperp_axis.step*(count+1u)))*vperp_axis.step;
            y = integral - c;
            t = G + y;
            c = (t - G) - y;
            G = t;
        }
        return G;
    }
} /* whfi */ } /* iki */

#endif /* KahanVDFMoments_h */
