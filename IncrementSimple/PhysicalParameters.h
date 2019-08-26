#pragma once
#ifndef PhysicalParameters_h
#define PhysicalParameters_h

#include <cmath>

namespace iki { namespace whfi {
    /**
     * Physical parameters to describe the velocity destribution function (VDF)
     */
    template <typename T>
    struct PhysicalParameters final {
        PhysicalParameters(T nc, T betta_c, T TcTh_ratio, T bulk_velocity_to_alfven): nc(nc), betta_c(betta_c), TcTh_ratio(TcTh_ratio), bulk_velocity_to_alfven(bulk_velocity_to_alfven)  {
            nh = T(1.) - nc;
            betta_root_c = std::sqrt(T(0.5)*betta_c);
            betta_root_h = std::sqrt(T(0.5)*betta_c/TcTh_ratio);
            bulk_to_term_c = bulk_velocity_to_alfven / betta_root_c * std::sqrt(T(1. / 1836.));
            bulk_to_term_h = -(nc/nh) * bulk_to_term_c * std::sqrt(TcTh_ratio);
        }
        
        //fundamental parameters
        T nc;                      //core particles density
        T TcTh_ratio;              //ratio of the core temperature to the halo temperature
        T betta_c;                 //ratio of the core thermal pressure to the magnetic pressure
        T bulk_velocity_to_alfven; //bulk speed in terms of alfven speed
        
        //derived parameters
        T nh;
        T betta_root_c, betta_root_h;     //square root of the betta parameters core and halo
        T bulk_to_term_c, bulk_to_term_h; //bulk velocity in terms of thermal speed
    };
} /* whfi */ } /* iki */

#endif /* PhysicalParameters_h */
