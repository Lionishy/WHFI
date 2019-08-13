#pragma once
#ifndef PhysicalParameters_h
#define PhysicalParameters_h

/**
 * Physical parameters to describe the velocity destribution function (VDF)
 */
template <typename T>
struct PhysicalParameters final {
    T nc, nh;                         //density of the core and the halo
    T betta_root_c, betta_root_h;     //square root of the betta parameter (magnetic pressure to thermal pressure ratio) over 2
    T bulk_to_term_c, bulk_to_term_h; //ratio of the bulk velocity to termal velocity
};

#endif /* PhysicalParameters_h */
