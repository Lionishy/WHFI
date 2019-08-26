#pragma once
#ifndef VDFTabulator_h
#define VDFTabulator_h

#include "VDF.h"
#include "Table.h"

#include <stdint.h>


namespace iki { namespace whfi {
    template <typename T>
    class VDFUniformGridTabulator final {
    public:
        VDFUniformGridTabulator(VDF<T> vdf): vdf(vdf) { }
        
        UniformGridTable<T,2u,1u> operator()(UniformGrid<T,2u> const grid) {
            auto table = UniformGridTable<T,2u,1u>(grid);
            for (uint64_t vparall_counter = 0u; vparall_counter != grid[1].size; ++vparall_counter) {
                for (uint64_t vperp_counter = 0u; vperp_counter != grid[0].size; ++vperp_counter) {
                    table[vperp_counter+vparall_counter*grid[0].size] =
                        vdf(
                            grid[0].begin + vperp_counter*grid[0].step
                            , grid[1].begin + grid[1].step*vparall_counter
                        );
                }
            }
            return table;
        }
        
    private:
        VDF<T> vdf;
    };
} /* whfi */ } /* iki */

#endif /* VDFTabulator_h */
