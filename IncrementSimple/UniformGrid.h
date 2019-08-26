//
//  UniformGrid.hpp
//  AdvancedTable
//
//  Created by mposimba on 21/08/2019.
//  Copyright © 2019 mposimba. All rights reserved.
//

#pragma once
#ifndef UniformGrid_h
#define UniformGrid_h

#include <stddef.h>
#include <stdint.h>

extern "C" {
    struct fAxis {
        uint64_t size; float begin, step;
    };
    
    struct fUniformGrid {
        uint8_t dim;
        fAxis *axes; //.dim number of elements
    };
    
    struct dAxis {
        uint64_t size; float begin, step;
    };
    
    struct dUniformGrid {
        uint8_t dim;
        dAxis *axes; //.dim number of elements
    };
}

#endif /* UniformGrid_h */
