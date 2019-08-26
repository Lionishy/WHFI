//
//  Grid.h
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
    struct fGrid {
        uint8_t dim;
        uint64_t *size;
        float *begin; //.dim number of elements
        float *steps; //sum of size[i] number of elements
    };

    struct dGrid {
        uint8_t dim;
        uint64_t *size;
        float *begin; //.dim number of elements
        float *steps; //sum of size[i] number of elements
    };
}

#endif /* UniformGrid_h */
