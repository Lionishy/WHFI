//
//  GridTable.h
//  AdvancedTable
//
//  Created by mposimba on 21/08/2019.
//  Copyright © 2019 mposimba. All rights reserved.
//
#pragma once
#ifndef GridTable_h
#define GridTable_h

#include "Grid.h"

extern "C" {
    struct fGridTable {
        fGrid grid;
        uint8_t value_scale;
        float *values; //.value_scale * product of grid.size[i]
    };

    struct dGridTable {
        dGrid grid;
        uint8_t value_scale;
        double *values; //.value_scale * product of grid.size[i]
    };
}

#endif /* GridTable_h */
