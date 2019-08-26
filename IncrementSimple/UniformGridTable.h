//
//  StepTable.h
//  AdvancedTable
//
//  Created by mposimba on 21/08/2019.
//  Copyright © 2019 mposimba. All rights reserved.
//

#pragma once
#ifndef StepTable_h
#define StepTable_h

extern "C" {
    struct fUniformGridTable {
        fUniformGrid grid;
        uint8_t value_scale;
        float *values; //.value_scale * product of grid.size[i]
    };
    
    struct dUniformGridTable {
        dUniformGrid grid;
        uint8_t value_scale;
        double *values; //.value_scale * product of grid.size[i]
    };
}

#endif /* StepTable_h */
