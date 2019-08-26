#pragma once
#ifndef VDFImport_h
#define VDFImport_h

#include "Table.h"
#include "TableIO.h"

#include <iostream>
#include <fstream>
#include <string>

namespace iki { namespace whfi {
    template <typename T>
    UniformGridTable<T,2u,1u> VDFImportUniformBinary(std::string VDFImportFile) {
        std::ifstream binary_stream;
        binary_stream.exceptions(std::ios::badbit|std::ios::failbit);
        
        binary_stream.open(VDFImportFile,std::ios::binary);
        auto grid = UniformGrid<T,2u>();
        read_axis_binary(binary_stream, grid[0]);
        read_axis_binary(binary_stream, grid[1]);
        
        auto table = UniformGridTable<T,2u,1u>(grid);
        read_table_binary(binary_stream, table);
        
        return table;
    }
} /* whfi */ } /* iki */

#endif /* VDFImport_h */
