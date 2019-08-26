#pragma once
#ifndef WDRImport_h
#define WDRImport_h

#include "Table.h"

#include <string>
#include <iostream>
#include <fstream>

namespace  iki { namespace whfi {
    template <typename T>
    ArgumentGridTable<T,1u,3u> WDRImportTableBinary(std::string WDRTableFilename) {
        ArgumentGrid<T,1u> grid;
        std::ifstream binary_stream;
        binary_stream.exceptions(std::ios::badbit|std::ios::failbit);
        
        binary_stream.open(WDRTableFilename,std::ios::binary);
        read_axis_binary(binary_stream, grid[0]);
        auto table = ArgumentGridTable<T, 1u, 3u>(grid);
        read_table_binary(binary_stream, table);
        return table;
    }
} /* whfi */ } /* iki */

#endif /* WDRImport_h */
