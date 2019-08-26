#pragma once
#ifndef VDFExport_h
#define VDFExport_h

#include "Table.h"
#include "TableIO.h"
#include "PhysicalParameters.h"
#include "VDF.h"
#include "VDFTabulator.h"

#include <iostream>
#include <fstream>
#include <iomanip>

namespace iki { namespace whfi {
    template <typename T>
    UniformGridTable<T,2u,1u> VDFExportUniformTable(PhysicalParameters<T> params, UniformGrid<T,2u> grid) {
        return VDFUniformGridTabulator<T>(VDF<T>(params))(grid);
    }
    
    template <typename T>
    void VDFExportUnifromAscii(std::string VDFTableFile, PhysicalParameters<T> params, UniformGrid<T,2u> grid) {
        std::ofstream ascii_stream;
        ascii_stream.exceptions(std::ios::badbit|std::ios::failbit);
        
        ascii_stream.open(VDFTableFile);
        ascii_stream << std::setprecision(7) << std::fixed << VDFExportUniformTable(params, grid);
    }
    
    template <typename T>
    void VDFExportUniformBinary(std::string VDFTableFile, PhysicalParameters<T> params, UniformGrid<T,2u> grid) {
        std::ofstream binary_stream;
        binary_stream.exceptions(std::ios::badbit|std::ios::failbit);
        
        binary_stream.open(VDFTableFile,std::ios::binary);
        write_table_binary(binary_stream, VDFExportUniformTable(params, grid));
    }
} /* whfi */ } /* iki */

#endif /* VDFExport_h */
