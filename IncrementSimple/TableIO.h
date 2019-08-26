//
//  TableIO.h
//  AdvancedTable
//
//  Created by mposimba on 21/08/2019.
//  Copyright © 2019 mposimba. All rights reserved.
//

#pragma once
#ifndef TableIO_h
#define TableIO_h

#include "Table.h"

#include <iostream>
#include <stdint.h>

template <typename T, uint8_t Dim, uint8_t Scale>
std::ostream& operator<<(std::ostream &ascii_os, iki::UniformGridTable<T,Dim, Scale> const &table) {
    uint64_t idx = 0u;
    for (auto it = table.begin(), end = table.end(); it != end;) {
        auto argument = table.argument_of(idx);
        for (auto const &arg : argument)
            ascii_os << arg << ' ';
        for (uint8_t val_idx = 0u; val_idx != Scale; ++val_idx, ++idx, ++it)
            ascii_os << *it << ' ';
        ascii_os << '\n';
    }
    return ascii_os << std::flush;
}

template <typename T, uint8_t Dim, uint8_t Scale>
std::istream& operator>>(std::istream &ascii_is, iki::UniformGridTable<T,Dim, Scale> &table) {
    for (auto it = table.begin(), end = table.end(); it != end;) {
        T arg;
        for (uint8_t arg_count = 0u; arg_count != Dim; ++arg_count)
            ascii_is >> arg;
        for (uint8_t val_count = 0u; val_count != Scale; ++val_count, ++it)
            ascii_is >> *it;
    }
    return ascii_is;
}

template <typename T, uint8_t Dim, uint8_t Scale>
std::ostream& operator<<(std::ostream &ascii_os, iki::ArgumentGridTable<T,Dim,Scale> const &table) {
    for (uint64_t idx = 0u; idx != table.values_count(); ) {
        for (auto arg : table.argument_of(idx))
            ascii_os << arg << ' ';
        for (uint8_t scale_count = 0u; scale_count != Scale; ++scale_count, ++idx)
            ascii_os << table[idx] << ' ';
        ascii_os << '\n';
    }
    return ascii_os << std::flush;
}

template <typename T, uint8_t Dim, uint8_t Scale>
std::istream& operator>>(std::istream &ascii_in, iki::ArgumentGridTable<T,Dim,Scale> &table) {
    
    return ascii_in;
}

namespace iki {
    template <typename T>
    std::ostream& write_axis_binary(std::ostream &binary_os, Axis<T> const &axis) {
        binary_os.write(reinterpret_cast<char const *>(std::addressof(axis.size)),sizeof(uint64_t));
        binary_os.write(reinterpret_cast<char const *>(std::addressof(axis.begin)),sizeof(T));
        binary_os.write(reinterpret_cast<char const *>(std::addressof(axis.step)),sizeof(T));
        return binary_os << std::flush;
    }
    
    template <typename T, uint8_t Dim>
    std::ostream& write_grid_binary(std::ostream &binary_os, UniformGrid<T,Dim> const &grid) {
        for (auto const &axis : grid)
            write_axis_binary(binary_os, axis);
        return binary_os << std::flush;
            
    }
    
    template <typename T, uint8_t Dim, uint8_t Scale>
    std::ostream& write_table_binary(std::ostream &binary_os, UniformGridTable<T,Dim, Scale> const &table) {
        write_grid_binary(binary_os,table.getGrid());
        binary_os.write(reinterpret_cast<char const *>(table.values_data()), sizeof(T)*table.values_count());
        return binary_os << std::flush;
    }
    
    template <typename T>
    std::istream& read_axis_binary(std::istream &binary_is, Axis<T> &axis) {
        binary_is.read(reinterpret_cast<char*>(std::addressof(axis.size)), sizeof(uint64_t));
        binary_is.read(reinterpret_cast<char*>(std::addressof(axis.begin)),sizeof(T));
        binary_is.read(reinterpret_cast<char*>(std::addressof(axis.step)),sizeof(T));
        return binary_is;
    }
    
    template <typename T, uint8_t Dim, uint8_t Scale>
    std::istream& read_table_binary(std::istream &binary_is, UniformGridTable<T,Dim,Scale> &table) {
        binary_is.read(reinterpret_cast<char*>(table.values_data()),sizeof(T)*table.values_count());
        return binary_is;
    }
    
    template <typename T>
    std::ostream& write_axis_binary(std::ostream &binary_os, ArgumentAxis<T> const &axis) {
        uint64_t axis_size = axis.size();
        binary_os.write(reinterpret_cast<char const*>(&axis_size),sizeof(uint64_t));
        binary_os.write(reinterpret_cast<char const*>(axis.arguments_data()),sizeof(T)*axis_size);
        return binary_os;
    }
    
    template <typename T, uint8_t Dim>
    std::ostream& write_grid_binary(std::ostream &binary_os, ArgumentGrid<T,Dim> const &grid) {
        for (auto const &axis : grid)
            write_axis_binary(binary_os,axis);
        return binary_os;
    }
    
    template <typename T, uint8_t Dim, uint8_t Scale>
    std::ostream& write_table_binary(std::ostream &binary_os, ArgumentGridTable<T,Dim,Scale> const &table) {
        write_grid_binary(binary_os, table.getGrid());
        binary_os.write(reinterpret_cast<char const*>(table.values_data()),sizeof(T)*table.values_count());
        return binary_os;
    }
    
    template <typename T>
    std::istream& read_axis_binary(std::istream &binary_is, ArgumentAxis<T> &axis) {
        {
            uint64_t axis_size;
            binary_is.read(reinterpret_cast<char*>(&axis_size),sizeof(uint64_t));
            axis.resize(axis_size);
        }
        
        binary_is.read(reinterpret_cast<char*>(axis.arguments_data()),sizeof(T)*axis.size());
        return binary_is;
    }
    
    template <typename T, uint8_t Dim, uint8_t Scale>
    std::istream& read_table_binary(std::istream &binary_is, ArgumentGridTable<T,Dim,Scale> &table) {
        binary_is.read(reinterpret_cast<char*>(table.values_data()),sizeof(T)*table.values_count());
        return binary_is;
    }
}

#endif /* TableIO_h */
