//
//  Table.h
//  AdvancedTable
//
//  Created by mposimba on 21/08/2019.
//  Copyright © 2019 mposimba. All rights reserved.
//

#ifndef Table_h
#define Table_h


#include "UniformGrid.h"
#include "UniformGridTable.h"

#include <iterator>
#include <vector>
#include <array>
#include <initializer_list>
#include <algorithm>
#include <utility>
#include <stdint.h>

namespace iki {
    template <typename T>
    struct Axis final {
        uint64_t size = 0u; T begin = T(0.), step = T(1.);
    };
    
    template <typename T, uint8_t Dim = 1u>
    class UniformGrid final {
    public:
        UniformGrid() { }
        template <typename Iterator>
        UniformGrid(Iterator begin, Iterator end) { std::copy(begin,end,std::begin(axes)); }
        UniformGrid(std::initializer_list<Axis<T>> al): UniformGrid(std::begin(al),std::end(al)) { }
        
        Axis<T>& operator[](uint8_t dim) { return axes[dim]; }
        Axis<T> operator[](uint8_t dim) const { return axes[dim]; }
        
        auto begin() { return axes.begin(); }
        auto begin() const { return axes.cbegin(); }
        auto end() { return axes.end(); }
        auto end() const { return axes.cend(); }
        
        uint64_t values_size() const {
            uint64_t full_size = 1;
            for (auto const &axis : axes)
                full_size *= axis.size;
            return full_size;
        }
    private:
        std::array<Axis<T>,Dim> axes;
    };
    
    template <typename T, uint8_t Dim = 1u, uint8_t Scale = 1u>
    class UniformGridTable final {
    public:
        UniformGridTable(UniformGrid<T,Dim> grid): grid(grid), values(Scale*grid.values_size()) { }
        
        auto begin() { return values.begin(); }
        auto begin() const { return values.cbegin(); }
        auto end() { return values.end(); }
        auto end() const { return values.cend(); }
        
        uint64_t values_count() const { return values.size(); }
        T operator[](uint64_t idx) const { return values[idx]; }
        T& operator[](uint64_t idx) { return values[idx]; }
        T* values_data() { return values.data(); }
        T const* values_data() const { return values.data(); }
        
        std::array<T,Dim> argument_of(uint64_t idx) const {
            idx /= Scale;
            std::array<T,Dim> argument;
            std::transform(grid.begin(),grid.end(),argument.begin(),[&idx](auto const &axis){
                auto arg = axis.begin + idx%axis.size*axis.step;
                idx /= axis.size;
                return arg;
            });
            return argument;
        }
        
        UniformGrid<T,Dim> getGrid() const { return grid; }
        UniformGrid<T,Dim>& getGrid() { return grid; }
        UniformGridTable<T,Dim,Scale>& setGrid(UniformGrid<T,Dim> const &new_grid) {
            grid = new_grid;
            values.resize(Scale*grid.values_size());
            return *this;
        }
        
    private:
        UniformGrid<T,Dim> grid;
        std::vector<float> values;
    };
    
    /* AgumentTable */
    template <typename T>
    class ArgumentAxis final {
    public:
        ArgumentAxis() { }
        ArgumentAxis(uint64_t size): arguments(size) { }
        template <typename Iterator>
        ArgumentAxis(Iterator begin, Iterator end) {
            std::copy(begin,end,back_inserter(arguments));
        }
        
        ArgumentAxis& resize(uint64_t new_size) { arguments.resize(new_size); return *this; }
        
        auto begin() { arguments.begin(); }
        auto end() { arguments.end(); }
        auto begin() const { arguments.cbegin(); }
        auto end() const { arguments.cend(); }
        
        uint64_t size() const { return arguments.size(); }
        T operator[](uint64_t idx) const { return arguments[idx]; }
        T& operator[](uint64_t idx) { return arguments[idx]; }
        T* arguments_data() { return arguments.data(); }
        T const* arguments_data() const { return arguments.data(); }
        
    private:
        std::vector<T> arguments;
    };
    
    template <typename T, uint8_t Dim>
    class ArgumentGrid final {
    public:
        ArgumentGrid() { }
        template <typename Iterator>
        ArgumentGrid(Iterator begin, Iterator end) { std::copy(begin,end,std::begin(axes)); }
        ArgumentGrid(std::initializer_list<ArgumentAxis<T>> al): ArgumentGrid(std::begin(al),std::end(al)) { }
        
        ArgumentAxis<T> operator[](uint8_t idx) const { return axes[idx]; }
        ArgumentAxis<T>& operator[](uint8_t idx) { return axes[idx]; }
        
        auto begin() { return axes.begin(); }
        auto begin() const { return axes.cbegin(); }
        auto end() { return axes.end(); }
        auto end() const { return axes.cend(); }
        
        uint64_t values_size() {
            uint64_t full_size = 1u;
            for (auto const &axis : axes)
                full_size *= axis.size();
            return full_size;
        }
        
    private:
        std::array<ArgumentAxis<T>,Dim> axes;
    };
    
    template <typename T, uint8_t Dim, uint8_t Scale>
    class ArgumentGridTable final {
    public:
        ArgumentGridTable(ArgumentGrid<T,Dim> grid): grid(grid), values(Scale*grid.values_size()) { }
        
        auto begin() { return values.begin(); }
        auto begin() const { return values.cbegin(); }
        auto end() { return values.end(); }
        auto end() const { return values.cend(); }
        
        uint64_t values_count() const { return values.size(); }
        T operator[](uint64_t idx) const { return values[idx]; }
        T& operator[](uint64_t idx) { return values[idx]; }
        T* values_data() { return values.data(); }
        T const* values_data() const { return values.data(); }
        
        std::array<T,Dim> argument_of(uint64_t idx) const {
            idx /= Scale;
            std::array<T,Dim> argument;
            std::transform(grid.begin(),grid.end(),argument.begin(),[&idx](auto const &axis){
                auto arg = axis[idx%axis.size()];
                idx /= axis.size();
                return arg;
            });
            return argument;
        }
        
        ArgumentGrid<T,Dim> getGrid() const { return grid; }
        ArgumentGrid<T,Dim>& getGrid() { return grid; }
        ArgumentGridTable<T,Dim,Scale>& setGrid(ArgumentGrid<T,Dim> const &new_grid) {
            grid = new_grid;
            values.resize(Scale*grid.values_size());
            return *this;
        }
    private:
        ArgumentGrid<T, Dim> grid;
        std::vector<T> values;
    };
}

#endif /* Table_h */
