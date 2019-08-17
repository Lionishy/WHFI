#pragma once
#ifndef AdvancedTable_H
#define AdvancedTable_H

#include <vector>

template <typename ArgT, typename ValT> 
struct MultiscalarTable final {
	unsigned arg_size, val_size;
	std::vector<ArgT> arguments;
	std::vector<ValT> values;
};

template <typename ArgT>
struct Grid final {
	unsigned size;
	ArgT begin, step;
};
	
template <typename ArgT, typename ValT>
struct GridMultiscalarTable final {
	std::vector<Grid<ArgT>> grid;
	unsigned val_size;
	std::vector<ValT> values;
};

#endif /*AdvancedTable_H*/