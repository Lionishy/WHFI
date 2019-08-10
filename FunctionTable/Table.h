#pragma once
#ifndef Table_H
#define Table_H

#include <vector>

/**
 * The structure to represent a table function as a number of pairs (argument,value)
 */
template<typename ArgT, typename ValT>
struct ArgValueTable final {
	std::vector<std::pair<ArgT, ValT>> table;
};

/**
 * The structure to represent a table with a fixed argument step and a set of values 
 */
template <typename ArgT, typename ValT>
struct StepArgumentTable final {
	ArgT arg0, darg;
	std::vector<ValT> table;
};

#endif /*Table_H*/