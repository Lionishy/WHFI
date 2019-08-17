#pragma once
#ifndef AdvancedTableIO_H
#define AdvancedTableIO_H

#include "AdvancedTable.h"

#include <iostream>

/**
 * Function to write a GridMultiscalarTable instance into an ascii stream
 */
template <typename ArgT, typename ValT>
std::ostream& write_table_ascii(GridMultiscalarTable<ArgT, ValT> const &grid_table, std::ostream &ascii_os) {
	//at first we need to write grid
	ascii_os << unsigned(grid_table.grid.size()) << '\n'; //one need to know how many lines to read for the grid data
	for (auto const &g : grid_table.grid)
		ascii_os << g.size << ' ' << g.begin << ' ' << g.step << '\n';

	//next we need to write data
	ascii_os << grid_table.val_size << '\n'; //one need to know how many values compose a value vector 
	for (unsigned counter = 0, end = unsigned(grid_table.values.size()); counter != end; ) {
		ascii_os << grid_table.values[counter++];
		ascii_os << (0 == counter % grid_table.val_size ? '\n' : ' ');
	}
	return ascii_os;
}

/**
 * Function to read a GridMultiscalarTable from an ascii stream
 * The stream is expected to be filled by the previous call 'write_table_ascii'
 */
template <typename ArgT, typename ValT>
std::istream& read_table_ascii(GridMultiscalarTable<ArgT, ValT> &grid_table, std::istream &ascii_is) {
	//at first we have to read the grid parameters
	{
		unsigned grid_size; ascii_is >> grid_size; //we need to know how many lines compose grid data
		grid_table.grid.resize(grid_size);
		for (auto &g : grid_table.grid)
			ascii_is >> g.size >> g.begin >> g.step;
	}
	//next we have to read the values
	{
		ascii_is >> grid_table.val_size; //we need to know how many values compose a one value vector
		//the grid has been already read, so we known the grid dimenstions and can calculate the number of values
		{
			unsigned full_size = 1u;
			for (auto const &g : grid_table.grid)
				full_size *= g.size;
			grid_table.values.resize(full_size * grid_table.val_size);
		}
		for (auto &v : grid_table.values)
			ascii_is >> v;
	}
	return ascii_is;
}

/**
 * Function to write a GridMultiscaleTable instance into a binary stream
 */
template <typename ArgT, typename ValT>
std::ostream& write_table_binary(GridMultiscalarTable<ArgT, ValT> const &grid_table, std::ostream &binary_os) {
	using namespace std;
	//at first we need to write the grid parameters
	{
		unsigned int grid_size = unsigned(grid_table.grid.size());
		binary_os.write(reinterpret_cast<char const *>(&grid_size), sizeof(unsigned)); //one need to know how many triplets to read
		for (auto &g : grid_table.grid) {
			binary_os.write(reinterpret_cast<char const *>(addressof(g.size)), sizeof(unsigned));
			binary_os.write(reinterpret_cast<char const *>(addressof(g.begin)), sizeof(ArgT));
			binary_os.write(reinterpret_cast<char const *>(addressof(g.step)), sizeof(ArgT));
		}
	}
	//next we need to write the values
	{
		//one need to know how many values compose a one value vector
		binary_os.write(reinterpret_cast<char const *>(addressof(grid_table.val_size)), sizeof(unsigned));
		//values can be written as a sequential set of bytes
		binary_os.write(reinterpret_cast<char const *>(grid_table.values.data()), sizeof(ValT) * grid_table.values.size());
	}
	return binary_os;
}

/**
 * Function to read a GridMultiscaleTable instance from a binary stream
 * The stream is expected to be filled by the previous call 'write_table_binary'
 */
template <typename ArgT, typename ValT>
std::istream& read_table_binary(GridMultiscalarTable<ArgT, ValT> &grid_table, std::istream &binary_is) {
	using namespace std;
	//at first we have to read the grid parameters
	{
		unsigned grid_size;
		binary_is.read(reinterpret_cast<char *>(&grid_size), sizeof(unsigned)); //we need to know how many triples to read
		grid_table.grid.resize(grid_size); //we need to resize vector of grid
		for (auto &g : grid_table.grid) {
			binary_is.read(reinterpret_cast<char *>(addressof(g.size)), sizeof(unsigned));
			binary_is.read(reinterpret_cast<char *>(addressof(g.begin)), sizeof(ArgT));
			binary_is.read(reinterpret_cast<char *>(addressof(g.step)), sizeof(ArgT));
		}
	}
	//next we have to read the values
	{
		//we need to know how many values combine a one value vector
		binary_is.read(reinterpret_cast<char *>(addressof(grid_table.val_size)), sizeof(unsigned));
		//the grid has been already read, so we known the grid dimenstions and can calculate the number of values
		{
			unsigned full_size = 1u;
			for (auto const &g : grid_table.grid)
				full_size *= g.size;
			grid_table.values.resize(full_size * grid_table.val_size);
		}
		//values can be read as a sequential set of bytes
		binary_is.read(reinterpret_cast<char *>(grid_table.values.data()), sizeof(ValT) * grid_table.values.size());
	}
	return binary_is;
}

/**
 * Function to write a MultiscalarTalbe instance into an ascii stream
 */
template <typename ArgT, typename ValT>
std::ostream &write_table_ascii(MultiscalarTable<ArgT, ValT> const &grid_table, std::ostream &ascii_os) {
	using namespace std;
	//at first we need to write the total number of values and arguments
	//and the number of arguments and values in a one (arg,val) pair
	ascii_os << unsigned(grid_table.arguments.size()) << ' ' << grid_table.arg_size << ' ' << unsigned(grid_table.values.size()) << ' ' << grid_table.val_size << '\n';
	//next we need to write all arguments
	for (auto const &arg : grid_table.arguments)
		ascii_os << arg << '\n';
	ascii_os << flush;
	//next we need to write all values
	for (auto const &val : grid_table.values)
		ascii_os << val << '\n';
	ascii_os << flush;
	return ascii_os;
}

/**
 * Function to read a MultiscalarTable instance from an ascii stream
 * The stream is expected to be filled by the previous call 'write_table_ascii'
 */
template <typename ArgT, typename ValT>
std::istream& read_table_ascii(MultiscalarTable<ArgT, ValT> &grid_table, std::istream &ascii_is) {
	using namespace std;
	//at first we have to read the total number of arguments and values
	//and the number of arguments and values in a one (arg,val) pair
	{
		unsigned argument_size, value_size;
		ascii_is >> argument_size >> grid_table.arg_size >> values_size >> grid_table.val_size;
		grid_table.arguments.resize(argument_size); //we have to resize table vectors
		grid_table.values.resize(value_size);
	}

	//next we need to read all arguments and values
	for (auto &arg : grid_table.arguments)
		ascii_is >> arg;
	for (auto &val : grid_table.values)
		ascii_is >> val;
	return ascii_is;
}

/**
 * Function to write a MultiscaleTable instance into a binary stream
 */
template <typename ArgT, typename ValT>
std::ostream& write_talbe_binary(MultiscalarTable<ArgT, ValT> const &grid_table, std::ostream &binary_os) {
	//at first we need to write the total number of arguments and values
	//and the number of arguments and values in a one (arg,val) pair
	{
		unsigned arguments_size = unsigned(grid_table.arguments.size()), values_size = unsigned(grid_table.values.size());
		binary_os.write(reinterpret_cast<char const *>(&arguments_size), sizeof(unsigned));
		binary_os.write(reinterpret_cast<char const *>(addressof(grid_table.arg_size)), sizeof(unsigned));
		binary_os.write(reinterpret_cast<char const *>(&values_size), sizeof(unsigned));
		binary_os.write(reinterpret_cast<char const *>(addressof(grid_table.val_size)), sizeof(unsigned));
	}

	//next we need to write all the arguments and the values
	binary_os.write(reinterpret_cast<char const *>(grid_table.arguments.data()), sizeof(ArgT) * grid_table.arguments.size());
	binary_os.write(reinterpret_cast<char const *>(grid_table.values.data()), sizeof(ValT) * grid_table.values.size());
	return binary_os;
}

/**
 * Function to read a MultiscaleTable instance from a binary stream
 * The stream is expected to be filled by the previous call 'write_table_binary'
 */
template <typename ArgT, typename ValT>
std::istream& write_talbe_binary(MultiscalarTable<ArgT, ValT> &grid_table, std::istream &binary_is) {
	//at first we have to read the total number of arguments and values
	//and the number of arguments and values in a one (arg,val) pair
	{
		unsigned arguments_size, values_size;
		binary_is.read(reinterpret_cast<char *>(&arguments_size), sizeof(unsigned));
		binary_is.read(reinterpret_cast<char *>(addressof(grid_table.arg_size)), sizeof(unsigned));
		binary_is.read(reinterpret_cast<char *>(&values_size), sizeof(unsigned));
		binary_is.read(reinterpret_cast<char *>(addressof(grid_table.val_size)), sizeof(unsigned));
		grid_table.arguments.resize(arguments_size);
		grid_table.values.resize(values_size);
	}

	//next we have to read all the arguments and the values
	//we can read all the values as a sequence
	binary_is.read(reinterpret_cast<char *>(grid_table.arguments.data()), sizeof(ArgT) * grid_table.arguments.size());
	binary_is.read(reinterpret_cast<char *>(grid_table.values.data()), sizeof(ValT) * grid_table.values.size());
	return binary_is;
}

#endif /*AdvancedTableIO_H*/