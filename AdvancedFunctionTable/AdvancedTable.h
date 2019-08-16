#pragma once
#ifndef AdvancedTable_H
#define AdvancedTable_H

#include <vector>
#include <iostream>

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

template <typename ArgT, typename ValT>
std::istream& read_table_ascii(GridMultiscalarTable<ArgT, ValT> &grid_table, std::istream &ascii_is) {
	//at first we read the grid parameters
	{
		unsigned grid_size; ascii_is >> grid_size;
		grid_table.grid.resize(grid_size);
		for (auto &g : grid_table.grid)
			ascii_is >> g.size >> g.begin >> g.step;
	}
	//next we read the values
	{
		ascii_is >> grid_table.val_size;
		//the grid has been already read, so we known the grid dimenstions
		{
			unsigned full_size = 1u;
			for (auto const &g : grid_table.grid)
				full_size *= g.size;
			grid_table.values.resize(full_size*grid_table.val_size);
		}
		for (auto &v : grid_table.values)
			ascii_is >> v;
	}
	return ascii_is;
}
	
template <typename ArgT, typename ValT>
std::ostream& write_table_binary(GridMultiscalarTable<ArgT, ValT> const &grid_table, std::ostream &binary_os) {
	using namespace std;
	//at first we write the grid parameters
	{
		unsigned int grid_size = unsigned(grid_table.grid.size());
		binary_os.write(reinterpret_cast<char const *>(&grid_size), sizeof(unsigned));
		for (auto &g : grid_table.grid) {
			binary_os.write(reinterpret_cast<char const *>(addressof(g.size)), sizeof(unsigned));
			binary_os.write(reinterpret_cast<char const *>(addressof(g.begin)), sizeof(ArgT));
			binary_os.write(reinterpret_cast<char const *>(addressof(g.step)), sizeof(ArgT));
		}
	}
	//next we write the values
	{
		binary_os.write(reinterpret_cast<char const *>(addressof(grid_table.val_size)), sizeof(unsigned));
		//values can be read as a sequential set of bytes
		binary_os.write(reinterpret_cast<char const*>(grid_table.values.data()), sizeof(ValT) * grid_table.values.size());
	}
	return binary_os;
}

template <typename ArgT, typename ValT>
std::istream& read_table_binary(GridMultiscalarTable<ArgT, ValT> &grid_table, std::istream &binary_is) {
	using namespace std;
	//at first we read the grid parameters
	{
		unsigned grid_size;
		binary_is.read(reinterpret_cast<char *>(&grid_size), sizeof(unsigned));
		grid_table.grid.resize(grid_size);
		for (auto &g : grid_table.grid) {
			binary_is.read(reinterpret_cast<char *>(addressof(g.size)), sizeof(unsigned));
			binary_is.read(reinterpret_cast<char *>(addressof(g.begin)), sizeof(ArgT));
			binary_is.read(reinterpret_cast<char *>(addressof(g.step)), sizeof(ArgT));
		}
	}
	//next we read the values
	{
		binary_is.read(reinterpret_cast<char *>(addressof(grid_table.val_size)), sizeof(unsigned));
		//the grid has been already read, so we known the grid dimenstions
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

#endif /*AdvancedTable_H*/