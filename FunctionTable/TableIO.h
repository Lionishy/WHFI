#pragma once
#ifndef TableIO_H
#define TableIO_H

#include "Table.h"

#include <iostream>

/**
 * Ouput functions
 */
/**
 * Function writes an instance of the ArgValueTalbe struct to an ascii stream
 * The first line is the number of (arg,val) pairs
 */
template <typename ArgT, typename ValT>
std::ostream& write_table_ascii(ArgValueTable<ArgT, ValT> const &table, std::ostream &ascii_out) {
	ascii_out << static_cast<unsigned int>(table.table.size()) << '\n';
	for (auto const &arg_value_pair : table.table)
		ascii_out << arg_value_pair.first << ' ' << arg_value_pair.second << '\n';
	return ascii_out << std::flush;
}

/**
 * Function writes an instance of the StepArgumentTable struct to an ascii stream 
 * The first line is the number of values
 * The second and the third lines are the initial argument and the argument step, correspondingly
 */
template <typename ArgT, typename ValT>
std::ostream& write_table_ascii(StepArgumentTable<ArgT, ValT> const &table, std::ostream &ascii_out) {
	ascii_out << static_cast<unsigned int>(table.table.size()) << '\n';
	ascii_out << table.arg0 << '\n' << table.darg << '\n';
	for (auto const &value : table.table)
		ascii_out << value << '\n';
	return ascii_out << std::flush;
}

/**
 * Function writes an instance of the ArgValueTable struct to a binary stream
 * The first several bytes encode the number of (arg,val) pairs
 * The data are written as represented in memory
 */
template <typename ArgT, typename ValT>
std::ostream& write_table_binary(ArgValueTable<ArgT, ValT> const &table, std::ostream &binary_out) {
	unsigned int size = static_cast<unsigned int>(table.table.size());
	binary_out.write(reinterpret_cast<char*>(&size), sizeof(unsigned int));
	for (auto const &arg_value_pair : table.table) {
		binary_out.write(reinterpret_cast<char const*>(std::addressof(arg_value_pair.first)), sizeof(ArgT));
		binary_out.write(reinterpret_cast<char const*>(std::addressof(arg_value_pair.second)), sizeof(ValT));
	}
	return binary_out << std::flush;
}

/**
 * Function writes an instance of the StepArgumentTable struct to a binary stream
 * The first several bytes encode: the number of the values to write, the initial argument and the step of the argument
 * The data are written as represented in memory
 */
template <typename ArgT, typename ValT>
std::ostream& write_table_binary(StepArgumentTable<ArgT, ValT> const &table, std::ostream &binary_out) {
	unsigned int size = static_cast<unsigned int>(table.table.size());
	binary_out.write(reinterpret_cast<char*>(&size), sizeof(unsigned int));
	binary_out.write(reinterpret_cast<char const*>(std::addressof(table.arg0)), sizeof(ArgT));
	binary_out.write(reinterpret_cast<char const*>(std::addressof(table.darg)), sizeof(ArgT));
	for (auto const &value : table.table)
		binary_out.write(reinterpret_cast<char*>(&value), sizeof(ValT));
	return binary_out << std::flush;
}

/**
 * Input functions
 */
/**
 * Function reads data from an ascii stream to an instance of the ArgValueTable struct
 * The data are expected to be written by the previous call of the 'write_table_ascii'
 */
template <typename ArgT, typename ValT>
std::istream& read_table_ascii(ArgValueTable<ArgT, ValT> &table, std::istream &ascii_in) {
	unsigned int size;
	ascii_in >> size;
	table.table.resize(size);
	for (auto &arg_value_pair : table.table)
		ascii_in >> arg_value_pair.first >> arg_value_pair.second;
	return ascii_in;
}

/**
 * Function reads data from an ascii stream to an instance of the ArgValueTable struct
 * The data are expected to be written by the previous call of the 'write_table_ascii'
 */
template <typename ArgT, typename ValT>
std::istream& read_table_ascii(StepArgumentTable<ArgT, ValT> &table, std::istream& ascii_in) {
	unsigned int size;
	ascii_in >> size;
	table.table.resize(size);
	for (auto &value : table.table)
		ascii_in >> value;
	return ascii_in;
}

/**
 * Function reads data from a binary stream to an instance of the ArgValueTable struct
 * The data are expected to be written by the previous call of the 'write_table_binary'
 */
template <typename ArgT, typename ValT>
std::istream& read_table_binary(ArgValueTable<ArgT, ValT>& table, std::istream& ascii_in) {
	unsigned int size;
	ascii_in.read(reinterpret_cast<char*>(&size), sizeof(unsigned int));
	table.table.resize(size);
	for (auto& arg_value_pair : table.table) {
		ascii_in.read(reinterpret_cast<char*>(std::addressof(arg_value_pair.first)), sizeof(ArgT));
		ascii_in.read(reinterpret_cast<char*>(std::addressof(arg_value_pair.second)), sizeof(ValT));
	}
	return ascii_in;
}

/**
 * Function reads data from a binary stream to an instance of the StepArgumentTable struct
 * The data are expected to be written by the previous call of the 'write_table_binary'
 */
template <typename ArgT, typename ValT>
std::istream& read_table_binary(StepArgumentTable<ArgT, ValT> &table, std::istream &ascii_in) {
	unsigned int size;
	ascii_in.read(reinterpret_cast<char*>(&size), sizeof(unsigned int));
	table.table.resize(size);
	ascii_in.read(reinterpret_cast<char*>(std::addressof(table.arg0)), sizeof(ArgT));
	ascii_in.read(reinterpret_cast<char*>(std::addressof(table.darg)), sizeof(ArgT));
	for (auto &value : table.table)
		ascii_in.read(reinterpret_cast<char*>(&value), sizeof(ValT));
	return ascii_in;
}

#endif