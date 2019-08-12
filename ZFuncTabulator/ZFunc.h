#pragma once
#ifndef ZFunc_H
#define ZFunc_H

#include "Tabulator.h"
#include "Runge.h"
#include "Table.h"

#include <vector>
#include <utility>
#include <algorithm>

/**
 * Function creates a vector of values of the plasma dispersion function dZ/dx = -x*Z - 1 for the arguments from 0 to the max_argument with the darg step
 * darg and max_argument are expected to be positive
 */
template <typename T>
std::vector<T> ZFunc_tabulator(T darg, T max_argument, unsigned loop_count = 1u) {
	using namespace std;
	unsigned step_count = unsigned(max_argument / darg)/loop_count + 1u;
	vector<T> table(step_count);
	table[0] = T(0);
	kahan_tabulator<T,T>(
		make_runge4thorder<T,T>([](T arg, T val)->T { return -arg * val - T(1.); })
		, begin(table)+1, end(table)
		, T(0), T(0), darg
		, loop_count
	);
	return table; //we count on NRVO
}

/**
 * Function calls 'ZFunc_tabulator' and transforms the values vector to the AgrValueTable structure
 * darg and max_argument follow the 'ZFunc_tabulator' formal parameters constrains
 */
template <typename T>
ArgValueTable<T,T> ZFuncArgValueTable(T darg, T max_argument, unsigned loop_count = 1u) {
	using namespace std;
	vector<T> values = ZFunc_tabulator<T>(darg, max_argument, loop_count);
	T write_darg = darg * loop_count; unsigned arg_counter = 0;
	vector<pair<T,T>> table(values.size());
	transform(begin(values), end(values), begin(table)
		, [&arg_counter, write_darg](T value) { return make_pair(write_darg*(arg_counter++), value);}
	);
	return { table };
}

/**
 * Function calls 'ZFunc_tabulator' and creates an instance of the StepArgumentTable structure
 * darg and max_argument follow the 'ZFunc_tabulator' formal parameters constrains
 */
template <typename T>
StepArgumentTable<T,T> ZFuncStepArgumentTable(T darg, T max_argument, unsigned loop_count = 1u) {
	using namespace std;
	vector<T> values = ZFunc_tabulator<T>(darg, max_argument, loop_count);
	return { T(0), darg, values };
}


#endif