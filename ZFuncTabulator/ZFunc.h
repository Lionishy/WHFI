#pragma once
#ifndef ZFunc_H
#define ZFunc_H

#include "Tabulator.h"
#include "Runge.h"

#include <vector>

template <typename T>
std::vector<T> ZFunc_tabulator(T darg, unsigned loop_count = 1u) {
	using namespace std;
	vector<T> table;
	kahan_tabulator<T,T>(
		make_runge4thorder<T,T>([](T arg, T val)->T { return -arg * val - 1.; })
		, begin(table), end(table)
		, T(0), T(0), darg
		, loop_count
	);
	return table; //we count on NRVO
}


#endif