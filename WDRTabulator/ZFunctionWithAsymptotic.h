#pragma once
#ifndef ZFunctionWithAsymptotic_H
#define ZFunctionWithAsymptotic_H

#include "Table.h"

#include <vector>
#include <memory>
#include <cmath>

/**
 * Class represents plasma dispersion function (ZFunc) as a combination of the known values table and the asymptotic expression for the large values
 */
template <typename T>
class ZFuncWithAsymptotic final {
public:
	ZFuncWithAsymptotic(std::shared_ptr<StepArgumentTable<T,T> const> known_values): known_values(known_values) { }

	T operator()(T arg) const {
		T farg = std::abs(arg);
		unsigned idx = unsigned(farg / known_values->darg);
		if ((idx + 1u) < known_values->table.size()) {
			return T(arg >= 0. ? 1. : -1.) * ((known_values->table[idx + 1u] - known_values->table[idx]) / known_values->darg * (farg - known_values->darg * idx) + known_values->table[idx]);
		}
		else {
			T over = T(1.) / arg, square = over * over;
			return -over * (T(1.) + square + T(3.) * square * square);
		}
	}

private:
	std::shared_ptr<StepArgumentTable<T,T> const> known_values;
};

#endif /*ZFunctionWithAsymptotic_H*/