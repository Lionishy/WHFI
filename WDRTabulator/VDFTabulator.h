#pragma once
#ifndef VDFTabulator_H
#define VDFTabulator_H

#include "Table.h"
#include "PhysicalParameters.h"
#include <cmath>

template <typename T>
class InitialVDF final {
public:
	InitialVDF(PhysicalParameters<T> p): p(p) { }
	ArgValueTable<T,T> slice_vperp(T vperp) const {
		T coeff_c = std::exp(-vperp * vperp * T(0.5)), coeff_h = std::exp(-vperp * vperp * T(0.5) * T(0.1));
		auto f_vparall = [coeff_c, coeff_h, this](T vparall) {
			return p.nc * coeff_c * std::exp(-T(0.5) * (vparall - p.bulk_to_term_c) * (vparall - p.bulk_to_term_c)) + p.nh * coeff_h * std::exp(-T(0.5) * (T(std::sqrt(0.1)) * vparall - p.bulk_to_term_h) * (T(std::sqrt(0.1)) * vparall - p.bulk_to_term_h));
		};

		ArgValueTable<T, T> table;
		T arg0 = T(-10.), step = T(1.e-2), end = T(10.), arg = arg0;
		unsigned step_counter = 0;
		do {
			table.table.push_back({ arg, f_vparall(arg) });
			arg = arg0 + step * ++step_counter;
		} while (arg < end);
		return table;
	}

private:
	PhysicalParameters<T> p;
};

template <typename T>
auto make_initialvdf(PhysicalParameters<T> p) {
	return InitialVDF(p);
}

#endif /*VDFTabulator_H*/