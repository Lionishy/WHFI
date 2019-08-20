#pragma once
#ifndef AnalyticalMoments_H
#define AnalyticalMoments_H

#include "PhysicalParameters.h"

#include <vector>
#include <cmath>

template <unsigned pow, typename T>
auto power(T x) {
	T res = T(1);
	for (unsigned idx = 0u; idx < pow; ++idx) { res *= x; }
	return res;
}

template <typename T>
class AnalyticalMoments final {
public:
	AnalyticalMoments(PhysicalParameters<T> p, T TcTh_ratio): p(p), TcTh_ratio(TcTh_ratio) { }

	std::vector<T> g(T vparall_begin, T vparall_step, unsigned size) const {
		auto g_vparall = [this](T vparall) {
			return p.nc * std::exp(-power<2>(vparall - p.bulk_to_term_c)/T(2)) + p.nh * std::sqrt(TcTh_ratio) * exp(-power<2>(vparall * std::sqrt(TcTh_ratio) - p.bulk_to_term_h)/T(2));
		};
		std::vector<T> table(size);
		for (unsigned count = 0u; count < size; ++count) 
			table[count] = g_vparall(vparall_begin + vparall_step * count);
		return table;
	}

	std::vector<T> G(T vparall_begin, T vparall_step, unsigned size) const {
		auto G_vparall = [this](T vparall) {
			return p.nc * std::exp(-power<2>(vparall - p.bulk_to_term_c) / T(2)) + p.nh * std::sqrt(T(1)/TcTh_ratio) * exp(-power<2>(vparall * std::sqrt(TcTh_ratio) - p.bulk_to_term_h) / T(2));
		};
		std::vector<T> table(size);
		for (unsigned count = 0u; count < size; ++count)
			table[count] = G_vparall(vparall_begin + vparall_step * count);
		return table;
	}

private:
	PhysicalParameters<T> p;
	T TcTh_ratio;
};

#endif /*AnalyticalMoments_H*/