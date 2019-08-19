#pragma once
#ifndef VDFTabulator_H
#define VDFTabulator_H

#include "AdvancedTable.h"
#include "PhysicalParameters.h"

#include <cmath>

template <typename T>
class InitialVDF final {
public:
	InitialVDF(PhysicalParameters<T> p, T TcTh_ratio, Grid<T> gridVperp, Grid<T> gridVparall): p(p), TcTh_ratio(TcTh_ratio), gridVperp(gridVperp), gridVparall(gridVparall) { }

	MultiscalarTable<T, T> as_multiscalar() const {
		MultiscalarTable<T, T> vdf_table;
		vdf_table.arg_size = 2u; vdf_table.val_size = 1u;
		vdf_table.arguments.resize(vdf_table.arg_size * gridVperp.size * gridVparall.size);
		vdf_table.values.resize(gridVperp.size * gridVparall.size);
		{
			T vperp = gridVperp.begin, vparall = gridVparall.begin;
			for (unsigned step_counter_vparall = 0u; step_counter_vparall < gridVparall.size; ++step_counter_vparall) {
				vparall = gridVparall.begin + gridVparall.step * step_counter_vparall;
				for (unsigned step_counter_vperp = 0u;  step_counter_vperp < gridVperp.size; ++step_counter_vperp) {
					vperp = gridVperp.begin + gridVperp.step * step_counter_vperp;
					unsigned idx = step_counter_vparall * gridVperp.size + step_counter_vperp;
					vdf_table.arguments[2*idx] = vperp;
					vdf_table.arguments[2*idx + 1] = vparall;
					vdf_table.values[idx] = vdf(vperp, vparall);
				}
			}
		}
		return vdf_table;
	}

private:
	T vdf(T vperp, T vparall) const {
		T coeff_c = std::exp(-vperp * vperp * T(0.5)), coeff_h = std::exp(-vperp * vperp * T(0.5)*TcTh_ratio );
		return p.nc * coeff_c * std::exp(-T(0.5) * (vparall - p.bulk_to_term_c) * (vparall - p.bulk_to_term_c)) + p.nh * coeff_h * std::exp(-T(0.5) * (std::sqrt(TcTh_ratio)*vparall - p.bulk_to_term_h) * (std::sqrt(TcTh_ratio) * vparall - p.bulk_to_term_h));
	}

	PhysicalParameters<T> p;
	T TcTh_ratio;
	Grid<T> gridVperp, gridVparall;
};

template <typename T>
auto make_initialvdf(PhysicalParameters<T> p, T TcTh_ratio, Grid<T> gridVperp, Grid<T> gridVparall) {
	return InitialVDF(p,TcTh_ratio,gridVperp,gridVparall);
}

#endif /*VDFTabulator_H*/