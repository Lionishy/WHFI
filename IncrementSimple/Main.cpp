#include "AdvancedTable.h"
#include "AdvancedTableIO.h"
#include "PhysicalParameters.h"
#include "AnalyticalMoments.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iterator>
#include <utility>
#include <exception>
#include <stdexcept>
#include <algorithm>
#include <cmath>

template <typename T, typename Iterator>
auto kahan_trapez_integration(Iterator begin, Iterator end, T step) {
	volatile T y = T(0), t = T(0), c = T(0); T f = T(0);
	if (begin == end && next(begin) != end) throw std::invalid_argument("To few values to integrate");
	for (auto prev = begin, curr = next(begin); curr != end; ++prev,++curr) {
		T next_val = T(0.5) * (*prev + *curr) * step;
		y = next_val - c;
		t = f + y; 
		c = (t - f) - y; 
		f = t;
	}
	return f;
}

template <typename T>
class VDF final {
public:
	VDF(PhysicalParameters<T> p, T TcTh_ratio): p(p), TcTh_ratio(TcTh_ratio) { }

	T operator()(T vperp, T vparall) const {
		T coeff_c = std::exp(-vperp * vperp * T(0.5)), coeff_h = std::exp(-vperp * vperp * T(0.5) * TcTh_ratio);
		return p.nc * coeff_c * std::exp(-T(0.5) * (vparall - p.bulk_to_term_c) * (vparall - p.bulk_to_term_c)) + p.nh * ((std::sqrt(TcTh_ratio)) * (std::sqrt(TcTh_ratio)) * (std::sqrt(TcTh_ratio))) * coeff_h * std::exp(-T(0.5) * (vparall * std::sqrt(TcTh_ratio) - p.bulk_to_term_h) * (vparall * std::sqrt(TcTh_ratio) - p.bulk_to_term_h));
	}

	std::vector<T> vparall_slice(T vparall, T vperp_begin, T vperp_step, unsigned size) const {
		std::vector<T> table(size);
		for (unsigned count = 0u; count < size; ++count) {
			table[count] = (vperp_begin + vperp_step * count)*this->operator()(vperp_begin + vperp_step * count, vparall);
		}
		return table;
	}

	std::vector<T> vperp_slice(T vperp, T vparall_begin, T vparall_step, unsigned size) const {
		std::vector<T> table(size);
		for (unsigned count = 0u; count < size; ++count) {
			table[count] = this->operator()(vperp,vparall_begin + vparall_step*count);
		}
		return table;
	}

private:
	PhysicalParameters<T> p;
	T TcTh_ratio;
};

template <typename T>
PhysicalParameters<T> calculate_parameters(T nc, T betta_c, T TcTh_ratio, T bulk_velocity) {
	PhysicalParameters<T> params;
	params.nc = nc; params.nh = T(1.) - nc; //nc + nh = 1.
	params.betta_root_c = sqrt(T(0.5) * betta_c);
	params.betta_root_h = sqrt(T(0.5) * betta_c / TcTh_ratio);
	params.bulk_to_term_c = bulk_velocity / params.betta_root_c * sqrt(T(1. / 1836.));
	params.bulk_to_term_h = -(params.nc / params.nh) * params.bulk_to_term_c * sqrt(TcTh_ratio);
	return params;
}

std::pair<std::vector<float>, std::vector<float>> vdf_perp_moments(MultiscalarTable<float, float> const &vdf_data) {
	using namespace std;
	//I know the dimensions at the moment 200x1000 the step is 1.e-2
	vector<float> g(2000), G(2000), g_current(1000), G_current(1000);
	for (unsigned g_idx = 0; g_idx < 2000; ++g_idx) {
		for (unsigned idx = 0; idx < 1000; ++idx) {
			float vperp, vdf;
			vperp = vdf_data.arguments[2*(1000*g_idx + idx)];
			vdf = vdf_data.values[1000*g_idx + idx];
			g_current[1000 - idx - 1] = vperp * vdf;
			G_current[1000 - idx - 1] = vperp * vperp * vperp * 0.5f * vdf;
		}
		g[g_idx] = kahan_trapez_integration(begin(g_current), end(g_current), 1.e-2f);
		G[g_idx] = kahan_trapez_integration(begin(G_current), end(G_current), 1.e-2f);
	}
	return { g,G };
}

template<typename T>
T resonant_term(std::vector<T> const &g, std::vector<T> const &G, T v0, T step, T vr, T k, T betta_root_c) {
	T max_res = v0 + step * (g.size()-1);
	if (std::abs(vr) >= max_res) return T(std::nan(""));
	unsigned idx = unsigned((vr - v0) / step);
	return (G[idx + 1u] - G[idx]) / step * T(k > 0. ? 1 : -1.) -T(0.5) * (g[idx] + g[idx + 1u]) / (std::abs(k) * betta_root_c);
}

int main() {
	using namespace std;
	auto p = calculate_parameters(0.85f, 1.f / 0.85f, 0.1f, -15.f);
	float const PI = 3.1415927f;

	MultiscalarTable<float, float> wdr_data;
	MultiscalarTable<float, float> vdf_data;
	{
		ifstream wdr_data_stream("./fwdr-15.txt");
		read_table_ascii(wdr_data,wdr_data_stream);
	}
	{
		ifstream vdf_data_stream("./fVDF-15.txt");
		read_table_ascii(vdf_data, vdf_data_stream);
	}
	
	if (false) {
		auto moments = AnalyticalMoments(p, 0.1f);
		auto gMoment = moments.g(-10.f, 1.e-2f, 2000u), GMoment = moments.G(-10.f, 1.e-2f, 2000u);
		auto gG = vdf_perp_moments(vdf_data);

		ofstream out("./gG.txt"); out << setprecision(8) << fixed;
		float arg0 = -10.f; unsigned arg_count = 0u;
		for (auto pos = 0; pos != gG.first.size(); ++pos) {
			out << (arg0 + 1.e-2 * arg_count++) << ' ';
			out << gG.first[pos] << ' ' << gG.second[pos] << ' ' << gMoment[pos] << ' ' << GMoment[pos] << '\n';
		}
	}
	
	if (true) {
		auto gG = vdf_perp_moments(vdf_data);
		std::vector<float> gamma(wdr_data.arguments.size());
		for (unsigned pos = 0u; pos != gamma.size(); ++pos) {
			if (std::isnan(wdr_data.values[3 * pos])) { gamma[pos] = 0.f; continue; }
			gamma[pos] = -std::sqrt(PI / 2.f) * resonant_term(gG.first, gG.second, -10.f, 1.e-2f, wdr_data.values[3 * pos + 2], wdr_data.arguments[pos], p.betta_root_c) / wdr_data.values[3 * pos + 1];
			cout << gamma[pos] << endl;
		}
		{
			ofstream gamma_out("./gamma-15.txt"); gamma_out << setprecision(8) << fixed;
			for (unsigned pos = 0u; pos != gamma.size(); ++pos) {
				gamma_out << wdr_data.arguments[pos] << ' ' << gamma[pos] << '\n';
			}
		}
	}

	if (false) {
		auto moments = AnalyticalMoments(p, 0.1f);
		auto gMoment = moments.g(-10.f, 1.e-2f, 2000u), GMoment = moments.G(-10.f, 1.e-2f, 2000u);

		std::vector<float> gamma(wdr_data.arguments.size());
		for (unsigned pos = 0u; pos != gamma.size(); ++pos) {
			gamma[pos] = -std::sqrt(PI / 2.f) * resonant_term(gMoment, GMoment, -10.f, 1.e-2f, wdr_data.values[3 * pos + 2], wdr_data.arguments[pos], p.betta_root_c) / wdr_data.values[3 * pos + 1];
		}
		{
			ofstream gamma_out("./gamma-7.txt"); gamma_out << setprecision(8) << fixed;
			for (unsigned pos = 0u; pos != gamma.size(); ++pos) {
				gamma_out << wdr_data.arguments[pos] << ' ' << gamma[pos] << '\n';
			}
		}
	}



	return 0;
}