#include "ZFunctionWithAsymptotic.h"
#include "Table.h"
#include "TableIO.h"
#include "LambdaR.h"
#include "LambdaRRootDerivative.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include <future>
#include <iomanip>

template <typename T>
std::vector<T> generate_grid(T begin, T end, T step) {
	using namespace std;
	vector<T> grid; grid.reserve(unsigned((end - begin) / step));
	unsigned step_counter = 0; auto arg = begin;
	while (arg <= end) {
		if (abs(arg) > T(1.e-7)) {
			grid.push_back(arg);
		}
		arg = step * ++step_counter + begin;
	}
	return grid;
}

template <typename T, typename FunT>
T step_solver(FunT F, T from, T to, T step) {
	unsigned step_count = 0; T curr = F(from), next = curr;
	do {
		next = F(step * ++step_count);
		if (curr * next <= T(0.))
			return T(0.5) * (step * (step_count - 1u) + step * step_count);
		curr = next;
	} while (step * step_count < to);
	return T(std::nan(""));
}

template <typename T>
ArgValueTable<T, T> real_part_wdr(std::string ZFuncTableFilename, PhysicalParameters<T> params, std::vector<T> k_grid = generate_grid<T>(T(-0.7),T(0.7),T(1.e-3))) {
	using namespace std;
	auto table_ptr = make_shared<StepArgumentTable<T,T>>();
	{
		ifstream zfunc_table_in;
		zfunc_table_in.exceptions(ios::badbit | ios::failbit);

		zfunc_table_in.open(ZFuncTableFilename, ios::in | ios::binary);
		read_table_binary(*table_ptr, zfunc_table_in);
	}
	
	vector<pair<T,future<T>>> omegas;
	auto kernel = [table_ptr, params](T k) {
		return step_solver(make_lambdar(ZFuncWithAsymptotic<T>(table_ptr), params, k), T(0.), T(1.), T(1.e-4));
	};
	for (auto k : k_grid)
		omegas.push_back({ k, async(launch::async,kernel,k) });

	ArgValueTable<T, T> result_table;
	for (auto& f : omegas)
		result_table.table.push_back({ f.first,f.second.get() });
	return result_table;
}

template <typename T>
ArgValueTable<T, T> omega_to_derivative_transform(std::string ZFuncTableFilename, ArgValueTable<T, T> const &k_omega_table, PhysicalParameters<T> p) {
	using namespace std;
	ArgValueTable<T, T> k_derivative;
	k_derivative.table.resize(k_omega_table.table.size());

	auto table_ptr = make_shared<StepArgumentTable<T, T>>();
	{
		ifstream zfunc_table_in;
		zfunc_table_in.exceptions(ios::badbit | ios::failbit);

		zfunc_table_in.open(ZFuncTableFilename, ios::in | ios::binary);
		read_table_binary(*table_ptr, zfunc_table_in);
	}

	ZFuncWithAsymptotic<T> Z(table_ptr);
	transform(begin(k_omega_table.table), end(k_omega_table.table), begin(k_derivative.table), 
		[Z,p](auto& k_omega) {
			return pair<T,T>( k_omega.first, make_lambdar_rootderivative(Z,p,k_omega.first)(k_omega.second) );
		}
	);
	return k_derivative;
}

/**
 * Function to ease construction of the PhysicalParameters struct from a lesser number of the input values
 * @param nc - share of the core particles
 * @param betta_c - core part betta (pressure ratio)
 * @param TcTh_ratio - ratio of the core temperature to the halo temperature
 * @param bulk_velocity - core bulk velocity in terms of Alfven speed
 */
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

int main() {
	using namespace std;
	try {
		auto p = calculate_parameters(0.85f, 1.f / 0.85f, 0.1f, -1.f);
		auto k_omega_table = real_part_wdr("./fZFunc.tbl", p);
		auto k_derivative_table = omega_to_derivative_transform("./fZFunc.tbl", k_omega_table, p);
		{
			ofstream wdr("./fwdr-1.txt"); wdr << setprecision(8) << fixed;
			write_table_ascii(k_omega_table, wdr);
		}
		{
			ofstream wdr_derivative("./fwdr_derive-1.txt"); wdr_derivative << setprecision(8) << fixed;
			write_table_ascii(k_derivative_table, wdr_derivative);
		}
	}
	catch (exception const& ex) {
		cout << ex.what() << endl;
	}



	return 0;
}