#include "ZFunctionWithAsymptotic.h"
#include "Table.h"
#include "TableIO.h"
#include "LambdaR.h"

#include <fstream>
#include <iostream>
#include <cmath>

int main() {
	using namespace std;


	/*PhysicalParameters<float> params;
	params.nc = 0.85f; params.nh = 0.15f; //nh+nc = 1.
	params.betta_root_c = sqrt(0.5f); params.betta_root_h = sqrt(10.f)*params.betta_root_c; //Tc/Th = 0.1
	params.bulk_to_term_c = -1.f/params.betta_root_c * sqrt(1.f/1836.f);
	params.bulk_to_term_h = (params.nc / params.nh) * params.bulk_to_term_c * sqrt(0.1f); //Tc/Th = 0.1

	auto fTable_ptr = make_shared<StepArgumentTable<float, float>>();
	{
		ifstream zfunc_table_in;
		zfunc_table_in.exceptions(ios::badbit | ios::failbit);

		zfunc_table_in.open("./fZfunc.tbl",ios::in|ios::binary);
		read_table_binary(*fTable_ptr, zfunc_table_in);
	}

	auto LR = make_lambdar(ZFuncWithAsymptotic<float>(fTable_ptr), params, -0.34f);
	unsigned step_count = 0; float step = 1.e-4f, LRcurr = LR(0.f), LRnext = LRcurr;
	do {
		LRnext = LR(step * ++step_count);
		if (LRcurr * LRnext <= 0.f) break;
		LRcurr = LRnext;
	} while (step * step_count < 1.f);
	cout << step * (step_count - 1u) << ' ' << step * step_count << endl;
	*/

	return 0;
}