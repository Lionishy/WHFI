#include "ZFunctionWithAsymptotic.h"
#include "Table.h"
#include "TableIO.h"

#include <fstream>
#include <iostream>

int main() {
	using namespace std;
	auto fTable_ptr = make_shared<StepArgumentTable<float, float>>();
	{
		ifstream zfunc_table_in;
		zfunc_table_in.exceptions(ios::badbit | ios::failbit);

		zfunc_table_in.open("./fZfunc.tbl",ios::in|ios::binary);
		read_table_binary(*fTable_ptr, zfunc_table_in);
	}
	ZFuncWithAsymptotic<float> fZFunc(fTable_ptr);

	{
		ArgValueTable<float, float> newTable;
		for (float arg = 0.; arg < 100.f; arg += 7.e-3f)
			newTable.table.push_back({ arg,fZFunc(arg) });

		ofstream zfunc_table_out;
		zfunc_table_out.exceptions(ios::badbit | ios::failbit);

		zfunc_table_out.open("./fZfunc.txt");
		write_table_ascii(newTable, zfunc_table_out);
	}


	return 0;
}