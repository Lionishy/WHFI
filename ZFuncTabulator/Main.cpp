#include "ZFunc.h"
#include "Table.h"
#include "TableIO.h"

#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

int main(int argc, char** argv) {
	using namespace std;

	if (argc < 2) {
		cout << "Too few parameters" << endl;
		return 0;
	}

	string command(argv[1]);
	if ("da" == command) { //double ascii argument value
		double const darg = 1.e-6, max_argument = 10.;
		auto table = ZFuncArgValueTable<double>(darg, max_argument, 1);
		ofstream zfunc_out_stream("./dZFunc.txt");
		zfunc_out_stream << setprecision(8) << fixed;
		write_table_ascii(table, zfunc_out_stream);
	}
	else if ("db" == command) { //double binary step argument 
		double const darg = 1.e-6, max_argument = 10.;
		auto table = ZFuncStepArgumentTable<double>(darg, max_argument, 1);
		ofstream zfunc_out_stream("./dZFunc.tbl", ios::out | ios::binary);
		write_table_binary(table, zfunc_out_stream);
	}
	else if ("fa" == command) { //float ascii argument value
		float const darg = 1.e-6f, max_argument = 10.f;
		auto table = ZFuncArgValueTable<float>(darg, max_argument, 1);
		ofstream zfunc_out_stream("./fZFunc.txt");
		zfunc_out_stream << setprecision(8) << fixed;
		write_table_ascii(table, zfunc_out_stream);
	}
	else if ("fb" == command) { //float binary step argument
		float const darg = 1.e-6f, max_argument = 10.f;
		auto table = ZFuncStepArgumentTable<float>(darg, max_argument, 1);
		ofstream zfunc_out_stream("./fZFunc.tbl", ios::out | ios::binary);
		write_table_binary(table, zfunc_out_stream);
	}
	else {
		cout << "Unknown command" << endl;
	}

	return 0;
}