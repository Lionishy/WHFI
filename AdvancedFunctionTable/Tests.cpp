#include "AdvancedTable.h"
#include "AdvancedTableIO.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>

MultiscalarTable<float, float> one_to_one_float_multiscalar_table_prepeare() {
	using namespace std;
	MultiscalarTable<float, float> table;
	table.arg_size = 1u; table.val_size = 1u;
	table.arguments.resize(10); table.values.resize(10);

	float arg0 = 0.f, step = 1.e-1f; unsigned counter = 0;
	generate(begin(table.arguments), end(table.arguments), [arg0, step, &counter] { return arg0 + step * counter++; });
	transform(begin(table.arguments), end(table.arguments), begin(table.values), [](auto arg) {return 10.f * arg; });
	return table;
}

bool one_to_one_multiscalar_table_ascii_write() {
	using namespace std;
	auto table = one_to_one_float_multiscalar_table_prepeare();
	stringstream ascii_stream; ascii_stream << setprecision(3) << fixed;
	write_table_ascii(table, ascii_stream);
	string expected(
"10 1 10 1\n\
0.000\n\
0.100\n\
0.200\n\
0.300\n\
0.400\n\
0.500\n\
0.600\n\
0.700\n\
0.800\n\
0.900\n\
0.000\n\
1.000\n\
2.000\n\
3.000\n\
4.000\n\
5.000\n\
6.000\n\
7.000\n\
8.000\n\
9.000\n\
");
	if (expected != ascii_stream.str()) {
		cerr << "Test 'one_to_one_multiscalar_table_ascii_write' failed!" << endl;
		cerr << "Expected: " << expected << endl;
		cerr << "Actual: " << ascii_stream.str() << endl;
	}
	return expected == ascii_stream.str();
}