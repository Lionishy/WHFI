#include "Table.h"
#include "TableIO.h"

#include <sstream>

int main() {
	using namespace std;
	stringstream ascii_stream;
	ArgValueTable<float,float> table;
	table.table = { {0.f,1.f}, {1.f,2.f}, {2.f, 3.f}, {3.f, 4.f} };
	write_table_ascii(table, ascii_stream);
	cout << ascii_stream.str() << endl;


	return 0;
}