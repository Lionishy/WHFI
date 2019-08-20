#include "AdvancedTable.h"
#include "AdvancedTableIO.h"

#include <fstream>

MultiscalarTable<float, float> g_integrated(MultiscalarTable<float, float> const &vdf_data) {

}

int main() {
	using namespace std;
	MultiscalarTable<float, float> wdr_data;
	MultiscalarTable<float, float> vdf_data;
	{
		ifstream wdr_data_stream("./fwdr-7.txt");
		read_table_ascii(wdr_data,wdr_data_stream);
	}
	{
		ifstream vdf_data_stream("./fVDF-7.txt");
		read_table_ascii(vdf_data, vdf_data_stream);
	}

	

	return 0;
}