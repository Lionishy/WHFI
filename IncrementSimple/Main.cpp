#include "ZFunc.h"
#include "PhysicalParameters.h"
#include "Table.h"
#include "WDRExport.h"
#include "WDRImport.h"
#include "VDFExport.h"
#include "VDFImport.h"
#include "KahanVDFMoments.h"
#include "AnalyticalMoments.h"
#include "MathHelper.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <utility>
#include <tuple>
#include <cmath>
#include <optional>

template <typename T>
std::optional<T> resonant_term(iki::UniformGridTable<T,1u,2u> const &moments, T k, T vr, T betta_root_c) {
    if (moments.getGrid()[0].begin + moments.getGrid()[0].step*(moments.getGrid()[0].size-1u) <= vr || moments.getGrid()[0].begin + moments.getGrid()[0].step >= vr) return std::nullopt;
    
    auto idx = uint64_t((vr - moments.getGrid()[0].begin)/moments.getGrid()[0].step);
    T g_interpolation, G_interpolation, v0 = moments.getGrid()[0].begin + moments.getGrid()[0].step*idx;
    {
        auto left = moments[2*idx], right = moments[2*(idx+1)];
        g_interpolation = (right-left)/moments.getGrid()[0].step*(vr-v0) + left;
    }
    
    {
        auto left = T(0.5)*(moments[2*(idx+1)+1u] - moments[2*(idx-1)+1u])/moments.getGrid()[0].step, right = T(0.5)*(moments[2*(idx+2)+1u] - moments[2*idx+1u])/moments.getGrid()[0].step;
        G_interpolation = (right-left)/moments.getGrid()[0].step*(vr-v0) + left;
    }
    
    return G_interpolation*iki::sign(k) - g_interpolation/(std::abs(k)*betta_root_c);
}

int main(int argc, const char * argv[]) {
    using namespace std;
    using namespace iki;
    using namespace iki::whfi;
    
    if (argc < 2) {
        cout << "Empty command" << endl;
        return 0;
    }
    std::string command(argv[1]);
	if ("zfunc" == command) {
		auto table = ZFunc_tabulator(1.e-6f, 20.f, 100);
		{
			ofstream binary_stream("./fZFunc.tbl",ios::binary);
			write_table_binary(binary_stream,table);
			cout << (binary_stream ? "Success!" : "Fail!") << endl;
		}

		return 0;
	}
    
    //fundamental parameters
    float nc, betta_c, TcTh_ratio, bulk_to_alfven;
    if (argc < 6) {
        cout << "Expected four numbers: nc betta_c TcTh_ratio bulk_to_alfven" << endl;
        return 0;
    }
    stringstream params_ios(string(argv[2]) + ' ' + string(argv[3]) + ' ' + string(argv[4]) + ' ' + string(argv[5]));
    params_ios >> nc >> betta_c >> TcTh_ratio >> bulk_to_alfven;
    if (!params_ios) {
        cout << "Can't read parameters" << endl;
        return 0;
    }
    
    PhysicalParameters<float> params(nc,betta_c/nc,TcTh_ratio,bulk_to_alfven);
    stringstream filename_suffix;
    filename_suffix << setprecision(2) << fixed << "-" << nc << "-" << betta_c << "-" << TcTh_ratio << "-" << (bulk_to_alfven < 0.f ? "m" : "") << std::fabs(bulk_to_alfven);
    
    string wdr_filename = string("./fWDR") + filename_suffix.str() + string(".tbl");
    string vdf_filename = string("./fVDF") + filename_suffix.str() + string(".tbl");
    
	try {
		if ("init" == command) {
			WDRExportBinary(params, -0.7f, 0.7f, 1.e-2f, "./fZFunc.tbl", wdr_filename);
			VDFExportUniformBinary(vdf_filename, params, { {64u,0.f,16.e-2f}, {512u,-12.3f,48.f * 1.e-3f} });
		}
		else if ("gamma" == command) {
			auto vdf_table = VDFImportUniformBinary<float>(vdf_filename);
			UniformGridTable<float, 1u, 2u> moments({ {512u,-12.3f,48.f * 1.e-3f} });
			auto it = vdf_table.begin(); auto vperp_size = vdf_table.getGrid()[0].size;
			for (uint64_t vparall_idx = 0u, end = moments.values_count(); vparall_idx != end;) {
				moments[vparall_idx++] = gUniformGrid(it, it + vperp_size, vdf_table.getGrid()[0]);
				moments[vparall_idx++] = GUniformGrid(it, it + vperp_size, vdf_table.getGrid()[0]);
				it += vperp_size;
			}

			auto wdr_table = WDRImportTableBinary<float>(wdr_filename);

			float const PI = 3.1415927f;
			vector<pair<float, float>> gamma;
			for (uint64_t wdr_idx = 0u; wdr_idx != wdr_table.getGrid()[0].size(); ++wdr_idx) {
				auto opt = resonant_term(moments, wdr_table.getGrid()[0][wdr_idx], wdr_table[3 * wdr_idx + 2], params.betta_root_c);
				if (!opt) continue;
				gamma.push_back({ wdr_table.getGrid()[0][wdr_idx],
					-std::sqrt(PI / 2.f) * *opt / wdr_table[3 * wdr_idx + 1]
				});
			}
			{
				ofstream ascii_stream;
				ascii_stream.exceptions(ios::badbit | ios::failbit);

				ascii_stream.open("./gamma" + filename_suffix.str() + ".txt");
				ascii_stream << setprecision(7) << fixed;
				for (auto const &p : gamma) {
					ascii_stream << p.first << ' ' << p.second << '\n';
				}
			}

			/*{
				auto analytical = AnalyticalMoments(params).g(-12.3f, 6.f * 1.e-3f, 4097u);
				vector<tuple<float,float,float>> g_moment;
				auto vperp_size = vdf_table.getGrid()[0].size;
				float arg0 = vdf_table.getGrid()[1].begin, arg_step = vdf_table.getGrid()[1].step;
				uint64_t arg_counter = 0u;
				for (auto it = vdf_table.begin(), end = vdf_table.end(); it != end; it += vperp_size)
					g_moment.push_back({
						arg0 + arg_step*(arg_counter++)
						, gUniformGrid(it,it+vperp_size,vdf_table.getGrid()[0])
						, analytical[arg_counter]
					});
				{
					ofstream ascii_stream("./g_moment.txt");
					ascii_stream << setprecision(7) << fixed;
					for_each(begin(g_moment),end(g_moment),[&ascii_stream](auto &p){
						ascii_stream << get<0>(p) << ' ' << get<1>(p) << ' ' << get<2>(p) << ' ' << std::fabs(get<1>(p)-get<2>(p)) << '\n';
					});
				}
			}

			{
				auto analytical = AnalyticalMoments(params).G(-12.3f, 6.f * 1.e-3f, 4097u);
				vector<tuple<float,float,float>> G_moment;
				auto vperp_size = vdf_table.getGrid()[0].size;
				float arg0 = vdf_table.getGrid()[1].begin, arg_step = vdf_table.getGrid()[1].step;
				uint64_t arg_counter = 0u;
				for (auto it = vdf_table.begin(), end = vdf_table.end(); it != end; it += vperp_size)
					G_moment.push_back({
						arg0 + arg_step*(arg_counter++)
						, GUniformGrid(it,it+vperp_size,vdf_table.getGrid()[0])
						, analytical[arg_counter]
					});
				{
					ofstream ascii_stream("./g2_moment.txt");
					ascii_stream << setprecision(7) << fixed;
					for_each(begin(G_moment),end(G_moment),[&ascii_stream](auto &p){
						ascii_stream << get<0>(p) << ' ' << get<1>(p) << ' ' << get<2>(p) << ' ' << std::fabs(get<1>(p)-get<2>(p)) << '\n';
					});
				}
			}*/
		}
		else {
			cout << "Unknonw command" << endl;
			return 0;
		}
	}
	catch (exception const &ex) {
		cerr << ex.what() << endl;
	}
        
        
    
    
    return 0;
}
