#pragma once
#ifndef WDRExport_H
#define WDRExport_H

#include "PhysicalParameters.h"
#include "Table.h"
#include "TableIO.h"
#include "ZFunctionWithAsymptotic.h"
#include "LambdaR.h"
#include "LambdaRRootDerivative.h"

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <utility>
#include <future>

namespace iki { namespace whfi {
    template <typename T>
    std::shared_ptr<UniformGridTable<T,1u,1u>> load_ZFuncTable(std::string filename) {
        using namespace std;
        using namespace iki;
		try {
			ifstream binary_stream;
			binary_stream.exceptions(ios::badbit | ios::failbit);
			binary_stream.open(filename, ios::binary);

			UniformGrid<T, 1u> grid;
			for (auto &axis : grid)
				read_axis_binary(binary_stream, axis);

			auto table_ptr = make_shared<UniformGridTable<T, 1u, 1u>>(grid);
			read_table_binary(binary_stream, *table_ptr);
			return table_ptr;
		}
		catch (exception &ex) {
			cerr << "Can't read ZFunc from " << filename << endl;
			throw;
		}
    }
    
    template <typename T>
    ArgumentGridTable<T,1u,3u> WDRExportTable(whfi::PhysicalParameters<T> params, T k_begin, T k_end, T k_step, std::string ZFuncTableFile) {
        using namespace  std;
        using namespace iki;
        
        auto ZFuncTable_ptr = load_ZFuncTable<T>(ZFuncTableFile);
        auto ZFunc = ZFuncWithAsymptotic<T>(ZFuncTable_ptr);
        
        auto lambda_r = whfi::make_lambdar(ZFunc,params);
        auto lambda_r_rootderivative = whfi::make_lambdar_rootderivative(ZFunc, params);
        auto resonant_velocity = [params](T omega, T k) { return (omega-T(1.))/(k*params.betta_root_c); };
        
        auto k_omega_kernell = [&lambda_r,&lambda_r_rootderivative,&resonant_velocity](T k)->optional<array<T,3>> {
            array<T,3> result;
            auto omega = step_solve<T>([&lambda_r,k](T omega) { return lambda_r(omega,k); }, T(1./1836.), T(1.), T(1.e-5));
            if (!omega) return nullopt;
            result[0] = *omega;
            result[1] = lambda_r_rootderivative(result[0],k);
            result[2] = resonant_velocity(result[0],k);
            return result;
        };
        
        vector<T> ks;
        {
            unsigned count = 0u; T arg0 = k_begin, end = k_end, arg = arg0;
            while (arg < end) {
                arg = arg0 + ++count*k_step;
                if (std::fabs(arg) < T(1.e-7)) continue;
                ks.push_back(arg);
            }
        }
        
        vector<pair<T,future<optional<array<T,3>>>>> async_result;
        transform(begin(ks),end(ks),back_inserter(async_result),[&k_omega_kernell](auto k)->pair<T,future<optional<array<T,3>>>> {
            return {k,async(launch::async,k_omega_kernell,k)};
        });
        
        vector<T> k_values; vector<T> values;
        for_each(begin(async_result),end(async_result),[&k_values,&values](auto &p) {
            auto &[k,future_opt] = p;
            auto opt = future_opt.get();
            if (!opt) return;
            k_values.push_back(k);
			for (auto val : *opt) {
				values.push_back(val);
			}
        });
        auto wdr_table = ArgumentGridTable<T,1u,3u>(ArgumentGrid<T,1u>({ArgumentAxis<T>(begin(k_values),end(k_values))}));
        copy(begin(values),end(values),begin(wdr_table));
        return wdr_table;
    }
    
    template <typename T>
    void WDRExportAscii(PhysicalParameters<T> params, T k_begin, T k_end, T k_step, std::string ZFuncTableFile, std::string WDRTableFile) {
        using namespace std;
        auto table = WDRExportTable(params,k_begin,k_end,k_step,ZFuncTableFile);
        {
            ofstream ascii_stream;
            ascii_stream.exceptions(ios::badbit|ios::failbit);
            
            ascii_stream.open(WDRTableFile);
            ascii_stream << setprecision(5) << fixed << table;
        }
    }
    
    template <typename T>
    void WDRExportBinary(PhysicalParameters<T> params, T k_begin, T k_end, T k_step, std::string ZFuncTableFile, std::string WDRTableFile) {
        using namespace std;
        auto table = WDRExportTable(params,k_begin,k_end,k_step,ZFuncTableFile);
        {
            ofstream binary_stream;
            binary_stream.exceptions(ios::badbit|ios::failbit);
            
            binary_stream.open(WDRTableFile,ios::binary);
            write_table_binary(binary_stream, table);
        }
    }
} /* whfi */ } /* iki */

#endif /* WDRExport_H */
