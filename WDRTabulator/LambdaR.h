#pragma once
#ifndef LambdaR_H
#define LambdaR_H

#include "PhysicalParameters.h"

/**
 * Class is a callable representing the real part of the whistler dispersion relation
 * Equation LambdaR(omega) = 0 describes the dispersion curve
 * To instante template one needs:
 *     FunT - a type of the one argument copyable and callable instance
 *     T - a type of the argument and value of the LambdaR callable
 * To consturct LambdaR one needs:
 *     {FunT} Z - a callable instance representing plasma dispersion function (ZFunc)
 *     {PhysicalParameters<T>} p - a bunch of mystical physical parameters discribing VDF
 *     {T} k - a wave number
 */
template <typename T, typename FunT>
class LambdaR final {
public:
	LambdaR(FunT Z, PhysicalParameters<T> p, T k): Z(Z), p(p), k(k) { }

	T operator()(T omega) const { // 1./(omega^2) is omitted deliberately 
		return (1. / 1836.) + k * k
			- p.nc * (omega / (k * p.betta_root_c) - p.bulk_to_term_c) * Z((omega - 1.) / (k * p.betta_root_c) - p.bulk_to_term_c)
			- p.nh * (omega / (k * p.betta_root_h) - p.bulk_to_term_h) * Z((omega - 1.) / (k * p.betta_root_h) - p.bulk_to_term_h);
	}

private:
	FunT Z;
	PhysicalParameters<T> p;
	T k;
};

/**
 * Function to construct LambdaR with the FunT parameter deduced automaticly 
 */
template <typename T, typename FunT>
auto make_lambdar(FunT Z, PhysicalParameters<T> p, T k) {
	 return LambdaR<T, FunT>(Z, p, k);
}

#endif /*LambdaR_H*/