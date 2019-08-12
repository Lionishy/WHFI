#pragma once
#ifndef ZFunctionWithAsymptotic_H
#define ZFunctionWithAsymptotic_H

#include "Table.h"

#include <vector>
#include <memory>

template <typename T>
class ZFuncWithAsymptotic final {
public:
	T operator()(T arg) const {
		
	}

private:
	std::shared_ptr<StepArgumentTable<T>> known_values;
};

#endif /*ZFunctionWithAsymptotic_H*/