#pragma once
#ifndef Runge_H
#define Runge_H

/**
 * Class Runge4thOrder represents Runge-Kutta step procedure of the 4th order accuracy to solve first order ODE
 * The ODE should be written as dy/dx = F(x,y)
 * F is a function of argument 'x' and the value 'y' of the function under consideration
 */
template <typename ArgT, typename ValT, typename FunT>
class Runge4thOrder final {
public:
	Runge4thOrder(FunT F): F(F) { }

	ValT operator() (ArgT step, ArgT arg, ValT value) const {
		ValT k1 = F(arg, value);
		ValT k2 = F(arg + step / ValT(2.), value + k1 * step / ValT(2.));
		ValT k3 = F(arg + step / ValT(2.), value + k2 * step / ValT(2.));
		ValT k4 = F(arg + step, value + k2 * step); 
		return step / ValT(6.) * (k1 + ValT(2.)*k2 + ValT(2.)*k3 + k4);
	}

private:
	FunT F;
};

/**
 * Function to instantiate Runge4thOrder object with two parameters only: argument type and value type
 */
template <typename ArgT, typename ValT, typename FunT>
auto make_runge4thorder(FunT F) { return Runge4thOrder<ArgT, ValT, FunT>(F); }

#endif