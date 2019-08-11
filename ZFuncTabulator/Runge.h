#pragma once
#ifndef Runge_H
#define Runge_H

/**
 * Class Runge4thOrder represents Runge-Kutta step procedure of the 4th order accuracy to solve first order ODE
 * The ODE should be written as dy/dx = F(x,y)
 * F is a function of argument 'x' and the value 'y' of the function under consideration
 */
template <typename FunT, typename ArgT, typename ValT> 
class Runge4thOrder final {
public:
	Runge4thOrder(FunT F): F(F) { }

	ValT operator() (ArgT step, ArgT arg, ValT value) const override {
		ValT k1 = F(arg, value);
		ValT k2 = F(arg + step / 2., value + k1 * step / 2.);
		ValT k3 = F(arg + step / 2., value + k2 * step / 2.);
		ValT k4 = F(arg + step, value + k2 * step); 
		return step / 6. * (k1 + k2 + k3 + k4);
	}

private:
	FunT F;
};

#endif