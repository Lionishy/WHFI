#pragma once
#ifndef Runge_H
#define Runge_H

namespace  iki {
    /**
     * Class Runge4thOrder represents Runge-Kutta step procedure of the 4th order accuracy to solve first order ODE
     * The ODE should be written as dy/dx = F(x,y)
     * F is a function of argument 'x' and the value 'y' of the function under consideration
     */
    template <typename T, typename FunT>
    class Runge4thOrder final {
    public:
        Runge4thOrder(FunT F): F(F) { }
        
        T operator() (T step, T arg, T value) const {
            T k1 = F(arg, value);
            T k2 = F(arg + step / T(2.), value + k1 * step / T(2.));
            T k3 = F(arg + step / T(2.), value + k2 * step / T(2.));
            T k4 = F(arg + step, value + k2 * step);
            return step / T(6.) * (k1 + T(2.)*k2 + T(2.)*k3 + k4);
        }
        
    private:
        FunT F;
    };
    
    /**
     * Function to instantiate Runge4thOrder object with two parameters only: argument type and value type
     */
    template <typename T, typename FunT>
    auto make_runge4thorder(FunT F) { return Runge4thOrder<T, FunT>(F); }
} /* iki */

#endif
