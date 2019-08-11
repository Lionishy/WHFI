#pragma once
#ifndef Tabulator_H
#define Tabulator_H

/**
 * Function takes advantage of the Kahan accumulation procedure in tabulating the function with any stepping forward method (Euler, Runge-Kutta etc.)
 */
template < typename ArgT, typename ValT, typename Step, typename InputIterator>
void kahan_tabulator(Step step, InputIterator result_begin, InputIterator result_end, ArgT arg0, ArgT f0, ArgT darg, unsigned int loop_size = 1u) {
	ValT f = f0, arg = arg0;
	for (; result_begin != result_end; ++result_begin) { //Kahan summation
		for (size_t loop_counter = 0; loop_counter != loop_size; ++loop_counter) { //In case we want to skip some irrelevant small substeps, 'loop_counter' is present
			volatile ValT y, t;      //Unoptimized temporary values for the Kahan accumulation
			volatile ValT c = { 0 }; //'c' stands for the 'low-bits error'
			//Algebraically, 'c' should always be zero
			//Beware overly-aggressive optimizing compilers! Here's why we are making use of the 'volatile' keyword with 'y', 't' and 'c' been declared

			y = step(darg, arg, f) - c;    //So far, so good: 'c' is zero
			t = f + y;       //Alas, 'sum' is big, 'y' small, so low-order digits of y are lost
			c = (t - f) - y; //(t - sum) recovers the high-order part of 'y'; subtracting 'y' recovers -(low part of 'y')
			f = t;  
			
			arg += darg;
		}
		*result_begin = f;
	}
}

#endif