double runge_step_ZFunc(double darg, double arg, double Z) {
	double k1 = -arg * Z - 1.;
	double k2 = -(arg + darg / 2.) * (Z + darg / 2. * k1) - 1;
	double k3 = -(arg + darg / 2.) * (Z + darg / 2. * k2) - 1;
	double k4 = -(arg + darg) * (Z + darg * k3) - 1;
	return darg / 6. * (k1 + k2 + k3 + k4);
}