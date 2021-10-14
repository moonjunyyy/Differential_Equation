#pragma once
#include <iostream>
#include <cmath>
#include <functional>
#include <Eigen/Dense>

class DifferntioalEquation
{
public:
	double dt = 0.01;
	std::function<double(double)> numeric_sol;

	DifferntioalEquation() { }

	void init();
	virtual void numericEq(){}
	virtual double update() { return 0.; }
	virtual void Output(std::ostream os, double endt) {}
};

class SpringDE : public DifferntioalEquation
{
public:
	double x_p, x_0, x_n;
	double b_0, v_0;
	double m, b, k;
	double D1 = 0, D2 = 0, C1 = 0, C2 = 0;

	SpringDE(double m, double b, double k, double n_1, double n);

	void init();
	void numericEq();
	double update();
	void Output(std::ostream& os, double endt);
};

class RLCCircuitDE : public DifferntioalEquation
{
public:
	double x_p, x_0, x_n;
	double v_0, dv_0;
	double a, b, c, v;
	double D1 = 0, D2 = 0, C1 = 0, C2 = 0;

	RLCCircuitDE(double R, double L, double C, double V,double v_0, double vprime_0);

	void init();
	void numericEq();
	double update();
	void Output(std::ostream& os, double endt);
};