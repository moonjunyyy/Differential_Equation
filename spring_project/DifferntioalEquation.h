#pragma once
#include <iostream>
#include <cmath>
#include <functional>
#include <Eigen/Dense>

class DifferntioalEquation
{
public:
	double x_0, x_1, x_2;
	double dt;
	double b_0, b_1;

	DifferntioalEquation();
	double init();
	virtual void analyticEq(std::ostream& os, double endTime);
	virtual void numericEq(std::ostream& os, double endTime);

	virtual double update();
};

class SpringDE : DifferntioalEquation
{
public:
	double m, b, k;
	SpringDE();
	void analyticEq(std::ostream& os, double endTime);
	void numericEq(std::ostream& os, double endTime);

	double update();
};
