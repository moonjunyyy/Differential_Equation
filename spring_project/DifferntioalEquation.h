#pragma once
#include <iostream>
#include <cmath>
class DifferntioalEquation
{
	double x_0, x_1, x_2;

public:
	DifferntioalEquation();
	double init();
	virtual double analyticEq(std::ostream& os, double endTime);
	virtual double numericEq(std::ostream& os, double endTime);
};

class SpringDE : DifferntioalEquation
{
	double m, b, k;
public:
	SpringDE();
	double analyticEq(std::ostream& os, double endTime);
	double numericEq(std::ostream& os, double endTime);
};
