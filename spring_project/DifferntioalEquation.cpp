#include "DifferntioalEquation.h"

void SpringDE::analyticEq(std::ostream& os, double endtime)
{
	
}

void SpringDE::numericEq(std::ostream& os, double endtime)
{
	std::function<double(double)> Sol;
	double D = (b * b) - (4 * m * k);
	if (D > 0)
	{
		D = sqrt(D);
		double D1 = (-b + D) / (2 * m), D2 = (-b - D) / (2 * m);

		Eigen::MatrixXd A(2, 2), B(2, 1), C;
		A(0, 0) = 1, A(0, 1) = 1, A(1, 0) = D1, A(1, 1) = D2;
		B(0, 0) = b_0, B(1, 0) = b_1;

		C = A.inverse() * B;
		Sol = [&](double x) ->double {C[0, 0] * exp(D1 * x) + C[1, 0] * exp(D2 * x); };
	}

	else if (D = 0)
	{
		D = -b / (2 * m);
		Eigen::MatrixXd A(2, 2), B(2, 1), C;
		A(0, 0) = 1, A(0, 1) = 1, A(1, 0) = D, A(1, 1) = 1;
		B(0, 0) = b_0, B(1, 0) = b_1;

		C = A.inverse() * B;
		Sol = [&](double x) -> double {C[0, 0] * exp(D * x) + C[1, 0] * x * exp(D * x); };
	}
	else
	{
		D = sqrt(abs(D)) / (2 * m);
		double alpha = -b / (2 * m);


	}
}