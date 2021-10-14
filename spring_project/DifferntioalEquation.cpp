#include "DifferntioalEquation.h"

void DifferntioalEquation::init() { return; }

SpringDE::SpringDE(double m, double b, double k, double n, double v)
{
	this->m = m, this->b = b; this->k = k, this->v_0 = v, this->b_0 = n;
	
	x_p = b_0 - (v_0 * dt);
	x_0 = b_0;
}

void SpringDE::init()
{
	x_p = b_0 - (v_0 * dt);
	x_0 = b_0;
	D1 = D2 = C1 = C2 = 0;
}

void SpringDE::numericEq()
{
	double D = (b * b) - (4 * m * k);
	if (D > 0)
	{
		D = sqrt(D);
		D1 = (-b + D) / (2 * m), D2 = (-b - D) / (2 * m);

		Eigen::MatrixXd A(2, 2), B(2, 1), C;

		A(0, 0) = 1,  A(0, 1) = 1;
		A(1, 0) = D1, A(1, 1) = D2;

		B(0, 0) = b_0, B(1, 0) = v_0;

		C = A.inverse() * B;

		C1 = C(0, 0);
		C2 = C(1, 0);

		numeric_sol = [&](double x) -> double {return C1 * exp(D1 * x) + C2 * exp(D2 * x); };
	}

	else if (D == 0)
	{
		D1 = -b / (2 * m);

		Eigen::MatrixXd A(2, 2), B(2, 1), C;

		A(0, 0) = 1,  A(0, 1) = 0;
		A(1, 0) = D1, A(1, 1) = 1;

		B(0, 0) = b_0, B(1, 0) = v_0;

		C = A.inverse() * B;

		C1 = C(0, 0);
		C2 = C(1, 0);

		numeric_sol = [&](double x) -> double {return C1 * exp(D1 * x) + C2 * x * exp(D1 * x); };
	}
	else
	{
		D2 = sqrt(abs(D)) / (2 * m);
		D1 = -b / (2 * m);

		Eigen::MatrixXd A(2, 2), B(2, 1), C;

		A(0, 0) = 1,  A(0, 1) = 0;
		A(1, 0) = D1, A(1, 1) = D2;

		B(0, 0) = b_0, B(1, 0) = v_0;

		C = A.inverse() * B;

		C1 = C(0, 0);
		C2 = C(1, 0);

		numeric_sol = [&](double x) -> double {return exp(D1 * x) * (C1 * cos(D2 * x) + C2 * sin(D2 * x)); };
	}
	return;
}

double SpringDE::update()
{
	x_n = (dt * dt / m) * (-k * x_0 - b * (x_0 - x_p) / dt) + (2 * x_0 - x_p);
	x_p = x_0; x_0 = x_n;
	return x_0;
}

void SpringDE::Output(std::ostream& os, double endt)
{
	for (double t = 0; t < endt; t += dt)
	{
		update();
		os << t << ", " << x_n << ", " << numeric_sol(t) << std::endl;
	}
}

RLCCircuitDE::RLCCircuitDE(double R, double L, double C, double V, double v_0, double vprime_0)
{
	this->a = 1, this->b = R / L, this->c = 1. / L / C, this->v = V;
	this->v_0 = v_0, this->dv_0 = vprime_0;

	x_0 = this->v_0;
	x_p = this->v_0;
}

void RLCCircuitDE::init()
{
	x_0 = this->v_0;
	x_p = this->v_0;
}

void RLCCircuitDE::numericEq()
{
	double D = (b * b) - (4 * a * c);
	if (D > 0)
	{
		D = sqrt(D);
		D1 = (-b + D) / (2 * a), D2 = (-b - D) / (2 * a);

		Eigen::MatrixXd A(2, 2), B(2, 1), C;

		A(0, 0) = 1, A(0, 1) = 1;
		A(1, 0) = D1, A(1, 1) = D2;

		B(0, 0) = v_0, B(1, 0) = dv_0;

		C = A.inverse() * B;

		C1 = C(0, 0);
		C2 = C(1, 0);

		numeric_sol = [&](double x) -> double {return v + C1 * exp(D1 * x) + C2 * exp(D2 * x); };
	}

	else if (D == 0)
	{
		D1 = -b / (2 * a);

		Eigen::MatrixXd A(2, 2), B(2, 1), C;

		A(0, 0) = 1, A(0, 1) = 0;
		A(1, 0) = D1, A(1, 1) = 1;

		B(0, 0) = v_0, B(1, 0) = dv_0;

		C = A.inverse() * B;

		C1 = C(0, 0);
		C2 = C(1, 0);

		numeric_sol = [&](double x) -> double {return v + C1 * exp(D1 * x) + C2 * x * exp(D1 * x); };
	}
	else
	{
		D2 = sqrt(abs(D)) / (2 * a);
		D1 = -b / (2 * a);

		Eigen::MatrixXd A(2, 2), B(2, 1), C;

		A(0, 0) = 1, A(0, 1) = 0;
		A(1, 0) = D1, A(1, 1) = D2;

		B(0, 0) = -v, B(1, 0) = dv_0;

		C = A.inverse() * B;

		C1 = C(0, 0);
		C2 = C(1, 0);

		numeric_sol = [&](double x) -> double {return v + exp(D1 * x) * (C1 * cos(D2 * x) + C2 * sin(D2 * x)); };
	}
	return;
}

double RLCCircuitDE::update()
{
	x_n = (dt * dt / a) * (v * c + -c * x_0 - b * (x_0 - x_p) / dt) + (2 * x_0 - x_p);
	x_p = x_0; x_0 = x_n;
	return x_0;
}

void RLCCircuitDE::Output(std::ostream& os, double endt)
{
	for (double t = 0; t < endt; t += dt)
	{
		update();
		os << t << ", " << x_n << ", " << numeric_sol(t) << std::endl;
	}
}
