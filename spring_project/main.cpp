#include <iostream>
#include <fstream>
#include <cmath>
#include "DifferntioalEquation.h"

int main(int argc, char* argv[])
{
	std::fstream fio("Simple_harmonic.csv", std::ios::out);

	SpringDE SDE1(0.1, 0.0, 10, 0, -20.);
	if (!fio) { std::cout << "Fail to open harmonic" << std::endl; goto fail_harmonic; }
	SDE1.dt = 0.01; // w = 10, T = 0.6283, T/20 = 0.0314
	SDE1.numericEq();
	SDE1.Output(fio, 3.);

	fio.close();
	
fail_harmonic:

	fio.open("OverDamp.csv", std::ios::out);

	SpringDE SDE2(0.1, 2.5, 10, 0, -15.);
	if (!fio) { std::cout << "Fail to open harmonic" << std::endl; goto fail_overdamp; }
	SDE2.dt = 0.01;
	SDE2.numericEq();
	SDE2.Output(fio, 1.);

	fio.close();

fail_overdamp:

	fio.open("UnderDamp.csv", std::ios::out);

	SpringDE SDE3(0.1, 0.2, 10, 1, 0);
	if (!fio) { std::cout << "Fail to open harmonic" << std::endl; goto fail_underdamp; }
	SDE3.dt = 0.01; // w = 3 sqrt(11) , T = 0.6315, T/20 = 0.0317
	SDE3.numericEq();
	SDE3.Output(fio, 5.);

	fio.close();

fail_underdamp:

	fio.open("RLC.csv", std::ios::out);

	RLCCircuitDE RLCDE(280., 100E-3, 0.4E-6, 48, 0, 0);
	if (!fio) { std::cout << "Fail to open RLC" << std::endl; goto fail_RLC; }
	RLCDE.dt = 0.00005; //
	RLCDE.numericEq();
	RLCDE.Output(fio, 0.01);

	fio.close();

fail_RLC:

	return 0;
}