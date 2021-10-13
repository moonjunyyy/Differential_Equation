#include <iostream>
#include <fstream>
#include <cmath>
#include "DifferntioalEquation.h"

int main(int argc, char* argv[])
{
	std::fstream fio("Simple_harmonic.csv", std::ios::out);

	SpringDE SDE1(0.1, 0.0, 10, 0, -20.);
	SDE1.numericEq();
	SDE1.Output(fio, 5.);

	fio.close();
	

	fio.open("OverDamp.csv", std::ios::out);

	SpringDE SDE2(0.1, 2.5, 10, 0, -15.);
	SDE2.numericEq();
	SDE2.Output(fio, 5.);

	fio.close();


	fio.open("UnderDamp.csv", std::ios::out);

	SpringDE SDE3(0.1, 0.2, 10, 1, 0);
	SDE3.numericEq();
	SDE3.Output(fio, 5.);

	fio.close();

	return 0;
}