#pragma once

#include <string>
#include <fstream>
#include <complex>

class CUtility
{
public:
	static bool file_exists(std::string filename);

	// Code for calculating Lambert W from https://github.com/IstvanMezo/LambertW-function
	static std::complex<double> LambertW(std::complex<double> z, int k = 0);

private:
	static std::complex<double> zexpz(std::complex<double> z) { return z * exp(z);	}
	//The derivative of z * exp(z) = exp(z) + z * exp(z)
	static std::complex<double> zexpz_d(std::complex<double> z) { return exp(z) + z * exp(z); }
	//The second derivative of z * exp(z) = 2. * exp(z) + z * exp(z)
	static std::complex<double> zexpz_dd(std::complex<double> z)	{ return 2. * exp(z) + z * exp(z); }

	//Determine the initial point for the root finding
	static std::complex<double> InitPoint(std::complex<double> z, int k);
};

