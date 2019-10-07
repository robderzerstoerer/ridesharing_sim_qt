#include <vector>
#include <complex>
#include <qmessagebox.h>
#include "utility.h"

#ifndef _EPSILON
#define MACRO_EPSILON 0.000001
#define _EPSILON
#endif

bool CUtility::file_exists(std::string filename)
{
	std::ifstream infile(filename);
	return infile.good();
}

std::complex<double> CUtility::InitPoint(std::complex<double> z, int k)
{
	const double pi{ 3.14159265358979323846 };
	const double e{ 2.71828182845904523536 };
	std::complex<double> I{ 0, 1 };
	std::complex<double> two_pi_k_I{ 0., 2. * pi * k };
	std::complex<double> ip{ log(z) + two_pi_k_I - log(log(z) + two_pi_k_I) };// initial point coming from the general asymptotic approximation
	std::complex<double> p{ sqrt(2. * (e * z + 1.)) };// used when we are close to the branch cut around zero and when k=0,-1

	if (abs(z - (-exp(-1.))) <= 1.) //we are close to the branch cut, the initial point must be chosen carefully
	{
		if (k == 0) ip = -1. + p - 1. / 3. * pow(p, 2) + 11. / 72. * pow(p, 3);
		if (k == 1 && z.imag() < 0.) ip = -1. - p - 1. / 3. * pow(p, 2) - 11. / 72. * pow(p, 3);
		if (k == -1 && z.imag() > 0.) ip = -1. - p - 1. / 3. * pow(p, 2) - 11. / 72. * pow(p, 3);
	}

	if (k == 0 && abs(z - .5) <= .5) ip = (0.35173371 * (0.1237166 + 7.061302897 * z)) / (2. + 0.827184 * (1. + 2. * z));// (1,1) Pade approximant for W(0,a)

	if (k == -1 && abs(z - .5) <= .5) ip = -(((2.2591588985 +
		4.22096 * I) * ((-14.073271 - 33.767687754 * I) * z - (12.7127 -
			19.071643 * I) * (1. + 2. * z))) / (2. - (17.23103 - 10.629721 * I) * (1. + 2. * z)));// (1,1) Pade approximant for W(-1,a)

	return ip;
}

std::complex<double> CUtility::LambertW(std::complex<double> z, int k)
{
	//For some particular z and k W(z,k) has simple value:
	if (z == 0.) return (k == 0) ? 0. : -INFINITY;
	if (z == -exp(-1.) && (k == 0 || k == -1)) return -1.;
	if (z == exp(1.) && k == 0) return 1.;

	//Halley method begins
	std::complex<double> w{ InitPoint(z, k) }, wprev{ InitPoint(z, k) }; // intermediate values in the Halley method
	const unsigned int maxiter = 30; // max number of iterations. This eliminates improbable infinite loops
	unsigned int iter = 0; // iteration counter
	double prec = 1.E-30; // difference threshold between the last two iteration results (or the iter number of iterations is taken)

	do
	{
		wprev = w;
		w -= 2. * ((zexpz(w) - z) * zexpz_d(w)) /
			(2. * pow(zexpz_d(w), 2) - (zexpz(w) - z) * zexpz_dd(w));
		iter++;
	} while ((abs(w - wprev) > prec) && iter < maxiter);
	return w;
}

double CUtility::two_node_av_scheduled_customers(ULL cap, double x, ULL B, double doubcap)
{
	double dcap = (double)cap;
	if (doubcap < 0.0)
		doubcap = dcap;
	if ((doubcap - x) <= MACRO_EPSILON)
		return 1e10;

	std::vector<std::complex<double>> zeroes;
	if (cap % 2 == 1)
	{
		for (int i = 0; i < cap; i++)
		{
			std::complex<double> arg = std::polar<double>(1.0, pi * ((double)i) / (dcap)) * x / (dcap)* exp(-x / (dcap));
			if (i % 2 == 0)
				arg = -arg;
			zeroes.push_back(-(dcap) / x * LambertW(arg, 0));
		}
	}
	else
	{
		for (int i = 0; i < (dcap / 2 - MACRO_EPSILON); i++)
		{
			std::complex<double> arg = std::polar<double>(1.0, pi * (2 * (double)i) / (dcap)) * x / (dcap)* exp(-x / (dcap));
			zeroes.push_back(-(dcap) / x * LambertW(-arg, 0));
			zeroes.push_back(-(dcap) / x * LambertW(arg, 0));
		}
	}

	std::complex<double> sum{ 0.0, 0.0 };
	for (int i = 0; i < cap; i++)
	{
		if (!(abs(zeroes[i].real() - 1.0) <= MACRO_EPSILON && abs(zeroes[i].imag()) <= MACRO_EPSILON))
			sum += std::complex<double>(1.0, 0.0) / (std::complex<double>(1.0, 0.0) - zeroes[i]);
	}

	if (abs(sum.imag()) >= MACRO_EPSILON)
	{
		QMessageBox Msgbox(QMessageBox::Icon::Critical, "Calculation error", "Calculation of two-node solution was unsuccessful.");
		Msgbox.exec();
	}

	return 2 / B * ((doubcap - (doubcap - x)* (doubcap - x))/ (2 * (doubcap - x)) + sum.real()) + x;
	//return sum.real();
}

double CUtility::two_node_stddev_scheduled_customers(ULL cap, double x, ULL B)
{
	double dcap = (double)cap;
	if ((dcap - x) <= MACRO_EPSILON)
		return 1e10;

	std::vector<std::complex<double>> zeroes;
	if (cap % 2 == 1)
	{
		for (int i = 0; i < cap; i++)
		{
			std::complex<double> arg = std::polar<double>(1.0, pi * ((double)i) / (dcap)) * x / (dcap)* exp(-x / (dcap));
			if (i % 2 == 0)
				arg = -arg;
			zeroes.push_back(-(dcap) / x * LambertW(arg, 0));
		}
	}
	else
	{
		for (int i = 0; i < (dcap / 2 - MACRO_EPSILON); i++)
		{
			std::complex<double> arg = std::polar<double>(1.0, pi * (2 * (double)i) / (dcap)) * x / (dcap)* exp(-x / (dcap));
			zeroes.push_back(-(dcap) / x * LambertW(-arg, 0));
			zeroes.push_back(-(dcap) / x * LambertW(arg, 0));
		}
	}


	std::complex<double> sum{ 0.0, 0.0 };
	for (int i = 0; i < cap; i++)
	{
		if (!(abs(zeroes[i].real() - 1.0) <= MACRO_EPSILON && abs(zeroes[i].imag()) <= MACRO_EPSILON))
			sum += zeroes[i] / ((std::complex<double>(1.0, 0.0) - zeroes[i])* (std::complex<double>(1.0, 0.0) - zeroes[i]));
	}

	if (abs(sum.imag()) >= MACRO_EPSILON)
	{
		QMessageBox Msgbox(QMessageBox::Icon::Critical, "Calculation error", "Calculation of two-node solution was unsuccessful.");
		Msgbox.exec();
	}

	return sqrt(2 / B * ((dcap * (dcap + 2 * x) + 6 * x * (dcap - x) * (dcap - x) - (dcap - x) * (dcap - x) * (dcap - x) * (dcap - x)) / (12 * (dcap - x) * (dcap - x)) - sum.real()) + x);
}