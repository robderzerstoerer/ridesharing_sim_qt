#include <vector>
#include <complex>
#include <random>
#include <qmessagebox.h>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
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

std::string CUtility::find_new_filename(std::string filename)
{
	bool b_found_name = false;
	int itry = 2;
	std::string newname;
	while (!b_found_name)
	{
		newname = filename.substr(0, filename.length() - 4) + "_" + std::to_string(itry) + ".dat";
		if (!CUtility::file_exists(newname))
		{
			b_found_name = true;
		}
		itry++;
	}

	return newname;
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

double CUtility::polylog_B(double p_full, int B)
{
	double sum = 0.0;

	for (int i = 1; i <= B; i++)
	{
		sum += sqrt((double)i) * (1 - p_full) * pow(p_full, i - 1);
	}

	return sum;
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

	return 2 / B * ((doubcap - (doubcap - x)* (doubcap - x))/ (2 * (doubcap - x)) + sum.real()) + (double (B-1))/B * x;
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

	return sqrt(2 / B * ((dcap * (dcap + 2 * x) + 6 * x * (dcap - x) * (dcap - x) - (dcap - x) * (dcap - x) * (dcap - x) * (dcap - x)) / (12 * (dcap - x) * (dcap - x)) - sum.real()) + (double(B - 1)) / B * x);
}

bool CUtility::read_file(std::string filename,
	std::map<std::string, std::vector<double>>& averages,
	std::map<std::string, std::vector<double>>& stddevs,
	std::map<std::string, std::vector<double>>& counts)
{
	averages.clear();
	stddevs.clear();
	counts.clear();

	std::string line;
	std::ifstream myfile(filename);
	std::string prevLine = "";
	bool NumberNext = false;
	if (myfile.is_open())
	{
		while (std::getline(myfile, line))
		{
			if (line == "PARAMS" || line == "MEASUREMENTS" || line == "")
			{
				NumberNext = false;
			}
			else if (NumberNext == true)
			{
				// extract n, av, stddev
				int find1 = line.find('\t');
				if (find1 != std::string::npos)
				{
					int find2 = line.find('\t', find1 + 1);
					counts[prevLine].push_back(std::stod(line.substr(0, find1)));
					averages[prevLine].push_back(std::stod(line.substr(find1 + 1, find2 - find1 - 1)));
					stddevs[prevLine].push_back(std::stod(line.substr(find2 + 1, line.length() - find2 - 1)));
				}
				else	// no standard deviation or multiple measurements stored
				{
					averages[prevLine].push_back(std::stod(line));
					stddevs[prevLine].push_back(-1.0);
					counts[prevLine].push_back(-1.0);
				}

				NumberNext = false;
			}
			else // new variable incoming
			{
				prevLine = line;
				NumberNext = true;
			}
		}
		myfile.close();

		return true;
	}
	else
		return false;
}

bool CUtility::read_eff_B_data_from_file(std::string eff_filename, std::vector<std::pair<ULL, double>>& out_vBEdata)
{
	std::ifstream myfile(eff_filename);

	bool NumberNext = false;
	std::string line;
	if (myfile.is_open())
	{
		while (std::getline(myfile, line))
		{
			if (line != "")
			{
				// file structure : B <tab> E
				size_t findtab = line.find('\t');
				ULL B = std::stoi(line.substr(0, findtab));
				double E = std::stod(line.substr(findtab + 1));
				if (out_vBEdata.size() > 0)
				{
					// do not save repeated entries
					if (out_vBEdata.back().first != B)
						out_vBEdata.push_back(std::make_pair(B, E));
				}
				else
					out_vBEdata.push_back(std::make_pair(B, E));
			}
		}
		myfile.close();
		return true;
	}
	else
		return false;
}


double CUtility::find_effective_B(std::vector<std::pair<ULL, double>> &vBEdata, double Eff)
{
	// If the measured efficiency is smaller than any efficiency in vBEdata,
	// return the smallest number of buses in vBEdata
	if (Eff < vBEdata[0].second)
	{
		return vBEdata[0].first;
	}

	// otherwise, we have to interpolate
	for (int i = 0; i < vBEdata.size(); i++)
	{
		if (i < vBEdata.size() - 1)
		{
			if (Eff > vBEdata[i].second && Eff < vBEdata[i + 1].second)
			{
				// linear interpolation
				double B_equiv = ((double)vBEdata[i].first) + ((double)(vBEdata[i + 1].first - vBEdata[i].first)) * (Eff - vBEdata[i].second) / (vBEdata[i + 1].second - vBEdata[i].second);
				return B_equiv;
			}
			else if (Eff < vBEdata[i].second && Eff > vBEdata[i + 1].second)
			{
				return (vBEdata[i].second + vBEdata[i + 1].second) / 2;
			}
		}
		else
		{
			// We have reached the end of our efficiency list. Return the last value (highest B).
			return (double)vBEdata[i].first;
		}
	}
}

bool CUtility::read_B_halves(std::string file, std::map<std::string, std::vector<std::pair<float, float>>> &out)
{
	out.clear();

	std::ifstream myfile(file);

	std::string line;
	std::string curtopo = "";

	if (myfile.is_open())
	{
		while (std::getline(myfile, line))
		{
			if (line.find_first_of("abcdefghijklmnopqrstuvwxyz") != std::string::npos)
			{
				curtopo = line;
			}
			else
			{
				// new (x, B1/2) pair
				size_t findspace = line.find(' ');
				float x = std::stof(line.substr(0, findspace));
				float bh = std::stof(line.substr(findspace + 1));
				out[curtopo].push_back(std::make_pair(x, bh));
			}
		}
		myfile.close();
		return true;
	}
	else
		return false;
}

double CUtility::find_B_half_x(std::vector <std::pair<float, float>>& bhalves, double x)
{
	for (int i = 0; i < bhalves.size(); i++)
	{
		if (bhalves[i].first > x)
		{
			if (i == 0)
				return bhalves[0].second;

			double factor = (x - bhalves[i - 1].first) / (bhalves[i].first - bhalves[i - 1].first);
			return bhalves[i - 1].second + factor * (bhalves[i].second - bhalves[i - 1].second);
		}
	}

	return bhalves[bhalves.size() - 1].second;
}

double CUtility::calc_p_full(int B, int seated_customers, int capacity)
{
	double sum = 0.0;

	for (int n = 1; n <= B; n++)
	{
		if ((seated_customers - n * capacity) < 0)
			break;
		if ((seated_customers - n * capacity) > (B - n) * (capacity - 1))
			continue;

		double bin1 = boost::math::binomial_coefficient<double>(B, n);
		int a = B - n;
		int b = seated_customers - n * capacity;
		int c = capacity - 1;
		double cust_perm = num_customer_permutations(a, b, c);
		sum += n * bin1 * cust_perm;
	}

	return sum / (B * num_customer_permutations(B, seated_customers, capacity));
}

double CUtility::num_customer_permutations(int n, int k, int cap)
{
	boost::multiprecision::cpp_dec_float_100 sum = 0.0;

	for (int j = 0; j <= n; j++)
	{
		boost::multiprecision::cpp_dec_float_100 bin1 = boost::math::binomial_coefficient<boost::multiprecision::cpp_dec_float_100>(n, j);

		boost::multiprecision::cpp_dec_float_100 bin2 = 0.0;
		int up = n + k - j * (cap + 1) - 1;
		if (up >= n - 1)
			bin2 = boost::math::binomial_coefficient<boost::multiprecision::cpp_dec_float_100>(up, n - 1);

		boost::multiprecision::cpp_dec_float_100 term = bin1 * bin2;
		if ((j % 2) == 1)
			term = -term;

		sum += term;
	}

	return sum.convert_to<double>();
}

double CUtility::fit_func_2d(double B_B_half, double p_full_inf)
{
	return B_B_half / (B_B_half + sqrt(B_B_half) / (sqrt(B_B_half) + FIT_A)) *
		(1.0 - 1.0 / B_B_half *
		 ((FIT_C + FIT_D * B_B_half + FIT_E * B_B_half * B_B_half) * (p_full_inf - exp(-1 / p_full_inf)) -
		  (FIT_G + FIT_H * B_B_half + FIT_I * B_B_half * B_B_half) * (p_full_inf * p_full_inf - exp(-1 / p_full_inf)) +
		  (FIT_J + FIT_K * B_B_half + FIT_L * B_B_half * B_B_half) * (sqrt(p_full_inf) - exp(-1 / p_full_inf)) - exp(-1 / p_full_inf)));
}