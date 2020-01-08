#pragma once

#include <string>
#include <fstream>
#include <complex>

#ifndef _INTEGER_TYPES
#define ULL uint64_t
#define LL int64_t
#define _INTEGER_TYPES
#endif

#ifndef _EPSILON
#define MACRO_EPSILON 0.000000000001
#define _EPSILON
#endif

constexpr double pi = 3.14159265358979323846;

#define FIT_A 0.254783105563
#define FIT_C -1.26347049141
#define FIT_D -0.559174166344
#define FIT_E 0.011734652797
#define FIT_G -10.6369376982
#define FIT_H -1.15375238312
#define FIT_I 0.0217529264488
#define FIT_J 0.471233330997
#define FIT_K 0.183602945701
#define FIT_L -0.00383227915831

class CUtility
{
public:
	static bool file_exists(std::string filename);

	static std::string find_new_filename(std::string filename);

	// Code for calculating Lambert W from https://github.com/IstvanMezo/LambertW-function
	static std::complex<double> LambertW(std::complex<double> z, int k = 0);

	static double polylog_B(double p_full, int B);

	// theoretical solution for minimal (two node) graph and 1 Bus
	static double two_node_av_scheduled_customers(ULL cap, double x, ULL B, double doubcap = -1.0);
	static double two_node_stddev_scheduled_customers(ULL cap, double x, ULL B);

	// reads Efficiency-data from pre-calculated file, which needs do have the two columns
	// "efficiency" and "number_of_buses" respectively. Returns if successful.
	static bool read_eff_B_data_from_file(std::string eff_filename, std::vector<std::pair<ULL, double>> &out_vBEdata);
	// calculates effective B from data of unlimited capacity simulation
	static double find_effective_B(std::vector<std::pair<ULL, double>> &vBEdata, double Eff);

	static bool read_B_halves(std::string file, std::map<std::string, std::vector<std::pair<float, float>>> &out);
	static double find_B_half_x(std::vector <std::pair<float, float>>& bhalves, double x);

	static bool read_file(std::string filename,
		std::map<std::string, std::vector<double>>& averages,
		std::map<std::string, std::vector<double>>& stddevs,
		std::map<std::string, std::vector<double>>& counts);

	static double calc_p_full(int B, int seated_customers, int capacity);

	static double fit_func_2d(double B_B_half, double p_full_inf);


private:
	static std::complex<double> zexpz(std::complex<double> z) { return z * exp(z);	}
	//The derivative of z * exp(z) = exp(z) + z * exp(z)
	static std::complex<double> zexpz_d(std::complex<double> z) { return exp(z) + z * exp(z); }
	//The second derivative of z * exp(z) = 2. * exp(z) + z * exp(z)
	static std::complex<double> zexpz_dd(std::complex<double> z)	{ return 2. * exp(z) + z * exp(z); }

	//Determine the initial point for the root finding
	static std::complex<double> InitPoint(std::complex<double> z, int k);

	static double num_customer_permutations(int n, int k, int cap);
};

