#ifndef MEASUREMENT_COLLECTOR_H
#define MEASUREMENT_COLLECTOR_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#ifndef _INTEGER_TYPES
#define ULL uint64_t
#define LL int64_t
#define _INTEGER_TYPES
#endif

#ifndef _EPSILON
#define MACRO_EPSILON 0.000000000001
#define _EPSILON
#endif

#include "customer.h"
#include "transporter.h"
#include "traffic_network.h"

class customer;
class transporter;

class measure {

public:

	measure();

	void new_measurement(long double val, double time = -1.0);

	/* We can enable calculating averages of this measure for all measurements within a
	certain time interval (partial_average_time_interval) */
	void enable_partial_average_calculation(double time_interval);
	void disable_partial_average_calculation() { calc_partial_averages = false; }

	unsigned long get_num_measurements() { return n; }
	long double get_average() { return av; }
	long double get_stddev() { return stddev; }
	std::vector<double>& get_partial_averages() { return partial_averages; }
	void reset();
	void print(std::ofstream& out);

private:

	unsigned long n;
	long double av;
	long double stddev;

	// Parameters which can be used if averages of this variable for every interval 
	// of 'partial_average_time_interval' time units are needed
	bool calc_partial_averages = false;
	double partial_average_time_interval;
	std::vector<double> partial_averages;
	int last_iTime;
	long double cur_average;
	int cur_n;
};

class measurement_collector
{
public:
	measurement_collector(traffic_network& param_network);
	virtual ~measurement_collector();

	void measure_request(customer& c, transporter& t);
	void measure_trip(std::pair<ULL, double> last_stop, std::pair<ULL, double> current_stop);
	void measure_system_status(std::vector<transporter>& transporter_list, double time);
	void reset();
	void print(std::ofstream& out, bool readable = false);

	measure drive_time;
	measure wait_time;
	measure delay_time;
	measure fraction_of_delayed_trips;

	measure drive_time_between_stops;
	measure drive_distance_between_stops;

	measure occupancy;
	measure scheduled_customers;
	measure number_of_planned_stops;
	measure planned_time_horizon;
	measure number_of_idle_transporters;

	measure p_full;
	measure p_full2; 

private:

	traffic_network& network;

	// Variables for measuring the number of scheduled customers in equidistant time intervals
	// (for verifying the validity of the final result)
	
};

#endif // MEASUREMENTS_H
