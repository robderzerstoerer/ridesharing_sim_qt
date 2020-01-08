#ifndef RIDSHARING_SIM_H
#define RIDSHARING_SIM_H

#ifndef _INTEGER_TYPES
#define ULL uint64_t
#define LL int64_t
#define _INTEGER_TYPES
#endif

#ifndef _EPSILON
#define MACRO_EPSILON 0.000000000001
#define _EPSILON
#endif

//class customer;
//class transporter;

#include <queue>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <functional>
#include <ctime>

#include "measurement_collector.h"
#include "traffic_network.h"
#include "customer.h"
#include "transporter.h"


enum ESIM_MODE 
{
	SIMPLE_SIM,
	EVERYTHING,
	B_AIM
};

// Parameters for running many simulations with different parameters
struct program_parameters
{
	ESIM_MODE sim_mode;

	std::vector<std::pair<std::string, ULL>> topology_list;  // (topology_name, number_of_nodes)

	std::vector<ULL> B_list;

	std::vector<int> capacity_list;

	std::vector<double> x_list;

	double E_aim;
	double lambda;

	ULL num_requests_per_bus_init;
	ULL num_requests_per_bus_sim;

	// Do we want to abort any point in the simulation after simulating a certain time?
	double scan_point_sim_time;  // in seconds, ~ 30s
	// when should we abort after simulating a longer time?
	double max_point_sim_time; // in seconds, ~ 5 min

	bool simulate_until_exact;

	bool calc_p_full;
	bool calc_cap_delay;  //unused

	std::string foldername;
	std::string filename;

	bool stop_thread;
};


// parameters for one simulation of definite topology, number of buses, capacity etc.
struct simulation_parameters
{
	program_parameters* parent;

	std::string topology;
	ULL number_of_nodes;

	ULL number_of_buses;

	int capacity;

	double normalized_request_rate;

	bool calc_p_full;
	bool calc_cap_delay;  //unused

	// Do we want to abort this simulation after simulating a certain time?
	double max_sim_time;  // in seconds
};


typedef std::priority_queue< std::pair<double, ULL>, std::vector<std::pair<double, ULL> >, std::greater< std::pair<double, ULL> > > transporter_event_queue_type;

class ridesharing_sim
{
public:
	ridesharing_sim(simulation_parameters* par_par);
	virtual ~ridesharing_sim();

	simulation_parameters* par;

	void set_normalized_request_rate(double param_normalized_request_rate);
	void set_request_rate(double param_request_rate);

	void reset_number_of_buses(ULL param_number_of_buses, int param_capacity);

	void reset_measurements();
	void print_measurements(std::ofstream& out);

	void init_network();

	bool init_new_sim(ULL param_number_of_buses,
		double param_normalized_request_rate,
		LL param_num_init_requests = -1);
	void run_sim_request_list(std::list< std::pair< double, std::pair<ULL, ULL> > > request_list);
	bool run_sim_requests(ULL max_requests);
	void run_sim_time(long double max_time);

	double execute_next_event();

	void print_params(std::ofstream& out, bool readable = false);

	traffic_network network;
	std::vector<transporter> transporter_list;

	double time;
	ULL total_requests;
	ULL total_serviced_requests;

	double start_of_measured_time;
	ULL start_of_measured_total_requests;
	ULL start_of_measured_serviced_requests;

	double start_of_output_time;

	void enable_measurements(double param_measurement_time_step);
	void disable_measurements();
	bool do_measurement;
	double next_measurement_time;
	double measurement_time_step;

	measurement_collector measurements;

	void enable_timeseries_output(double param_output_time_step, std::string output_filename); // bisher nicht aufgerufen
	void disable_timeseries_output();
	bool do_timeseries_output;
	double next_output_time;
	double output_time_step;

	void output(double output_time);
	std::ofstream out;



	transporter_event_queue_type transporter_event_queue;
	double next_request_time;

	double normalized_request_rate;
	double request_rate;

	std::mt19937_64 random_generator;
	std::exponential_distribution<double> exp_dist;

	std::clock_t start_clock;
};

#endif // RIDSHARING_SIM_H
