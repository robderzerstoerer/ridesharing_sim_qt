#ifndef TRANSPORTER_H
#define TRANSPORTER_H

#include <cstdlib>
#include <iostream>
#include <deque>
#include <list>

#ifndef _INTEGER_TYPES
#define ULL uint64_t
#define LL int64_t
#define _INTEGER_TYPES
#endif

#ifndef _EPSILON
#define MACRO_EPSILON 0.000000000001
#define _EPSILON
#endif

#include "measurement_collector.h"
#include "customer.h"
#include "traffic_network.h"

class measurement_collector;
class customer;
class transporter;

void best_offer_st_alone(
	int param_origin,
	int param_destination,
	double param_request_time,
	float delay_delta,
	int capacity,
	float current_time,
	int current_pos,	// actually only needed if bus idle (stops_length == 0), otherwise pass anything to this parameter
	int* stops_nodeindices,
	float* stops_times,
	int* occupancies_before_stops,
	int stops_length,
	float** netdistances,
	float best_dropoff_time,
	float* out_dropoff_time,
	int* out_pickup_insertion,
	int* out_dropoff_insertion);


struct stop
{
	ULL node_index;

	std::list<customer>::iterator c_it;
	double planned_time;

	bool is_pickup;
	bool is_dropoff;

	ULL type;

	stop(ULL param_node_index, std::list<customer>::iterator param_customer, double param_time, bool param_is_pickup, bool param_is_dropoff, ULL param_type)
		: node_index(param_node_index),
		c_it(param_customer),
		planned_time(param_time),
		is_pickup(param_is_pickup),
		is_dropoff(param_is_dropoff),
		type(param_type)
	{};
};

struct offer
{
	ULL transporter_index;
	transporter* best_transporter;

	double pickup_time;
	std::list< stop >::iterator pickup_insertion;

	double dropoff_time;
	std::list< stop >::iterator dropoff_insertion;

	bool is_better_offer;

	offer(ULL param_transporter_index, transporter* param_best_transporter, double param_pickup_time, std::list< stop >::iterator param_pickup_insertion, double param_dropoff_time, std::list< stop >::iterator param_dropoff_insertion)
		: transporter_index(param_transporter_index),
		best_transporter(param_best_transporter),
		pickup_time(param_pickup_time),
		pickup_insertion(param_pickup_insertion),
		dropoff_time(param_dropoff_time),
		dropoff_insertion(param_dropoff_insertion),
		is_better_offer(false)
	{};

	offer() : dropoff_time(std::numeric_limits<ULL>::max()) {};
};

class transporter
{
public:
	transporter(ULL param_index, ULL param_location, int param_capacity, std::mt19937_64 param_random_generator);
	void reset(ULL param_index, ULL param_location, int param_capacity);
	virtual ~transporter();

	ULL get_index() { return(index); }
	LL get_capacity() { return(capacity); }
	void set_capacity(ULL cap) { capacity = cap;  }
	double get_velocity() { return(velocity); }

	ULL get_current_location() { return(current_location); }
	double get_current_time() { return(current_time); }

	LL get_occupancy() { return(occupancy); }
	LL get_number_of_scheduled_customers() { return(assigned_customers.size()); }
	bool is_idle() { return(idle); }

	std::pair<ULL, double> get_next_route_node() { return(current_route.front()); }	//should always be the same as current position and time (!!!)

	ULL get_number_of_planned_stops() {
		assert(assigned_stops.size() == 2 * assigned_customers.size() - occupancy);
		return(assigned_stops.size());
	};
	double get_planned_time_horizon(double time) {
		if (assigned_stops.empty())
			return(0);
		else
			return(assigned_stops.rbegin()->planned_time - time);
	}
	std::list< stop >::iterator get_next_stop() { return(assigned_stops.begin()); }
	std::list< stop >::iterator no_stop() { return(assigned_stops.end()); }

	double execute_event(double time, traffic_network& n, measurement_collector& m, ULL& total_serviced_requests, bool do_measurement);	//execute event, returns next event time (if any)
	double handle_event_by_type(double time, traffic_network& n, stop& current_stop);
	double new_route(std::deque< std::pair<ULL, double> > param_new_route);

	
	offer best_offer(ULL param_origin, ULL param_destination, double param_request_time, traffic_network& n, offer& current_best_offer, bool calc_p_full);

	double assign_customer(double assignment_time, customer c, offer& o, traffic_network& n);

protected:

private:
	ULL index;

	LL capacity;		//negative value --> infinite capacity
	double velocity;

	//next node to be at and when the transporter arrives
	ULL current_location;
	double current_time;
	ULL node_of_last_stop;
	double time_of_last_stop;


	std::deque< std::pair<ULL, double> > current_route;	//list of nodes on the route to the next stop
	std::list< stop > assigned_stops;

	std::list<customer> assigned_customers;

	LL occupancy;
	bool idle;

	std::mt19937_64& random_generator;
};

#endif // TRANSPORTER_H
