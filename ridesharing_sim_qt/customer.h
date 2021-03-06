#ifndef CUSTOMER_H
#define CUSTOMER_H

#include <cassert>

#ifndef _INTEGER_TYPES
#define ULL uint64_t
#define LL int64_t
#define _INTEGER_TYPES
#endif

#ifndef _EPSILON
#define MACRO_EPSILON 0.000000000001
#define _EPSILON
#endif

#include "traffic_network.h"
#include "transporter.h"

class transporter;
class offer;

class customer
{
public:
	customer(ULL param_origin, ULL param_destination, double param_time, traffic_network& n, transporter& t, offer& o);
	virtual ~customer();

	double calculate_allowed_delay(traffic_network& n, transporter& t);

	//set on construction, cannot be set manually
	ULL get_origin() { return origin; }
	ULL get_destination() { return destination; }
	double get_request_time() { return request_time; }
	double get_allowed_dropoff_delay() { return allowed_dropoff_delay; }
	double get_allowed_pickup_delay() { return allowed_pickup_delay; }
	double get_offer_pickup_time() { return offer_pickup_time; }
	double get_offer_dropoff_time() { return offer_dropoff_time; }

	//set during simulation
	double get_pickup_time() { return pickup_time; }
	void set_pickup_time(double new_pickup_time) { pickup_time = new_pickup_time; }
	double get_dropoff_time() { return dropoff_time; }
	void set_dropoff_time(double new_dropoff_time) { dropoff_time = new_dropoff_time; }

	//modify the allowed delay
	void add_delay(double additional_delay);

protected:

private:
	ULL origin;				//set on construction
	ULL destination;		//set on construction
	double request_time;	//set on construction

	double offer_pickup_time;
	double offer_dropoff_time;

	double pickup_time;
	double dropoff_time;
	double allowed_dropoff_delay;
	double allowed_pickup_delay;
};

#endif // CUSTOMER_H
