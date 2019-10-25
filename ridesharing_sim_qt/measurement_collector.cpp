#include "measurement_collector.h"

measure::measure()
{
	// Variables for measuring the number of scheduled customers in equidistant time intervals
	last_iTime = 0;
	double cur_average = 0.0;
	int cur_n = 0;
}

void measure::new_measurement(long double val, double time)
{
	long double delta1, delta2;

	n += 1;
	delta1 = val - av;
	av += delta1 / n;
	delta2 = val - av;
	stddev += delta1 * delta2;

	if (calc_partial_averages && time > 0)
	{
		int iTime = trunc(time / partial_average_time_interval);
		if (iTime > last_iTime)
		{
			partial_averages.push_back(cur_average);
			cur_n = 0;
			cur_average = 0.0;
			last_iTime = iTime;
		}

		cur_n += 1;
		double delta = val - cur_average;
		cur_average += delta / cur_n;
	}
}

void measure::reset()
{
	n = 0;
	av = 0;
	stddev = 0;

	// reset help variables
	last_iTime = 0;
	double cur_average = 0.0;
	int cur_n = 0;
}

void measure::enable_partial_average_calculation(double time_interval)
{
	calc_partial_averages = true; 
	partial_average_time_interval = time_interval;

	// reset help variables
	last_iTime = 0;
	double cur_average = 0.0;
	int cur_n = 0;
}

void measure::print(std::ofstream& out)
{
	if (n > 0)
		out << n << '\t' << av << '\t' << sqrt(stddev / n) << '\t';
	else
		out << 0 << '\t' << 0 << '\t' << 0 << '\t';
};

measurement_collector::measurement_collector(traffic_network& param_network) : network(param_network)
{

}

measurement_collector::~measurement_collector()
{
	//dtor
}

void measurement_collector::measure_request(customer& c, transporter& t)
{
	wait_time.new_measurement(c.get_pickup_time() - c.get_request_time());
	drive_time.new_measurement(c.get_dropoff_time() - c.get_pickup_time());
	delay_time.new_measurement(c.get_dropoff_time() - c.get_pickup_time() - network.get_network_distance(c.get_origin(), c.get_destination()) / t.get_velocity());

	if (abs(c.get_dropoff_time() - c.get_pickup_time() - network.get_network_distance(c.get_origin(), c.get_destination()) / t.get_velocity()) > 10 * MACRO_EPSILON)
		fraction_of_delayed_trips.new_measurement(1);
	else
		fraction_of_delayed_trips.new_measurement(0);
}

void measurement_collector::measure_trip(std::pair<ULL, double> last_stop, std::pair<ULL, double> current_stop)
{
	drive_distance_between_stops.new_measurement(network.get_network_distance(last_stop.first, current_stop.first));
	drive_time_between_stops.new_measurement(current_stop.second - last_stop.second);
}

//void measurement_collector::measure_stops(...)
//{
//	new_measurement(drive_time_between_stops, c.get_pickup_time() - c.get_request_time());
//	new_measurement(drive_distance_between_stops, c.get_pickup_time() - c.get_request_time());
//}

void measurement_collector::measure_system_status(std::vector<transporter>& transporter_list, double time)
{
	ULL idle_transporters = 0;
	for (transporter& t : transporter_list)
	{
		occupancy.new_measurement(t.get_occupancy());
		scheduled_customers.new_measurement(t.get_number_of_scheduled_customers(), time);
		number_of_planned_stops.new_measurement(t.get_number_of_planned_stops());
		planned_time_horizon.new_measurement(t.get_planned_time_horizon(time));
		if (t.is_idle())
			++idle_transporters;
		if (t.get_occupancy() == t.get_capacity())
			p_full_snap.new_measurement(1);
		else
			p_full_snap.new_measurement(0);
	}
	number_of_idle_transporters.new_measurement(idle_transporters);
}

void measurement_collector::reset()
{
	wait_time.reset();
	drive_time.reset();
	delay_time.reset();
	fraction_of_delayed_trips.reset();

	drive_time_between_stops.reset();
	drive_distance_between_stops.reset();

	occupancy.reset();
	scheduled_customers.reset();
	number_of_planned_stops.reset();
	planned_time_horizon.reset();
	number_of_idle_transporters.reset();

	p_full.reset();
	p_full2.reset();
	p_full_snap.reset();

	// important: scheduled customers have to be measured every 100 time units
	// (this is used for validating the total average result)
	scheduled_customers.enable_partial_average_calculation(100.0);

	p_full.enable_partial_average_calculation(100.0);
	p_full2.enable_partial_average_calculation(100.0);
}

void measurement_collector::print(std::ofstream& out, bool readable)
{
	if (!readable)
	{
		wait_time.print(out);
		drive_time.print(out);
		delay_time.print(out);
		fraction_of_delayed_trips.print(out);

		drive_time_between_stops.print(out);
		drive_distance_between_stops.print(out);

		occupancy.print(out);
		scheduled_customers.print(out);
		number_of_planned_stops.print(out);
		planned_time_horizon.print(out);
		number_of_idle_transporters.print(out);

		p_full.print(out);
		p_full2.print(out);
		p_full_snap.print(out);

		out << std::endl;
	}
	else
	{
		out << "MEASUREMENTS" << std::endl << std::endl;

		out << "wait_time" << std::endl;
		wait_time.print(out);
		out << std::endl;

		out << "drive_time" << std::endl;
		drive_time.print(out);
		out << std::endl;

		out << "delay_time" << std::endl;
		delay_time.print(out);
		out << std::endl;

		out << "fraction_of_delayed_trips" << std::endl;
		fraction_of_delayed_trips.print(out);
		out << std::endl;

		out << "drive_time_between_stops" << std::endl;
		drive_time_between_stops.print(out);
		out << std::endl;

		out << "drive_distance_between_stops" << std::endl;
		drive_distance_between_stops.print(out);
		out << std::endl;

		out << "occupancy" << std::endl;
		occupancy.print(out);
		out << std::endl;

		out << "scheduled_customers" << std::endl;
		scheduled_customers.print(out);
		out << std::endl;

		out << "number_of_planned_stops" << std::endl;
		number_of_planned_stops.print(out);
		out << std::endl;

		out << "planned_time_horizon" << std::endl;
		planned_time_horizon.print(out);
		out << std::endl;

		out << "number_of_idle_transporters" << std::endl;
		number_of_idle_transporters.print(out);
		out << std::endl;

		out << "p_full" << std::endl;
		p_full.print(out);
		out << std::endl;

		out << "p_full2" << std::endl;
		p_full2.print(out);
		out << std::endl;

		out << "p_full_snap" << std::endl;
		p_full_snap.print(out);
		out << std::endl;

		out << std::endl;
	}
}