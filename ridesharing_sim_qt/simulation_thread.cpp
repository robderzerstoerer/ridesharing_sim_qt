#include "simulation_thread.h"
#include <qmessagebox.h>
#include <ctime>
#include "utility.h"

simulation_thread::simulation_thread(QObject* parent)
	: QThread(parent)
{
}

simulation_thread::~simulation_thread()
{
}


void simulation_thread::run()
{
	par.stop_thread = false;

	// variables for output
	QVector<double> vx, vC;

	for (std::pair<std::string, ULL> topology_n : par.topology_list)
	{
		for (int cap : par.capacity_list)
		{
			std::vector<double> x_list;
			if (par.capacity_list.size() == 1)
				x_list = par.x_list;
			else
			{
				for (double x = 0.5; x + 0.001 < cap; x += 0.5)
				{
					x_list.push_back(x);
				}
			}

			for (double x : x_list)
			{
				if (par.stop_thread)
					break;
				
				simulation_parameters sim_par;
				sim_par.parent = &par;
				sim_par.normalized_request_rate = x;
				sim_par.topology = topology_n.first;
				sim_par.number_of_nodes = topology_n.second;
				sim_par.calc_cap_delay = par.calc_cap_delay;
				sim_par.calc_p_full = par.calc_p_full; // will maybe be changed if Efficiency data has to be calculated (see below)
				sim_par.capacity = cap;

				// full data output
				std::ofstream out;
				std::string out_filename;
				// condensed data for plotting etc
				std::ofstream outplot;
				std::string out_plotfilename;
				// data from pre-calculated efficiency for unlimited capacity
				std::vector<std::pair<ULL, double>> vBEdata;
				std::vector<double> wait_time_data_unlim;

				if (par.capacity_list.size() == 1 && par.x_list.size() == 1)
					out_filename = par.filename;
				else
				{
					out_filename = par.foldername + 
						topology_n.first +
						"_N_" + std::to_string(topology_n.second) +
						"_cap_" + std::to_string(cap) +
						"_x_" + std::to_string(x + 0.0001).substr(0, 3) +
						".dat";
				}

				if (!CUtility::file_exists(out_filename))
					out.open(out_filename.c_str());
				else
				{
					// look for alternative filename
					std::string newname = CUtility::find_new_filename(out_filename);
					out.open(newname.c_str());

					std::string errormessage = "File " + out_filename + " already exists.\nUsing " + newname + " instead.";
					emit ErrorMessage(errormessage.c_str());

					out_filename = newname;
				}

				std::string plotfilename = out_filename + "eff.dat";
				if (!CUtility::file_exists(plotfilename))
					outplot.open(plotfilename.c_str());
				else
				{
					std::string errormessage = "File " + plotfilename + " already exists.";
					emit ErrorMessage(errormessage.c_str());
					return;
				}

				if (par.capacity_list.size() == 1 && par.x_list.size() == 1)
				{
					if (par.calc_p_full)
					{
						// read pre-calculated data for efficiency as a function of B for the present x in this simulation
						std::string eff_filename =
							"eff_" + topology_n.first +
							"_N_" + std::to_string(topology_n.second) +
							"_unl" +
							"_x_" + std::to_string(x + 0.0001).substr(0, 3) +
							".dateff.dat";

						if (!CUtility::read_eff_B_data_from_file(eff_filename, vBEdata))
						{
							std::string errormessage = "Could not open file " + eff_filename;
							emit ErrorMessage(errormessage.c_str());
						}
						else
						{
							std::map<std::string, std::vector<double>> averages, stddevs, counts;
							CUtility::read_file(eff_filename.substr(0, eff_filename.length() - 7), averages, stddevs, counts);
							wait_time_data_unlim = averages["wait_time"];
						}
					}

					simulate_B_list(sim_par, out, outplot, vBEdata, wait_time_data_unlim);
				}
				else
				{
					if (par.simulate_everything)
					{
						// We have to evaluate the efficiencies with unlimited capacities first
						sim_par.capacity = -1;
						sim_par.calc_p_full = false;

						std::ofstream out_eff;
						std::ofstream out_eff_plot;
						// We are simulating for multiple x. Therefore, we need to get a new filename for each x
						std::string eff_filename =
							par.foldername +
							"eff_" + topology_n.first +
							"_N_" + std::to_string(topology_n.second) +
							"_unl" +
							"_x_" + std::to_string(x + 0.0001).substr(0, 3) +
							".dat";

						std::string eff_plotfilename = eff_filename + "eff.dat";

						// if this file doesn't already exist, we have to calculate it
						if (!CUtility::file_exists(eff_plotfilename))
						{
							out_eff.open(eff_filename.c_str());
							out_eff_plot.open(eff_plotfilename.c_str());

							simulate_B_list(sim_par, out_eff, out_eff_plot);

							out_eff.close();
							out_eff_plot.close();
						}

						// read pre-calculated data for efficiency as a function of B for the present x in this simulation
						if (!CUtility::read_eff_B_data_from_file(eff_plotfilename, vBEdata))
						{
							std::string errormessage = "Could not open file " + eff_filename;
							emit ErrorMessage(errormessage.c_str());
						}
					}
					else if (par.calc_p_full)
					{
						// we already have pre-calculated data for efficiency. Read it.
						std::string vBEdata_filename =
							par.foldername +
							"eff_" + topology_n.first +
							"_N_" + std::to_string(topology_n.second) +
							"_unl" +
							"_x_" + std::to_string(x + 0.0001).substr(0, 3) +
							".dateff.dat";

						if (!CUtility::read_eff_B_data_from_file(vBEdata_filename, vBEdata))
						{
							std::string errormessage = "Could not open file " + vBEdata_filename;
							emit ErrorMessage(errormessage.c_str());
						}
					}

					// now the actual calculation with limited capacity. Here, also we need to calculate the p_full(2)
					sim_par.capacity = cap;
					sim_par.calc_p_full = par.calc_p_full;

					simulate_B_list(sim_par, out, outplot, vBEdata, wait_time_data_unlim);

				}

				out.close();
				outplot.close();

				emit ProcessTextChanged("Done simulating");
			}
		}
	}

}

void simulation_thread::simulate_B_list(
	simulation_parameters& sim_par, 
	std::ofstream& out, std::ofstream& outplot, 
	std::vector<std::pair<ULL, double>>& vBEdata,
	std::vector<double> wait_time_data_unlim)
{
	// Qt variables for plotting
	QVector<double> vB, vE;
	
	// first: scan with smaller simulation time
	sim_par.max_sim_time = par.scan_point_sim_time;
	int index_first_successful = 0;
	int index_last_successful = par.B_list.size() - 1;
	bool b_prev_successful = false;
	int current_index = 0;
	for (ULL B : par.B_list)
	{
		if (par.stop_thread)
			return;

		sim_par.number_of_buses = B;

		QString newText = "Simulating " + QString(sim_par.topology.c_str()) +
			", N=" + QString(std::to_string(sim_par.number_of_nodes).c_str()) +
			", x=" + QString(std::to_string(sim_par.normalized_request_rate).c_str()) +
			", cap=" + QString(std::to_string(sim_par.capacity).c_str()) +
			", B=" + QString(std::to_string(sim_par.number_of_buses).c_str());
		emit ProcessTextChanged(newText);
		emit GraphChanged(vB, vE);

		// simulate
		double Eff = single_simulation(sim_par, out, outplot, vBEdata, wait_time_data_unlim);

		if (!b_prev_successful && Eff > 0.0)
		{
			// first successful scan
			index_first_successful = current_index;
			b_prev_successful = true;
		}
		else if (b_prev_successful && Eff < 0.0)
		{
			// last successful scan
			index_last_successful = current_index - 1;
			break;
		}

		if (Eff >= 0.0)
		{
			vB.push_back(sim_par.number_of_buses);
			vE.push_back(Eff);
		}

		current_index++;
	}

	// if none of the scans was successful: try 5th last element with longer sim time
	if (!b_prev_successful && par.B_list.size() >= 5)
	{
		sim_par.max_sim_time = par.max_point_sim_time;

		if (par.stop_thread)
			return;

		sim_par.number_of_buses = par.B_list[par.B_list.size() - 5];

		QString newText = "Simulating " + QString(sim_par.topology.c_str()) +
			", N=" + QString(std::to_string(sim_par.number_of_nodes).c_str()) +
			", x=" + QString(std::to_string(sim_par.normalized_request_rate).c_str()) +
			", cap=" + QString(std::to_string(sim_par.capacity).c_str()) +
			", B=" + QString(std::to_string(sim_par.number_of_buses).c_str());
		emit ProcessTextChanged(newText);
		emit GraphChanged(vB, vE);

		// simulate
		double Eff = single_simulation(sim_par, out, outplot, vBEdata, wait_time_data_unlim);

		if (Eff > 0.0)
		{
			vB.push_front(sim_par.number_of_buses);
			vE.push_front(Eff);
			index_first_successful = par.B_list.size() - 5;
			index_last_successful  = par.B_list.size() - 5;
		}
	}

	// then: if some points have been skipped, try these with a longer simulation time
	sim_par.max_sim_time = par.max_point_sim_time;
	for (int i = index_first_successful - 1; i > 0; i--)
	{
		if (par.stop_thread)
			return;

		sim_par.number_of_buses = par.B_list[i];

		QString newText = "Simulating " + QString(sim_par.topology.c_str()) +
			", N=" + QString(std::to_string(sim_par.number_of_nodes).c_str()) +
			", x=" + QString(std::to_string(sim_par.normalized_request_rate).c_str()) +
			", cap=" + QString(std::to_string(sim_par.capacity).c_str()) +
			", B=" + QString(std::to_string(sim_par.number_of_buses).c_str());
		emit ProcessTextChanged(newText);
		emit GraphChanged(vB, vE);

		// simulate
		double Eff = single_simulation(sim_par, out, outplot, vBEdata, wait_time_data_unlim);

		if (Eff < 0.0)
			break;
		else
		{
			vB.push_front(sim_par.number_of_buses);
			vE.push_front(Eff);
		}
	}

	for (int i = index_last_successful + 1; i < par.B_list.size(); i++)
	{
		if (par.stop_thread)
			return;

		sim_par.number_of_buses = par.B_list[i];

		QString newText = "Simulating " + QString(sim_par.topology.c_str()) +
			", N=" + QString(std::to_string(sim_par.number_of_nodes).c_str()) +
			", x=" + QString(std::to_string(sim_par.normalized_request_rate).c_str()) +
			", cap=" + QString(std::to_string(sim_par.capacity).c_str()) +
			", B=" + QString(std::to_string(sim_par.number_of_buses).c_str());
		emit ProcessTextChanged(newText);
		emit GraphChanged(vB, vE);

		// simulate
		double Eff = single_simulation(sim_par, out, outplot, vBEdata, wait_time_data_unlim);

		if (Eff < 0.0)
			break;
		else
		{
			vB.push_back(sim_par.number_of_buses);
			vE.push_back(Eff);
		}
	}
}

double simulation_thread::single_simulation(
	simulation_parameters& sim_par,
	std::ofstream& out, std::ofstream& outplot,
	std::vector<std::pair<ULL, double>>& vBEdata,
	std::vector<double> wait_time_data_unlim)
{
	ridesharing_sim sim(&sim_par);
	sim.init_network();

	if (!sim.init_new_sim(sim_par.number_of_buses, sim_par.normalized_request_rate, std::max((ULL)50000, par.num_requests_per_bus_init * sim_par.number_of_buses)))
		return -1.0;

	//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
	sim.enable_measurements((2.34546789 * sim_par.number_of_buses) / sim.request_rate);

	if (par.simulate_until_exact)
	{
		bool exact_result = false;
		while (!exact_result)
		{
			if (par.stop_thread)
				break;

			//sim.measurements.scheduled_customers.reset();

			//simulate (and measure) for 1000 requests per bus (at least 50000 requests)
			if (!sim.run_sim_requests(std::max((ULL)50000, par.num_requests_per_bus_sim * sim_par.number_of_buses)))
				return -1.0;

			double total_C_av = sim.measurements.scheduled_customers.get_average();
			double last_partial_C_av = sim.measurements.scheduled_customers.get_partial_averages().back();
			if (abs(last_partial_C_av - total_C_av) < 0.1 * total_C_av)
				exact_result = true;
			else
			{
				QString newText = QString::fromStdString(std::to_string(total_C_av)) + ", " + QString::fromStdString(std::to_string(last_partial_C_av));
				emit ProcessTextChanged(newText);
			}
		}
	}
	else
	{
		//simulate (and measure) for 100 requests per bus (at least 10000 requests)
		if (!sim.run_sim_requests(std::max((ULL)20000, par.num_requests_per_bus_sim * sim_par.number_of_buses)))
			return -1.0;
	}

	//output results
	sim.print_params(out, true);
	sim.print_measurements(out);
	out << std::endl << std::endl;

	double total_C_av = sim.measurements.scheduled_customers.get_average();
	double Eff = sim_par.normalized_request_rate / total_C_av;

	outplot << sim_par.number_of_buses << '\t' << Eff;

	double last_partial_C_av = sim.measurements.scheduled_customers.get_partial_averages().back();
	outplot << '\t' << "x" << '\t' << sim_par.normalized_request_rate / last_partial_C_av;

	if (sim_par.calc_p_full)
	{
		double total_p_full_av = sim.measurements.p_full.get_average();
		outplot << '\t' << "p" << '\t' << total_p_full_av;
		double last_partial_p_full_av = sim.measurements.p_full.get_partial_averages().back();
		outplot << '\t' << "px" << '\t' << last_partial_p_full_av;

		double total_p_full2_av = sim.measurements.p_full2.get_average();
		outplot << '\t' << "p2" << '\t' << total_p_full2_av;
		double last_partial_p_full2_av = sim.measurements.p_full2.get_partial_averages().back();
		outplot << '\t' << "p2x" << '\t' << last_partial_p_full2_av;

		double total_p_full_snap_av = sim.measurements.p_full_snap.get_average();
		outplot << '\t' << "p_snap" << '\t' << total_p_full_snap_av;

		/* now:
			1. Make a decision which calculated efficiency is the correct one (the average or the last value)
			2. Make a decision which calculated p_full(2) is the correct one (the average or the last value)
			3. Find the number of buses that would have the same efficiency in a simulation with unlimited capacity
			4. Find the effective number of buses calculated by B*(1-p_full) and B*(1-p_full2)
			and output them in this order.
		*/
		double real_eff;
		double real_p_full;
		double real_p_full2;

		/*
		if (abs(last_partial_C_av - total_C_av) < 0.06 * total_C_av)
			real_eff = par.normalized_request_rate / total_C_av;
		else
			real_eff = par.normalized_request_rate / last_partial_C_av;

		if (abs(last_partial_p_full_av - total_p_full_av) < 0.06 * total_p_full_av)
			real_p_full = total_p_full_av;
		else
			real_p_full = last_partial_p_full_av;

		if (abs(last_partial_p_full2_av - total_p_full2_av) < 0.06 * total_p_full2_av)
			real_p_full2 = total_p_full2_av;
		else
			real_p_full2 = last_partial_p_full2_av;
			*/

		real_eff = sim_par.normalized_request_rate / total_C_av;
		real_p_full = total_p_full_av;
		real_p_full2 = total_p_full2_av;

		double B_equiv = 0.0;
		if (vBEdata.size() > 0) // make sure we have the data
			B_equiv = CUtility::find_effective_B(vBEdata, real_eff);
		double B_eff_p_full = sim_par.number_of_buses * (1 - real_p_full);
		double B_eff_p_full2 = sim_par.number_of_buses * (1 - real_p_full2);

		outplot << '\t' << B_equiv << '\t' << B_eff_p_full << '\t' << B_eff_p_full2; 

		/*
		int found_B = -1;
		for (int f = 0; f < par.B_list.size(); f++)
		{
			if (par.B_list[f] == sim_par.number_of_buses)
			{
				found_B = f;
				break;
			}
		}
		if (found_B != -1)
		{
			double C_unlim = total_C_av - wait_time_data_unlim[found_B] * sim.network.get_mean_dropoff_distance() / sim_par.normalized_request_rate * CUtility::polylog_B(total_p_full2_av, sim_par.number_of_buses);
			outplot << '\t' << CUtility::find_effective_B(vBEdata, sim_par.normalized_request_rate / C_unlim);
		}
		else
		{

		} */
	}

	outplot << '\n';

	return Eff;
}