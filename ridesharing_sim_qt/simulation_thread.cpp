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


void sort_B_E_vectors(std::vector<int>& vB, std::vector<double>& vE)
{
	// selection sort
	for (int i = 0; i < vB.size(); i++)
	{
		int smallest_B = 1000000;
		int smallest_B_index = -1;

		for (int j = i; j < vB.size(); j++)
		{
			if (vB[j] < smallest_B)
			{
				smallest_B = vB[j];
				smallest_B_index = j;
			}
		}

		// swap
		vB[smallest_B_index] = vB[i];
		vB[i] = smallest_B;
		double smallest_E = vE[smallest_B_index];
		vE[smallest_B_index] = vE[i];
		vE[i] = smallest_E;
	}
}

void simulation_thread::run()
{
	par.stop_thread = false;

	// Qt variables for output
	QVector<double> vx, vC;

	std::map<std::string, std::vector<std::pair<float, float>>> bhalves;
	if (par.sim_mode == B_AIM)
	{
		if (!CUtility::read_B_halves("bhalves.txt", bhalves))
		{
			std::string errormessage = "Could not open file bhalves.txt";
			emit ErrorMessage(errormessage.c_str());
		}
	}

	for (std::pair<std::string, ULL> topology_n : par.topology_list)
	{
		if (par.sim_mode == B_AIM)
		{
			std::string topokey = topology_n.first + "_N_" + std::to_string(topology_n.second);
			if (bhalves.find(topokey) == bhalves.end())
			{
				std::string errormessage = "Topology " + topokey + " not calculated";
				emit ErrorMessage(errormessage.c_str());
				continue;
			}

			double bhalf = CUtility::find_B_half_x(bhalves[topokey], 7.5);
			double x_guess = par.E_aim / (1 - par.E_aim) + FIT_A * (FIT_A / 2 - sqrt(par.E_aim / (1 - par.E_aim) + 0.25));
			ULL B_guess = MAX(int(x_guess * bhalf + 0.5), 1);
			ULL last_B_guess = B_guess;
			double last_Eff = -1.0;
			std::vector<int> B_theta(par.capacity_list[0] + 1, 0);
			std::vector<double> E_theta(par.capacity_list[0] + 1, 0);

			traffic_network n;
			n.init(topology_n.first, topology_n.second);
			double lav = n.get_mean_dropoff_distance();

			for (int cap : par.capacity_list)
			{
				simulation_parameters sim_par;
				sim_par.parent = &par;
				sim_par.topology = topology_n.first;
				sim_par.number_of_nodes = topology_n.second;
				sim_par.calc_cap_delay = false;
				sim_par.calc_p_full = false;
				sim_par.capacity = cap;
				sim_par.max_sim_time = par.max_point_sim_time;

				// output dummies
				std::ofstream out;
				std::ofstream outplot;

				bool finished = false;
				int last_res = 0;
				int step_size = MAX(int(B_guess * 0.1), 1);
				std::vector<int> checked_B_list;
				std::vector<double> checked_E_list;
				int min_B = 0;
				int max_B = 10000000;

				while (!finished)
				{
					if (par.stop_thread)
						break;
					double x = par.lambda * lav / (1.0 * B_guess);  // replace 1.0 with transporter velocity if neccessary
					if (x >= cap)
					{
						min_B = MAX(min_B, B_guess);
						if (last_res != 0 && step_size > 1)
							step_size--;
						B_guess = B_guess + step_size;
						continue;
					}

					std::vector<int>::iterator it = std::find(checked_B_list.begin(), checked_B_list.end(), B_guess);
					if (it != checked_B_list.end())
					{
						if (step_size == 1)
						{
							// done, now find out which B gives closest efficiency to aim
							finished = true;
							sort_B_E_vectors(checked_B_list, checked_E_list);
							int closest_B = -1;
							double closest_E = -1.0;

							for (int i = 0; i < checked_E_list.size(); i++)
							{
								if (abs(checked_E_list[i] - par.E_aim) < abs(closest_E - par.E_aim))
								{
									closest_B = checked_B_list[i];
									closest_E = checked_E_list[i];
								}
							}

							// write solutions
							B_theta[cap] = closest_B;
							E_theta[cap] = closest_E;
							break;
						}
						else
						{
							step_size = (int) (step_size / 2.0);
							B_guess = (int) ((last_B_guess + B_guess) / 2 + 0.5);
							continue;
						}
					}
					else
					{
						if (B_guess <= min_B || B_guess >= max_B)
						{
							step_size = (int)(step_size / 2.0);
							B_guess = (int)((last_B_guess + B_guess) / 2 + 0.5);
							continue;
						}

						// not yet done iterating
						double B_half_x = CUtility::find_B_half_x(bhalves[topokey], x);

						sim_par.normalized_request_rate = x;
						sim_par.number_of_buses = B_guess;

						QString newText = "Simulating " + QString(sim_par.topology.c_str()) +
							", N=" + QString(std::to_string(sim_par.number_of_nodes).c_str()) +
							", x=" + QString(std::to_string(sim_par.normalized_request_rate).c_str()) +
							", cap=" + QString(std::to_string(sim_par.capacity).c_str()) +
							", B=" + QString(std::to_string(sim_par.number_of_buses).c_str());
						emit ProcessTextChanged(newText);

						double Eff = single_simulation(sim_par, out, outplot);

						checked_B_list.push_back(B_guess);
						checked_E_list.push_back(Eff);

						if (last_res == 0) //no previous calculation
						{
							if (Eff < par.E_aim)
							{
								min_B = B_guess;
								last_B_guess = B_guess;
								last_res = -1;
								last_Eff = Eff;
								B_guess = B_guess + step_size;
							}
							else
							{
								max_B = B_guess;
								if (B_guess == 1)
								{
									B_theta[cap] = B_guess;
									E_theta[cap] = Eff;
									break;
								}
								last_B_guess = B_guess;
								last_res = 1;
								last_Eff = Eff;
								B_guess = B_guess - step_size;
							}
						}
						else if (last_res == -1) // previous calculation gave too low efficiency
						{
							if (Eff < par.E_aim)
							{
								min_B = MAX(min_B, B_guess);
								last_B_guess = B_guess;
								last_res = -1;
								last_Eff = Eff;
								B_guess = B_guess + step_size;
							}
							else
							{
								max_B = MIN(max_B, B_guess);
								last_B_guess = B_guess;
								last_res = 1;
								last_Eff = Eff;
								step_size = MAX(int(step_size / 2.0 + 0.5), 1);
								B_guess = B_guess - step_size;
							}
						}
						else if (last_res == 1) // previous calculation gave too high efficiency
						{
							if (Eff > par.E_aim)
							{
								max_B = MIN(max_B, B_guess);
								if (B_guess == 1)
								{
									B_theta[cap] = B_guess;
									E_theta[cap] = Eff;
									break;
								}
								last_B_guess = B_guess;
								last_res = 1;
								last_Eff = Eff;
								B_guess = B_guess - step_size;
							}
							else
							{
								min_B = MAX(min_B, B_guess);
								last_B_guess = B_guess;
								last_res = 1;
								last_Eff = Eff;
								step_size = MAX(int(round(step_size / 2)), 1);
								B_guess = B_guess + step_size;
							}
						}
					}
				}
			}

			// output
			std::ofstream out;
			std::string out_filename = topology_n.first +
				"_N_" + std::to_string(topology_n.second) +
				"_aim_" + std::to_string(par.E_aim + 0.0001).substr(0, 5) +
				"_lambda_" + std::to_string(par.lambda + 0.0001).substr(0, 5) +
				".dat";

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

			std::string msg = "";
			for (int i = 1; i < B_theta.size(); i++)
			{
				out << i << '\t' << B_theta[i] << '\t' << E_theta[i] << '\n';
				msg += std::to_string(i) + " " + std::to_string(B_theta[i]) + " " + std::to_string(E_theta[i]) + '\n';
			}
			emit ErrorMessage(msg.c_str());

			out.close();

			emit ProcessTextChanged("Done simulating");
		}
		else
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
						}

						simulate_B_list(sim_par, out, outplot, vBEdata);
					}
					else
					{
						if (par.sim_mode == EVERYTHING)
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

						simulate_B_list(sim_par, out, outplot, vBEdata);

					}

					out.close();
					outplot.close();

					emit ProcessTextChanged("Done simulating");
				}
			}
		}
	}

}

void simulation_thread::simulate_B_list(
	simulation_parameters& sim_par, 
	std::ofstream& out, std::ofstream& outplot, 
	std::vector<std::pair<ULL, double>>& vBEdata)
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
		double Eff = single_simulation(sim_par, out, outplot, vBEdata);

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
		double Eff = single_simulation(sim_par, out, outplot, vBEdata);

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
		double Eff = single_simulation(sim_par, out, outplot, vBEdata);

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
		double Eff = single_simulation(sim_par, out, outplot, vBEdata);

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
	std::vector<std::pair<ULL, double>>& vBEdata)
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
	if (out.is_open())
	{
		sim.print_params(out, true);
		sim.print_measurements(out);
		out << std::endl << std::endl;
	}

	double total_C_av = sim.measurements.scheduled_customers.get_average();
	double Eff = sim_par.normalized_request_rate / total_C_av;

	if (outplot.is_open())
	{
		outplot << sim_par.number_of_buses << '\t' << Eff;

		double last_partial_C_av = sim.measurements.scheduled_customers.get_partial_averages().back();
		outplot << '\t' << "x" << '\t' << sim_par.normalized_request_rate / last_partial_C_av;
	}

	if (sim_par.calc_p_full && outplot.is_open())
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

	}

	if (outplot.is_open())
		outplot << '\n';

	return Eff;
}