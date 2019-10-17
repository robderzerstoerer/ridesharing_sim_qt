#include "simulation_thread.h"
#include <qmessagebox.h>
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
	std::ofstream out;
	if (!CUtility::file_exists(par.filename))
		out.open(par.filename.c_str());
	else
	{
		std::string errormessage = "File " + par.filename + " already exists.";
		emit ErrorMessage(errormessage.c_str());
		return;
	}

	//initialize simulation class
	ridesharing_sim sim(&par);
	sim.init_network();

	if (par.number_of_bus_calculations > 1)
	{
		QVector<double> vB, vE;
		QVector<double> vx, vC, vB_equiv_x_inf, vB_equiv_x;
		std::vector<std::pair<ULL, double>> vBEdata;

		// read pre-calculated data for efficiency as a function of B for the present x in this simulation
		std::string eff_filename = "eff_torus_100_unl_x_" + std::to_string(par.normalized_request_rate).substr(0, 3) + ".dateff.dat";
		
		if (!CUtility::read_eff_B_data_from_file(eff_filename, vBEdata))
		{
			std::string errormessage = "Could not open file " + eff_filename;
			emit ErrorMessage(errormessage.c_str());
		}

		std::ofstream outplot;
		std::string plotfilename = par.filename + "eff.dat";
		if (!CUtility::file_exists(plotfilename))
			outplot.open(plotfilename.c_str());
		else
		{
			std::string errormessage = "File " + plotfilename + " already exists.";
			emit ErrorMessage(errormessage.c_str());
			return;
		}

		emit ProcessTextChanged("Starting simulation");
		emit GraphChanged(vB, vE);
		//ui.textEdit_6->setPlainText(ui.textEdit_6->toPlainText() + "\nStarting simulation\n");

		double index_low_B = log(par.number_of_buses_from);
		double index_high_B = log(par.number_of_buses_to);
		for (double iter_B = index_low_B; iter_B <= (index_high_B + 0.0001); iter_B += (index_high_B - index_low_B) / (par.number_of_bus_calculations - 1))
		{
			if (par.stop_thread)
				break;
			ULL number_of_buses = (int)(exp(iter_B) + 0.5);
			QString newText = "Simulating " + QString::fromStdString(std::to_string(number_of_buses)) + " buses\n";
			emit ProcessTextChanged(newText);
			emit GraphChanged(vB, vE);
			
			if (par.number_of_request_rates > 1)
			{
				double index_low_x = log(par.normalized_request_rate_from);
				double index_high_x = log(par.normalized_request_rate_to);
				for (double iter_x = index_low_x; iter_x <= (index_high_x + 0.0001); iter_x += (index_high_x - index_low_x) / (par.number_of_request_rates - 1))
				{
					if (par.stop_thread)
						break;
					double normalized_request_rate = exp(iter_x);
					QString newText = "Simulating B = " + QString::fromStdString(std::to_string(number_of_buses)) + ", x = " + QString::fromStdString(std::to_string(normalized_request_rate));
					emit ProcessTextChanged(newText);
					emit GraphChanged(vx, vC);
					//ui.textEdit_6->setPlainText(newText.c_str());

					sim.init_new_sim(number_of_buses, normalized_request_rate, 500 * number_of_buses);

					//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
					sim.enable_measurements((1.0 * number_of_buses) / sim.request_rate);
					//simulate (and measure) for 100 requests per bus (at least 10000 requests)
					sim.run_sim_requests(std::max((ULL)10000, 500 * number_of_buses));

					//output results
					if (par.save)
					{
						sim.print_params(out, true);
						sim.print_measurements(out);
						out << std::endl << std::endl;
					}

					vB.push_back((double) number_of_buses);
					double Eff = normalized_request_rate / sim.measurements.scheduled_customers.get_average();
					vE.push_back(Eff);
					vx.push_back((double)normalized_request_rate);
					vC.push_back(sim.measurements.scheduled_customers.get_average());
					
					// find equivalent number of buses in unlimited-capacity-simulation
					vB_equiv_x_inf.push_back(CUtility::find_effective_B(vBEdata, Eff));

					outplot << vB.back() << '\t' << vE.back() << '\t' << vx.back() << '\t' << vC.back() << '\t' << vB_equiv_x_inf.back() << '\n';
				}
			}
			else
			{
				sim.init_new_sim(number_of_buses, par.normalized_request_rate, std::max((ULL)50000, 100 * number_of_buses));
				
				//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
				sim.enable_measurements((2.0 * number_of_buses) / sim.request_rate);
				
				if (par.simulate_until_exact)
				{
					bool exact_result = false;
					while (!exact_result)
					{
						if (par.stop_thread)
							break;

						sim.measurements.scheduled_customers.reset();

						//simulate (and measure) for 100 requests per bus (at least 10000 requests)
						sim.run_sim_requests(std::max((ULL)50000, 500 * number_of_buses));

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
					sim.run_sim_requests(std::max((ULL)10000, 500 * number_of_buses));
				}

				//output results
				if (par.save)
				{
					sim.print_params(out, true);
					sim.print_measurements(out);
					out << std::endl << std::endl;
				}

				vB.push_back((double)number_of_buses);
				double total_C_av = sim.measurements.scheduled_customers.get_average();
				double Eff = par.normalized_request_rate / total_C_av;
				vE.push_back(Eff);

				outplot << vB.back() << '\t' << vE.back();

				double last_partial_C_av = sim.measurements.scheduled_customers.get_partial_averages().back();
				outplot << '\t' << "x" << '\t' << par.normalized_request_rate / last_partial_C_av;

				if (par.calc_p_full)
				{
					double total_p_full_av = sim.measurements.p_full.get_average();
					outplot << '\t' << "p" << '\t' << total_p_full_av;
					double last_partial_p_full_av = sim.measurements.p_full.get_partial_averages().back();
					outplot << '\t' << "px" << '\t' << last_partial_p_full_av;

					double total_p_full2_av = sim.measurements.p_full2.get_average();
					outplot << '\t' << "p2" << '\t' << total_p_full2_av;
					double last_partial_p_full2_av = sim.measurements.p_full2.get_partial_averages().back();
					outplot << '\t' << "p2x" << '\t' << last_partial_p_full2_av;

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

					real_eff = par.normalized_request_rate / total_C_av;
					real_p_full = total_p_full_av;
					real_p_full2 = total_p_full2_av;

					double B_equiv = CUtility::find_effective_B(vBEdata, real_eff);
					double B_eff_p_full = number_of_buses * (1 - real_p_full);
					double B_eff_p_full2 = number_of_buses * (1 - real_p_full2);

					outplot << '\t' << B_equiv << '\t' << B_eff_p_full << '\t' << B_eff_p_full2;
				}

				outplot << '\n';
			}
		}

		//plt::loglog(vB.toStdVector(), vE.toStdVector());
		//plt::show();

		outplot.close();

		emit ProcessTextChanged("Done simulating");

	}
	else if (par.number_of_request_rates > 1)
	{
		QVector<double> vx, vC;

		emit ProcessTextChanged("Starting simulation");
		emit GraphChanged(vx, vC);
		//ui.textEdit_6->setPlainText(ui.textEdit_6->toPlainText() + "\nStarting simulation\n");

		std::ofstream outplot;
		std::string plotfilename = par.filename + "plot.dat";
		if (!CUtility::file_exists(plotfilename))
			outplot.open(plotfilename.c_str());
		else
		{
			std::string errormessage = "File " + plotfilename + " already exists.";
			emit ErrorMessage(errormessage.c_str());
			return;
		}

		double index_low = log(par.normalized_request_rate_from);
		double index_high = log(par.normalized_request_rate_to);
		for (double iter = index_low; iter <= (index_high + 0.0001); iter += (index_high - index_low) / (par.number_of_request_rates - 1))
		{
			if (par.stop_thread)
				break;
			double normalized_request_rate = exp(iter);
			QString newText = "Simulating x = " + QString::fromStdString(std::to_string(normalized_request_rate));
			emit ProcessTextChanged(newText);
			emit GraphChanged(vx, vC);
			//ui.textEdit_6->setPlainText(newText.c_str());

			sim.init_new_sim(par.number_of_buses, normalized_request_rate, std::max((ULL)30000, 1000 * par.number_of_buses));

			//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
			sim.enable_measurements((1.0 * par.number_of_buses) / sim.request_rate);

			if (par.simulate_until_exact)
			{
				bool exact_result = false;
				while (!exact_result)
				{
					if (par.stop_thread)
						break;

					sim.measurements.scheduled_customers.reset();

					//simulate (and measure) for 100 requests per bus (at least 10000 requests)
					sim.run_sim_requests(std::max((ULL)50000, 1000 * par.number_of_buses));

					double total_C_av = sim.measurements.scheduled_customers.get_average();
					double last_partial_C_av = sim.measurements.scheduled_customers.get_partial_averages().back();
					if (abs(last_partial_C_av - total_C_av) < 0.1 * total_C_av)
						exact_result = true;
				}
			}
			else
			{
				//simulate (and measure) for 100 requests per bus (at least 10000 requests)
				sim.run_sim_requests(std::max((ULL)50000, 1000 * par.number_of_buses));
			}

			//output results
			if (par.save)
			{
				sim.print_params(out, true);
				sim.print_measurements(out);
				out << std::endl << std::endl;
			}

			vx.push_back((double)normalized_request_rate);
			vC.push_back(sim.measurements.scheduled_customers.get_average());

			outplot << vx.back() << '\t' 
				<< vC.back() << '\t' 
				<< sim.measurements.scheduled_customers.get_stddev() << '\t'
				<< CUtility::two_node_av_scheduled_customers(5, normalized_request_rate, 2) << '\t'
				<< CUtility::two_node_stddev_scheduled_customers(5, normalized_request_rate, 2);

			if (par.calc_p_full)
			{
				double total_p_full_av = sim.measurements.p_full.get_average();
				outplot << '\t' << "p1" << '\t' << total_p_full_av;
				double last_partial_p_full_av = sim.measurements.p_full.get_partial_averages().back();
				outplot << '\t' << "px" << '\t' << last_partial_p_full_av;

				double total_p_full2_av = sim.measurements.p_full2.get_average();
				outplot << '\t' << "p2" << '\t' << total_p_full2_av;
				double last_partial_p_full2_av = sim.measurements.p_full2.get_partial_averages().back();
				outplot << '\t' << "p2x" << '\t' << last_partial_p_full2_av;
			}

			double last_partial_C_av = sim.measurements.scheduled_customers.get_partial_averages().back();
			

			outplot << '\t' << "x" << '\t' << last_partial_C_av;

			

			outplot << '\n';
		}

		//plt::loglog(vB.toStdVector(), vE.toStdVector());
		//plt::show();

		outplot.close();

		emit ProcessTextChanged("Done simulating");
	}
	else
	{
		emit ProcessTextChanged("Starting simulation: ");

		sim.init_new_sim(par.number_of_buses, par.normalized_request_rate, 10 * par.number_of_buses);

		//turn on measurements with a given step size, measure every (number of buses) requests for a total of ~ 100 measurements
		sim.enable_measurements((1.0 * par.number_of_buses) / sim.request_rate);
		//simulate (and measure) for 100 requests per bus (at least 10000 requests)
		sim.run_sim_requests(std::max((ULL)10000, 100 * par.number_of_buses));

		emit ProcessTextChanged("Done simulating!");

		//output results
		sim.print_params(out, true);
		sim.print_measurements(out);

		out << std::endl << std::endl;
	}

	out.close();
}